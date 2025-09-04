function run_casper_simulation
% RUN_CASPER Run PSA simulation with multiple cycles and progress display.
%
tic;
[sim, bed_states, tank_states] = createPSASimulation();
for b = 1:sim.num_beds
    sim.diagnostics(b) = struct( ...
        'time',[], 'Pin',[], 'Pout',[], 'y1_in',[], 'y1_out',[], 'u_in',[], 'u_out',[]);
end

num_cycles=sim.n_cycles; 

% Build initial state  Y0
Y0 = [];
for b = 1:sim.num_beds
    bed = bed_states{b};
    Y0 = [Y0; bed.C; reshape(bed.C_species, [],1); bed.T; reshape(bed.q, [],1)];
end

opts = odeset('RelTol', sim.RelTol, ...
              'AbsTol', sim.AbsTol, ...
              'MaxStep', sim.dt_max);
% opts = odeset('RelTol', sim.RelTol, ...
%               'AbsTol', sim.AbsTol, ...
%               'Stats', 'on', ...
%               'MaxStep', sim.dt_max, ...
%               'OutputFcn', @odeProgress);
% 
% ---- NEW: cycle accumulators (beds & tanks) ----
bed_cycle = init_bed_accum(sim);
tank_cycle = init_tank_accum(sim);
T_all = [];
Y_all = [];
bed_states_all = cell(num_cycles, 1);  % Store final bed states per cycle
tank_states_all = cell(num_cycles, 1); % (optional) for tanks too

for cycle = 1:num_cycles
    fprintf('\n=== Beginning cycle %d of %d ===\n', cycle, num_cycles);
    % ---- reset per-cycle accumulators each cycle ----
    bed_cycle = init_bed_accum(sim);
    tank_cycle = init_tank_accum(sim);

    for idx = 1:length(sim.step_times)-1
        t0 = sim.step_times(idx);
        tf = sim.step_times(idx+1);
        current_step = sim.step(idx);
        
        fprintf('  Cycle %d — Step %d: t = %.2f to %.2f s\n', cycle, idx, t0, tf);
        
        % Integration across this step
        [Tseg, Yseg] = ode15s(@(t,Y) multi_bed_ode(t, Y, sim, current_step), ...
                              [t0, tf], Y0, opts);
        
        % T_all = [T_all; Tseg];
        % Y_all = [Y_all; Yseg];
        Y0 = Yseg(end,:)';  % Update initial state for next step
        % Boundary-only logging (cheap)
        beds = unpack_bed_state_vector(Y0, sim);
        for b = 1:sim.num_beds
            P = beds{b}.C .* sim.R .* beds{b}.T;
            sim.diagnostics(b).time(end+1)  = tf;
            sim.diagnostics(b).Pin(end+1)   = P(1);
            sim.diagnostics(b).Pout(end+1)  = P(end);
            sim.diagnostics(b).y1_in(end+1) = beds{b}.y(1,1);
            sim.diagnostics(b).y1_out(end+1)= beds{b}.y(end,1);
        end
            % ---- NEW: integrate signed boundary flows over this step & print ----
            step_flow = integrate_step_flows_signed(sim, current_step, Tseg, Yseg);
            print_step_flow_report_signed(sim, idx, step_flow);  % neat, aligned
            % accumulate to cycle totals
            bed_cycle = accum_bed_signed(bed_cycle, step_flow.beds);
            tank_cycle = accum_tank(tank_cycle, step_flow.tanks);
            plot_diagnostics(sim, Y0);    % keep this lightweight
            drawnow limitrate
        end

    % ---- End of cycle: print cycle totals (beds + tanks) ----
    print_cycle_totals(sim, cycle, bed_cycle, tank_cycle);

    end
    bed_states_all{cycle} = unpack_bed_state_vector(Y0, sim);
    tank_states_all{cycle} = tank_states;  % If tanks are updated per step
toc;


save('cycle_sim.mat', 'sim', 'bed_states_all', 'tank_states_all');

fprintf('\n✅ Simulation complete. Final state saved to cycle_sim.mat\n');


% Example output: final concentration profile of Bed A
Ct_end = Y0(1:sim.num_nodes);
figure; plot(sim.z_nodes, Ct_end, 'LineWidth', 2);
xlabel('z (m)'); ylabel('Total Concentration (mol/m³)');
title('Final Bed A Concentration Profile'); grid on;
end


function status = combinedOutputFcn(t, y, flag)
    status1 = odeProgress(t, y, flag);
    status2 = odeplot(t, y, flag);
    status = status1;
    status = status1 || status2;
end



function S = init_bed_accum(sim)
% Per-bed accumulators for signed boundary totals at z0 and zL
S = struct();
for b = 1:sim.num_beds
    S(b).Fz0_total = 0;                    % mol (signed, +z positive)
    S(b).FzL_total = 0;                    % mol (signed, +z positive)
    S(b).Fz0_species = zeros(sim.n_species,1);
    S(b).FzL_species = zeros(sim.n_species,1);
end
end

function S = init_tank_accum(sim)
% Per-tank totals (positive into tank, negative out of tank)
S = struct();
for k = 1:numel(sim.tanks)
    S(k).name = sim.tanks(k).name;
    S(k).n_total  = 0;                     % mol (signed; + into tank)
    S(k).n_species = zeros(sim.n_species,1);
end
end

function Acc = accum_bed_signed(Acc, StepBeds)
for b = 1:numel(Acc)
    Acc(b).Fz0_total   = Acc(b).Fz0_total   + StepBeds(b).Fz0_total;
    Acc(b).FzL_total   = Acc(b).FzL_total   + StepBeds(b).FzL_total;
    Acc(b).Fz0_species = Acc(b).Fz0_species + StepBeds(b).Fz0_species;
    Acc(b).FzL_species = Acc(b).FzL_species + StepBeds(b).FzL_species;
end
end

function Acc = accum_tank(Acc, StepTanks)
for k = 1:numel(Acc)
    idx = find(strcmpi({StepTanks.name}, Acc(k).name), 1);
    if ~isempty(idx)
        Acc(k).n_total   = Acc(k).n_total   + StepTanks(idx).n_total;
        Acc(k).n_species = Acc(k).n_species + StepTanks(idx).n_species;
    end
end
end

function step_flow = integrate_step_flows_signed(sim, current_step, Tseg, Yseg)
% Integrate signed molar flows at z0 and zL for each bed over the step.
% Also integrate tank totals (positive into tank, negative out).

    nB = sim.num_beds; nS = sim.n_species;
    beds_acc  = init_bed_accum(sim);
    tanks_acc = init_tank_accum(sim);

    % Previous samples for trapezoid
    prev = struct();
    prev.has = false;

    for k = 1:numel(Tseg)
        t  = Tseg(k);
        Yk = Yseg(k,:).';

        % Interface bed properties (P,T,y at each end)
        raw = cell(nB,1);
        for b = 1:nB
            raw{b} = get_bed_interface_props(Yk, b, sim);
        end
        % Boundary conditions / velocities
        [bc0, bcL] = flow_network(t, sim, raw, current_step);

        % Instantaneous signed molar rates at each end
        inst = struct();
        for b = 1:nB
            % Concentrations at bed boundary use receiving boundary P,T.
            C0 = raw{b}.z0.P/(sim.R*raw{b}.z0.T);
            CL = raw{b}.zL.P/(sim.R*raw{b}.zL.T);

            n0 = C0 * bc0{b}.u * sim.A_bed;          % mol/s (signed, +z)
            nL = CL * bcL{b}.u * sim.A_bed;          % mol/s (signed, +z)

            s0 = bc0{b}.y(:) * n0;                   % species mol/s (signed)
            sL = bcL{b}.y(:) * nL;

            inst(b).n0 = n0; inst(b).nL = nL;
            inst(b).s0 = s0; inst(b).sL = sL;
        end

        % Trapezoidal integrate (k-1,k)
        if prev.has
            dt = Tseg(k) - Tseg(k-1);
            w  = 0.5*dt;
            for b = 1:nB
                beds_acc(b).Fz0_total   = beds_acc(b).Fz0_total   + w*(prev.inst(b).n0 + inst(b).n0);
                beds_acc(b).FzL_total   = beds_acc(b).FzL_total   + w*(prev.inst(b).nL + inst(b).nL);
                beds_acc(b).Fz0_species = beds_acc(b).Fz0_species + w*(prev.inst(b).s0 + inst(b).s0);
                beds_acc(b).FzL_species = beds_acc(b).FzL_species + w*(prev.inst(b).sL + inst(b).sL);
            end

            % Tanks (interpret endpoints in current_step)
            for b = 1:nB
                bed_name = sprintf('Bed%c','A'+b-1);

                % z0: destination tank gets bed outflow across z0 (which is negative along +z)
                if isfield(current_step.(bed_name).z0,'destination')
                    dst = current_step.(bed_name).z0.destination; if iscell(dst), dst = dst{1}; end
                    ti = tank_index_by_name(sim, dst);
                    if ~isempty(ti)
                        % instantaneous signed flow into tank at k-1 and k:
                        % bed outflow across z0 is (-n0_positive). Flow into tank is -min(n0,0)
                        fkm1 = -min(prev.inst(b).n0, 0);
                        fk   = -min(inst(b).n0, 0);
                        tanks_acc(ti).n_total   = tanks_acc(ti).n_total   + w*(fkm1 + fk);
                        skm1 = -min(prev.inst(b).s0, 0);  % elementwise
                        sk   = -min(inst(b).s0, 0);
                        tanks_acc(ti).n_species = tanks_acc(ti).n_species + w*(skm1 + sk);
                    end
                end
                % z0: source tank supplies bed inflow (positive along +z)
                if isfield(current_step.(bed_name).z0,'source')
                    src = current_step.(bed_name).z0.source; if iscell(src), src = src{1}; end
                    ti = tank_index_by_name(sim, src);
                    if ~isempty(ti)
                        fkm1 =  max(prev.inst(b).n0, 0);
                        fk   =  max(inst(b).n0, 0);
                        tanks_acc(ti).n_total   = tanks_acc(ti).n_total   - w*(fkm1 + fk);  % negative = out of tank
                        skm1 =  max(prev.inst(b).s0, 0);
                        sk   =  max(inst(b).s0, 0);
                        tanks_acc(ti).n_species = tanks_acc(ti).n_species - w*(skm1 + sk);
                    end
                end

                % zL: destination tank gets bed outflow across zL (positive along +z)
                if isfield(current_step.(bed_name).zL,'destination')
                    dst = current_step.(bed_name).zL.destination; if iscell(dst), dst = dst{1}; end
                    ti = tank_index_by_name(sim, dst);
                    if ~isempty(ti)
                        fkm1 =  max(prev.inst(b).nL, 0);
                        fk   =  max(inst(b).nL, 0);
                        tanks_acc(ti).n_total   = tanks_acc(ti).n_total   + w*(fkm1 + fk);
                        skm1 =  max(prev.inst(b).sL, 0);
                        sk   =  max(inst(b).sL, 0);
                        tanks_acc(ti).n_species = tanks_acc(ti).n_species + w*(skm1 + sk);
                    end
                end
                % zL: source tank supplies bed inflow across zL (negative along +z)
                if isfield(current_step.(bed_name).zL,'source')
                    src = current_step.(bed_name).zL.source; if iscell(src), src = src{1}; end
                    ti = tank_index_by_name(sim, src);
                    if ~isempty(ti)
                        fkm1 = -min(prev.inst(b).nL, 0);
                        fk   = -min(inst(b).nL, 0);
                        tanks_acc(ti).n_total   = tanks_acc(ti).n_total   - w*(fkm1 + fk);  % negative = out of tank
                        skm1 = -min(prev.inst(b).sL, 0);
                        sk   = -min(inst(b).sL, 0);
                        tanks_acc(ti).n_species = tanks_acc(ti).n_species - w*(skm1 + sk);
                    end
                end
            end
        end

        % store prev
        prev.inst = inst;
        prev.has  = true;
    end

    step_flow.beds  = beds_acc;
    step_flow.tanks = tanks_acc;
end

function idx = tank_index_by_name(sim, name)
    idx = [];
    if ~(ischar(name) || isstring(name)), return; end
    for k = 1:numel(sim.tanks)
        if strcmpi(string(name), sim.tanks(k).name)
            idx = k; return;
        end
    end
end

function [Fin, Fout, Sin, Sout] = end_fluxes(sim, props, bc)
% Compute instantaneous molar flows at one bed end
% Fin/Fout: scalars mol/s; Sin/Sout: n_species-by-1 mol/s components
    C = props.P / (sim.R * props.T);         % mol/m^3 at bed boundary
    n_dot = C * bc.u * sim.A_bed;            % mol/s (+ into bed)
    if n_dot >= 0
        Fin = n_dot; Fout = 0;
        Sin = bc.y(:) * n_dot;  Sout = zeros(sim.n_species,1);
    else
        Fin = 0;     Fout = -n_dot;
        Sin = zeros(sim.n_species,1);  Sout = bc.y(:) * (-n_dot);
    end
end

function x = trap_inc(value_k, Tseg, k)
% trapezoidal increment using sample k and (k-1)
    if k==1, x=0; return; end
    dt = Tseg(k) - Tseg(k-1);
    x  = 0.5 * (value_k + 0) * dt;  % value at k-1 unknown here -> recomputed upstream
    % NOTE: We call trap_inc with instantaneous value at k and k-1 implicitly
end

function vec = trap_inc_vec(value_k, Tseg, k)
    if k==1, vec=zeros(numel(value_k),1); return; end
    dt = Tseg(k) - Tseg(k-1);
    vec = 0.5 * (value_k + 0) * dt;
end

function print_step_flow_report_signed(sim, step_idx, step_flow)
    fprintf('    — Step %d signed boundary flows (mol) —\n', step_idx);

    % Beds
    fprintf('      %-6s  %-2s  %14s    %s\n', 'Bed', 'End', 'Total (mol)', 'Species (mol) [first 6]');
    for b = 1:sim.num_beds
        bedlabel = sprintf('Bed%c','A'+b-1);
        % z0 line
        sp = step_flow.beds(b).Fz0_species;
        fprintf('      %-6s  %-2s  %14.6e    ', bedlabel, 'z0', step_flow.beds(b).Fz0_total);
        print_species_vec(sp, 6);
        % zL line
        sp = step_flow.beds(b).FzL_species;
        fprintf('      %-6s  %-2s  %14.6e    ', bedlabel, 'zL', step_flow.beds(b).FzL_total);
        print_species_vec(sp, 6);
        fprintf('\n');
    end

    % Tanks
    fprintf('      \n');
    fprintf('      %-14s  %14s    %s\n', 'Tank', 'Total (mol)', 'Species (mol) [first 6]');
    for k = 1:numel(step_flow.tanks)
        fprintf('      %-14s  %14.6e    ', step_flow.tanks(k).name, step_flow.tanks(k).n_total);
        print_species_vec(step_flow.tanks(k).n_species, 6);
    end
end

function print_species_vec(v, maxshow)
    v = v(:);
    m = min(numel(v), maxshow);
    for j = 1:m
        fprintf('%10.3e', v(j));
    end
    if numel(v) > m
        fprintf('  ...');
    end
    fprintf('\n');
end

function print_cycle_totals(sim, cycle_idx, bed_cycle, tank_cycle)
    fprintf('\n=== Cycle %d totals (signed) ===\n', cycle_idx);

    % Beds
    fprintf('  %-6s  %-2s  %14s    %s\n', 'Bed', 'End', 'Total (mol)', 'Species (mol) [first 6]');
    for b = 1:sim.num_beds
        bedlabel = sprintf('Bed%c','A'+b-1);
        sp = bed_cycle(b).Fz0_species;
        fprintf('  %-6s  %-2s  %14.6e    ', bedlabel, 'z0', bed_cycle(b).Fz0_total); print_species_vec(sp, 6);
        sp = bed_cycle(b).FzL_species;
        fprintf('  %-6s  %-2s  %14.6e    ', bedlabel, 'zL', bed_cycle(b).FzL_total); print_species_vec(sp, 6);
        fprintf('\n');
    end

    % Tanks
    fprintf('  \n');
    fprintf('  %-14s  %14s    %s\n', 'Tank', 'Total (mol)', 'Species (mol) [first 6]');
    for k = 1:numel(tank_cycle)
        fprintf('  %-14s  %14.6e    ', tank_cycle(k).name, tank_cycle(k).n_total);
        print_species_vec(tank_cycle(k).n_species, 6);
    end
end
