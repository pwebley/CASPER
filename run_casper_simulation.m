function run_casper_simulation(run_name)
% RUN_CASPER_SIMULATION  PSA simulation with per-step reporting + logging.
%
% Usage:
%   run_casper_simulation                % auto run name (timestamp)
%   run_casper_simulation('baseline_A')  % custom run name
%
% Produces a self-contained results file:
%   casper_run_<run_name>__YYYYMMDD_HHMMSS.mat
%
% Requires on path:
%   createPSASimulation, multi_bed_ode, unpack_bed_state_vector,
%   flow_network, get_bed_interface_props, integrate_step_flows_signed,
%   (and your other CASPER functions)

if nargin < 1 || isempty(run_name)
    run_name = strtrim(input('Enter run name (e.g. Baseline_O2N2): ', 's'));
    if isempty(run_name), run_name = 'run'; end
end
% sanitize (spaces -> underscores; keep alnum, dash, underscore)
run_name = regexprep(run_name, '\s+', '_');
run_name = regexprep(run_name, '[^A-Za-z0-9_\-]', '_');

% ensure runs/ exists
runs_dir = 'runs';
if ~exist(runs_dir, 'dir'), mkdir(runs_dir); end

tic;
[sim, bed_states, tank_states] = createPSASimulation();

% light diagnostics (unchanged)
for b = 1:sim.num_beds
    sim.diagnostics(b) = struct('time',[],'Pin',[],'Pout',[], ...
        'y1_in',[],'y1_out',[],'u_in',[],'u_out',[]);
end

num_cycles = sim.n_cycles;

% Build initial state Y0
Y0 = [];
for b = 1:sim.num_beds
    bed = bed_states{b};
    Y0 = [Y0; bed.C; reshape(bed.C_species, [], 1); bed.T; reshape(bed.q, [], 1)];
end

% ODE options
opts = odeset('RelTol', sim.RelTol, 'AbsTol', sim.AbsTol, 'MaxStep', sim.dt_max);

% --- NEW: logger ---
logger = casper_logger_start(run_name, sim);
sim.log_boundary_history = true;

% --- storage for (optional) raw time/state ---
T_all = []; Y_all = [];

for cycle = 1:num_cycles
    fprintf('\n=== Beginning cycle %d of %d ===\n', cycle, num_cycles);

    % init per-cycle accumulators
    bed_cycle  = init_bed_accum(sim);
    tank_cycle = init_tank_accum(sim);

    % mark the start of the cycle in logger
    logger = casper_logger_mark_cycle_start(logger, cycle, Y0);

    for idx = 1:length(sim.step_times)-1
        t0 = sim.step_times(idx);
        tf = sim.step_times(idx+1);
        current_step = sim.step(idx);

        fprintf('  Cycle %d — Step %d: t = %.2f to %.2f s\n', cycle, idx, t0, tf);

        % Integrate
        solfcn = @(t,Y) multi_bed_ode(t, Y, sim, current_step);
        [Tseg, Yseg] = feval(sim.ode_solver, solfcn, [t0, tf], Y0, opts);

        % append raw (optional)
        T_all = [T_all; Tseg];
        Y_all = [Y_all; Yseg];
        % Keep per-step trajectories for postprocessing
        step_struct = struct();
        step_struct.idx  = idx;             % step number within cycle
        step_struct.step = current_step;    % the config used
        step_struct.Tseg = Tseg(:);
        step_struct.Yseg = Yseg;

        if ~exist('last_cycle','var') || ~isstruct(last_cycle)
            last_cycle = struct('steps', step_struct, ...
                'step_times', sim.step_times(:)-sim.step_times(1), ...
                'hist', struct());
        else
            last_cycle.steps(end+1) = step_struct;
        end

        % update state for next step
        Y0 = Yseg(end, :)';
        if sim.log_boundary_history
            % Store for last_cycle.hist
            step_struct.step = current_step;
            step_struct.Tseg = Tseg;      % keep trajectories
            step_struct.Yseg = Yseg;

            if ~exist('last_cycle','var')
                last_cycle = struct('steps', step_struct, 'hist', struct(), 'step_times', sim.step_times(:)-sim.step_times(1));
            else
                last_cycle.steps(end+1) = step_struct;
            end
        end

        % light boundary logging
        beds = unpack_bed_state_vector(Y0, sim);
        for b = 1:sim.num_beds
            P = beds{b}.C .* sim.R .* beds{b}.T;
            sim.diagnostics(b).time(end+1)  = tf;
            sim.diagnostics(b).Pin(end+1)   = P(1);
            sim.diagnostics(b).Pout(end+1)  = P(end);
            sim.diagnostics(b).y1_in(end+1) = beds{b}.y(1,1);
            sim.diagnostics(b).y1_out(end+1)= beds{b}.y(end,1);
        end

        % ---- integrate signed boundary flows & print
        step_flow = integrate_step_flows_signed(sim, current_step, Tseg, Yseg);
        print_step_flow_report_signed(sim, idx, step_flow);

        % accumulate to cycle totals
        bed_cycle  = accum_bed_signed(bed_cycle, step_flow.beds);
        tank_cycle = accum_tank(tank_cycle, step_flow.tanks);

        % ---- logger: store step summaries (time histories + Y_end)
        logger = casper_logger_step(logger, sim, idx, Tseg, Yseg, current_step, step_flow);

        % plots (your lightweight function)
        try
            % plot_diagnostics(sim, Y0);
            drawnow limitrate
        catch
        end
    end

    % end-of-cycle reporting
    print_cycle_totals(sim, cycle, bed_cycle, tank_cycle);

    % logger: add cycle summary (turn on store_profiles if you want evolution)
    logger = casper_logger_cycle_end(logger, cycle, bed_cycle, tank_cycle, Y0, ...
                                     false, []);  % set true,[steps] to track evolution
end
toc;

% Optional: keep originals, plus last state
bed_states_all{cycle}  = unpack_bed_state_vector(Y0, sim); %#ok<NASGU>
tank_states_all{cycle} = tank_states; %#ok<NASGU>
Y_end = Y0; %#ok<NASGU>
save(run_name, 'sim', 'bed_states_all', 'tank_states_all', 'last_cycle', '-v7.3');

% legacy save (optional)
save('cycle_sim.mat','sim','bed_states_all','tank_states_all','Y_end');

% Timestamped name (from your logger)
casper_logger_save(logger);  % keeps casper_run_<run_name>__YYYYMMDD_HHMMSS.mat
% casper_logger_save already wrote: logger.run_info.saved_as (timestamped)
ts_file = logger.run_info.saved_as;

% Also keep convenient copies:
non_ts_file = fullfile(runs_dir, sprintf('casper_run_%s_latest.mat', run_name));
named_file  = fullfile(runs_dir, sprintf('casper_run_%s.mat',        run_name));

% Auto-write rich report (next to the timestamped run file, in runs/)
try
    report_dir = 'runs';
    if ~exist(report_dir,'dir'), mkdir(report_dir); end
    report_file = fullfile(report_dir, sprintf('casper_report_%s__%s.txt', logger.run_info.run_name, logger.run_info.timestamp));
    casper_write_report(sim, logger, report_file);
catch ME
    warning('casper:reportFail', 'Report generation failed: %s', ME.message);
end

% try/catch to avoid halting if copy fails
try
    copyfile(ts_file, non_ts_file);
    copyfile(ts_file, named_file);
    fprintf('📦 Saved copies:\n  %s\n  %s\n', non_ts_file, named_file);
catch ME
    warning('casper:copyFail', 'Could not copy run file to runs/: %s', ME.message);
end

% Also save convenient names:

copyfile(logger.run_info.saved_as, non_ts_file);

copyfile(logger.run_info.saved_as, named_file);

fprintf('📦 Saved copies:\n  %s\n  %s\n', non_ts_file, named_file);

fprintf('\n✅ Simulation complete.\n');
end


function status = combinedOutputFcn(t, y, flag)
    status1 = odeProgress(t, y, flag);
    status2 = odeplot(t, y, flag);
    status = status1;
    status = status1 || status2;
end



function S = init_bed_accum(sim)
S = struct();
for b = 1:sim.num_beds
    S(b).Fz0_total = 0; S(b).FzL_total = 0;
    S(b).Fz0_species = zeros(sim.n_species,1);
    S(b).FzL_species = zeros(sim.n_species,1);
end
end

function S = init_tank_accum(sim)
S = struct();
for k = 1:numel(sim.tanks)
    S(k).name      = sim.tanks(k).name;
    S(k).n_total   = 0;
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
% ========================= LOGGER (new) ========================= %

function logger = casper_logger_start(run_name, sim)
ts = datestr(now, 'yyyymmdd_HHMMSS');
logger.run_info.run_name  = char(run_name);
logger.run_info.timestamp = ts;
logger.run_info.saved_as  = sprintf('casper_run_%s__%s.mat', run_name, ts);
logger.sim = sim;

logger.last_cycle.steps = struct([]);
logger.last_cycle.Y_start_cycle = [];
logger.last_cycle.cycle_index = NaN;

logger.cycle_summaries = struct('cycle_index',{},'bed_cycle',{},'tank_cycle',{});

logger.profile_snapshots.enabled = false;
logger.profile_snapshots.step_indices = [];
logger.profiles_over_cycles = struct('cycle_index',{},'step_index',{},'Y_end',{});
end

function logger = casper_logger_mark_cycle_start(logger, cycle_idx, Y0)
logger.last_cycle.cycle_index = cycle_idx;
logger.last_cycle.steps = struct([]);
logger.last_cycle.Y_start_cycle = Y0(:);
end

function logger = casper_logger_step(logger, sim, step_idx, Tseg, Yseg, current_step, step_flow)
% Store decimated boundary histories and end-of-step state for one step.

    % --- Decimate to ~400 samples to keep files light ---
    keep = max(1, floor(numel(Tseg)/400));
    Ti = Tseg(1:keep:end);
    Yi = Yseg(1:keep:end, :).';   % [ny x nt]

    nB  = sim.num_beds;
    nS  = sim.n_species;
    Nt  = numel(Ti);

    hist = struct();
    hist.t = Ti(:);

    % Preallocate per-bed, per-end containers
    for b = 1:nB
        bedlabel = sprintf('Bed%c', 'A'+b-1);

        % z0
        hist.(bedlabel).z0.P = zeros(Nt,1);
        hist.(bedlabel).z0.T = zeros(Nt,1);
        hist.(bedlabel).z0.u = zeros(Nt,1);
        hist.(bedlabel).z0.n = zeros(Nt,1);
        hist.(bedlabel).z0.y = zeros(Nt,nS);

        % zL
        hist.(bedlabel).zL.P = zeros(Nt,1);
        hist.(bedlabel).zL.T = zeros(Nt,1);
        hist.(bedlabel).zL.u = zeros(Nt,1);
        hist.(bedlabel).zL.n = zeros(Nt,1);
        hist.(bedlabel).zL.y = zeros(Nt,nS);
    end

    % Time samples
    for k = 1:Nt
        % Build interface states for ALL beds at this time sample
        raw_all = cell(nB,1);
        for bb = 1:nB
            raw_all{bb} = get_bed_interface_props(Yi(:,k), bb, sim);
        end

        % Compute BCs for all beds at this time
        [bc0_all, bcL_all] = flow_network(Ti(k), sim, raw_all, current_step);

        % Fill per-bed histories
        for b = 1:nB
            bedlabel = sprintf('Bed%c', 'A'+b-1);
            raw_b = raw_all{b};
            bc0b  = bc0_all{b};
            bcLb  = bcL_all{b};

            % Concentration at boundary uses the RECEIVING boundary T
            % (consistent with signed flow definitions)
            C0 = raw_b.z0.P/(sim.R*bc0b.T);
            CL = raw_b.zL.P/(sim.R*bcLb.T);

            n0 = C0 * bc0b.u * sim.A_bed;   % mol/s (+ along +z)
            nL = CL * bcLb.u * sim.A_bed;

            % z0
            hist.(bedlabel).z0.P(k)   = raw_b.z0.P;
            hist.(bedlabel).z0.T(k)   = bc0b.T;
            hist.(bedlabel).z0.u(k)   = bc0b.u;
            hist.(bedlabel).z0.n(k)   = n0;
            hist.(bedlabel).z0.y(k,:) = bc0b.y(:).';

            % zL
            hist.(bedlabel).zL.P(k)   = raw_b.zL.P;
            hist.(bedlabel).zL.T(k)   = bcLb.T;
            hist.(bedlabel).zL.u(k)   = bcLb.u;
            hist.(bedlabel).zL.n(k)   = nL;
            hist.(bedlabel).zL.y(k,:) = bcLb.y(:).';
        end
    end

    % End-of-step packed state
    Y_end = Yseg(end,:).';

    % Stash into logger
    logger.last_cycle.steps(step_idx).Tseg      = Ti;
    logger.last_cycle.steps(step_idx).hist      = hist;
    logger.last_cycle.steps(step_idx).Y_end     = Y_end;
    logger.last_cycle.steps(step_idx).step_flow = step_flow;
end

function logger = casper_logger_cycle_end(logger, cycle_idx, bed_cycle, tank_cycle, Y_end, store_profiles, which_steps)
i = numel(logger.cycle_summaries) + 1;
logger.cycle_summaries(i).cycle_index = cycle_idx;
logger.cycle_summaries(i).bed_cycle = bed_cycle;
logger.cycle_summaries(i).tank_cycle = tank_cycle;

if nargin>=6 && store_profiles
    logger.profile_snapshots.enabled = true;
    logger.profile_snapshots.step_indices = which_steps(:).';
    for s = which_steps(:).'
        j = numel(logger.profiles_over_cycles) + 1;
        logger.profiles_over_cycles(j).cycle_index = cycle_idx;
        logger.profiles_over_cycles(j).step_index  = s;
        logger.profiles_over_cycles(j).Y_end       = Y_end(:);
    end
end
end

function casper_logger_save(logger)
save(logger.run_info.saved_as, '-struct', 'logger');
fprintf('📦 Saved run to %s\n', logger.run_info.saved_as);
end