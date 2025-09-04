
function run_casper_simulation_diag
% RUN_CASPER_SIMULATION_DIAG
% Drop-in driver that adds robust per-step diagnostics and printing.
%
% It uses your existing createPSASimulation, flow_network, get_bed_interface_props,
% bed_ode, and unpack_bed_state_vector functions.
%
% Usage: run_casper_simulation_diag

    %% Build sim and initial states
    [sim, bed_states, ~] = createPSASimulation();
    Y0 = pack_initial_state(sim, bed_states);
    t_now = 0.0;
         
    
    %% Step schedule (support sim.step / sim.step_times OR sim.steps)
    steps = [];
    if isfield(sim, 'steps') && ~isempty(sim.steps)
        steps = sim.steps;
    elseif isfield(sim, 'step') && ~isempty(sim.step)
        % Build a uniform steps struct: name, duration, config
        if isfield(sim,'step_times') && numel(sim.step_times) >= numel(sim.step)+1
            t_edges = sim.step_times(:);
            durs = diff(t_edges);
        else
            % If no step_times provided, use equal durations over sim.cycle_time
            durs = repmat(sim.cycle_time/numel(sim.step), numel(sim.step), 1);
        end
        steps = repmat(struct('name',"", 'duration',0, 'config',[]), numel(sim.step), 1);
        for i = 1:numel(sim.step)
            steps(i).name = sprintf("Step%d", i);
            steps(i).duration = durs(i);
            steps(i).config = sim.step(i); % pass through the per-bed config
        end
    else
        % Fallback to a single step using sim.cycle_time if nothing defined
        steps = struct('name',"Step1",'duration',sim.cycle_time,'config',[]);
    end
    nSteps = numel(steps);
    nCycles = sim.n_cycles;
% --- DEBUG PROBE (optional): check RHS/BC at specific times before full run
t_list = [0, 0.1:0.1:2.0];        % probe 0→2 s densely
casper_probe_RHS(sim, Y0, steps, t_list);


    %% ODE options with combined OutputFcn (progress + live diagnostics)
    diag = initDiagnostics(sim);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-9, ...
                  'OutputFcn', @(t,y,flag) combinedOutputFcn(t,y,flag, sim));
    opts = odeset(opts,'OutputFcn',[]);

    fprintf('\n=== CASPER: %d bed(s), %d species, %d node(s) ===\n', sim.num_beds, sim.n_species, sim.num_nodes);
    for cyc = 1:nCycles
        fprintf('\n=== Beginning cycle %d of %d ===\n', cyc, nCycles);

        for s = 1:nSteps
            step = steps(s);
            tspan = [t_now, t_now + step.duration];
            fprintf('  Cycle %d — %s (step %d/%d): t = %.2f → %.2f s\n', ...
                    cyc, step.name, s, nSteps, tspan(1), tspan(2));

            % Prepare step configuration for the network
            if isstruct(step.config) && ~isempty(fieldnames(step.config))
                current_step = step.config;
            else
                % Minimal empty config if not provided
                current_step = struct();
            end

            % Solve this step
            f = @(t,y) multi_bed_ode(t, y, sim, current_step);
            [t_sol, Y_sol] = ode23t(f, tspan, Y0, opts); %#ok<ASGLU>
            Y_end = Y_sol(end, :)';

            % Compute BCs at end of step to print meaningful summary
            raw = cell(sim.num_beds, 1);
            for b = 1:sim.num_beds
                raw{b} = get_bed_interface_props(Y_end, b, sim);
            end
            try
                [bc_z0, bc_zL] = flow_network(tspan(2), sim, raw, current_step);
            catch
                % If flow_network depends on live state only, still continue printing
                bc_z0 = cell(sim.num_beds, 1);
                bc_zL = cell(sim.num_beds, 1);
                for b = 1:sim.num_beds, bc_z0{b} = struct(); bc_zL{b} = struct(); end
            end

            % Per-step printout
            print_step_summary(cyc, s, step.name, sim, Y_end, tspan, bc_z0, bc_zL);

            % Carry to next step
            Y0 = Y_end;
            t_now = tspan(2);
        end
    end

    fprintf('\n✅ Simulation complete.\n');

    % Optional: quick plots
    try
        plot_diagnostics(sim);
    catch
    end

    %% ===== Nested helpers =====
    function status = combinedOutputFcn(t, y, flag, sim_local)
        % Text progress
        status1 = odeProgress(t, y, flag);
        % Append basic live diagnostics every accepted step
        if isempty(flag)
            % Nothing heavy here to keep the solver snappy
            % (Detailed logging is done from the per-step summary above)
        end
        status = status1 || false;
    end

end

function Y = pack_initial_state(sim, bed_states)
    % PACK_INITIAL_STATE — Packs {Ct, Ci(1..N-1), T, q(1..N)} for each bed
    NS = sim.num_nodes; N = sim.n_species;
    varlen = NS + NS*(N-1) + NS + NS*N;
    Y = zeros(sim.num_beds * varlen, 1);
    for b = 1:sim.num_beds
        bed = bed_states{b};
        Ct = bed.C(:);
        y = bed.y;  % NS x N
        Ci = y(:,1:N-1) .* Ct; Ci = reshape(Ci, [], 1);
        T  = bed.T(:);
        qi = bed.q; qi = reshape(qi, [], 1);
        idx0 = (b-1)*varlen;
        Y(idx0+1 : idx0+NS) = Ct;
        Y(idx0+NS+1 : idx0+NS+NS*(N-1)) = Ci;
        Y(idx0+NS+NS*(N-1)+1 : idx0+NS+NS*(N-1)+NS) = T;
        Y(idx0+NS+NS*(N-1)+NS+1 : idx0+varlen) = qi;
    end
end

function diag = initDiagnostics(sim)
    % INITDIAGNOSTICS — keeps space for future expansion
    for b = 1:sim.num_beds
        diag(b).time = [];
        diag(b).Pin = []; diag(b).Pout = [];
        diag(b).y1_in = []; diag(b).y1_out = [];
        diag(b).u_in = []; diag(b).u_out = [];
    end
end
