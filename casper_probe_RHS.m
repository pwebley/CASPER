function casper_probe_RHS(sim, Y0, steps, t_list)
% Integrates from 0 to max(t_list) with dynamic step configs,
% forces outputs at t_list, then checks RHS/BC at each probe time.

t_list = unique(t_list(:).');        % row vector, sorted
t_end  = t_list(end);

% Dynamic RHS that selects the proper step config for each t
f = @(t,y) local_rhs(t,y,sim,steps);

% Solver options (mild caps to help through kinks)
opts = odeset('RelTol',1e-6,'AbsTol',1e-8, ...
              'InitialStep',1e-4,'MaxStep',0.1, ...
              'OutputFcn',[],'Stats','on');

% IMPORTANT: pass TSPAN as the vector of probe times
[t_out, Y_out] = ode23t(f, t_list, Y0, opts);  %#ok<*ASGLU>

% For each requested time, check RHS/BC using the *returned* state
for i = 1:numel(t_out)
    ti = t_out(i);
    Yi = Y_out(i,:).';
    [k, cfg] = casper_step_for_time(sim, steps, ti);
    fprintf('\n-- Probe @ t = %.6g s (step %d) --\n', ti, k);
    casper_check_RHS_at_t_state(sim, Yi, cfg, ti);
end

    function dYdt = local_rhs(t, Y, sim_local, steps_local)
        % Pick correct step for current t
        [~, cfg_now] = casper_step_for_time(sim_local, steps_local, t);

        % Standard multi-bed assembly
        NS = sim_local.num_nodes; N = sim_local.n_species;
        varlen = NS + NS*(N-1) + NS + NS*N;
        dYdt = zeros(size(Y));
        % Build BCs once per call
        raw = cell(sim_local.num_beds,1);
        for b = 1:sim_local.num_beds
            raw{b} = get_bed_interface_props(Y, b, sim_local);
        end
        [bc_z0, bc_zL] = flow_network(t, sim_local, raw, cfg_now);

        for b = 1:sim_local.num_beds
            s = (b-1)*varlen+1; e = b*varlen;
            dYdt(s:e) = bed_ode(t, Y(s:e), sim_local, bc_z0{b}, bc_zL{b}, b);
        end
    end
end
