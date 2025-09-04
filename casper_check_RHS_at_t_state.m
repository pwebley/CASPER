function casper_check_RHS_at_t_state(sim, Y_state, step_config, t_now)
% Rebuild BCs at t_now, then evaluate dy for each bed and check finiteness.

% Build BCs exactly like solver
raw = cell(sim.num_beds,1);
for b = 1:sim.num_beds
    raw{b} = get_bed_interface_props(Y_state, b, sim);
end
[bc_z0, bc_zL] = flow_network(t_now, sim, raw, step_config);

% Unpack per-bed slices and test dy
NS = sim.num_nodes; N = sim.n_species;
varlen = NS + NS*(N-1) + NS + NS*N;
any_bad = false;

for b = 1:sim.num_beds
    s = (b-1)*varlen+1; e = b*varlen;
    yb = Y_state(s:e);

    if any(~isfinite(yb))
        fprintf('❌ Non-finite state at t=%.6g in Bed %c\n', t_now, 'A'+b-1);
        any_bad = true;
    end

    dy = bed_ode(t_now, yb, sim, bc_z0{b}, bc_zL{b}, b);

    if any(~isfinite(dy))
        fprintf('❌ Non-finite dy at t=%.6g in Bed %c\n', t_now, 'A'+b-1);
        any_bad = true;
    end
end

if ~any_bad
    fprintf('✅ RHS and BCs finite at t=%.6g.\n', t_now);
end
end
