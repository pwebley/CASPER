function casper_check_RHS_at_t0(sim, Y0, step_config)
% Verifies multi_bed_ode at t=0 is finite, and pinpoints NaNs/Infs.

t0 = 0.0;

% 1) Build BCs exactly like the solver does
raw = cell(sim.num_beds,1);
for b = 1:sim.num_beds
    raw{b} = get_bed_interface_props(Y0, b, sim);
end
[bc_z0, bc_zL] = flow_network(t0, sim, raw, step_config);

% 2) Evaluate each bed’s RHS and check finiteness
NS = sim.num_nodes; N = sim.n_species;
varlen = NS + NS*(N-1) + NS + NS*N;

any_bad = false;
for b = 1:sim.num_beds
    s = (b-1)*varlen+1; e = b*varlen;
    yb = Y0(s:e);
    dy = bed_ode(t0, yb, sim, bc_z0{b}, bc_zL{b}, b);

    if any(~isfinite(yb)),  fprintf('❌ Non-finite initial state in Bed %c\n','A'+b-1); any_bad=true; end
    if any(~isfinite(dy)),  fprintf('❌ Non-finite dy at t=0 in Bed %c\n','A'+b-1);    any_bad=true; end
end

if ~any_bad
    fprintf('✅ RHS and BCs finite at t=0.\n');
end
end
