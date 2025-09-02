function dYdt = multi_bed_ode(t, Y, sim, current_step)
% MULTI_BED_ODE Compute the time derivatives for all beds simultaneously.

num_beds = sim.num_beds;
nodes = sim.num_nodes;
states_per_bed = nodes * (1+2*sim.n_species);

dYdt = zeros(size(Y));

% Generate interface properties from state vector
raw = cell(num_beds, 1);
for b = 1:num_beds
    raw{b} = get_bed_interface_props(Y, b, sim);
end
%=== 📌 INSERT DIAGNOSTIC TRACKING HERE ===
if isfield(sim, 'diagnostics')
    for b = 1:num_beds
        % Unpack
        NS = sim.num_nodes;
        N = sim.n_species;
        offset = (b-1) * states_per_bed;

        idx_Ct = offset + (1:NS);
        idx_Ci = offset + NS + (1:NS*(N-1));
        idx_T  = offset + NS + NS*(N-1) + (1:NS);
        idx_qi = offset + NS*(N+1) + (1:NS*N);

        Ct = Y(idx_Ct);
        Ci_partial = reshape(Y(idx_Ci), NS, N-1);
        Ci = [Ci_partial, Ct - sum(Ci_partial,2)];
        T = Y(idx_T);
        qi = reshape(Y(idx_qi), NS, N);

        sim.diagnostics(b).time(end+1) = t;
        sim.diagnostics(b).Ct_1(end+1) = Ct(1);
        sim.diagnostics(b).Ct_N(end+1) = Ct(end);
        sim.diagnostics(b).y1_1(end+1) = Ci(1,1) / Ct(1);
        sim.diagnostics(b).y1_N(end+1) = Ci(end,1) / Ct(end);
        sim.diagnostics(b).P_1(end+1) = Ct(1) * sim.R * T(1);
        sim.diagnostics(b).q1_1(end+1) = qi(1,1);
        sim.diagnostics(b).q1_N(end+1) = qi(end,1);
        sim.diagnostics(b).q2_1(end+1) = qi(1,2);
        sim.diagnostics(b).q2_N(end+1) = qi(end,2);
        sim.diagnostics(b).T_1(end+1) = T(1);
        sim.diagnostics(b).T_N(end+1) = T(end);
    end
end
%=== END INSERT ===
% Compute boundary conditions at z=0 and z=L
[bc_z0, bc_zL] = flow_network(t, sim, raw, current_step);

% Loop over beds
for b = 1:num_beds
    idx_start = (b-1)*states_per_bed + 1;
    idx_end = b*states_per_bed;
    y_bed = Y(idx_start:idx_end);

    dy = bed_ode(t, y_bed, sim, bc_z0{b}, bc_zL{b}, b);
    dYdt(idx_start:idx_end) = dy;
end

end
