function dYdt = multi_bed_ode(t, Y, sim, current_step)
% MULTI_BED_ODE Compute the time derivatives for all beds simultaneously.

num_beds = sim.num_beds;
nodes = sim.num_nodes;
states_per_bed = nodes * (1+2*n_species);

dYdt = zeros(size(Y));

% Generate interface properties from state vector
raw = cell(num_beds, 1);
for b = 1:num_beds
    raw{b} = get_bed_interface_props(Y, b, sim);
end

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
