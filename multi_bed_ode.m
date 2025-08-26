function dYdt = multi_bed_ode(t, Y, sim, current_step)
% MULTI_BED_ODE system of ODEs for all beds in the PSA network.

num_beds = sim.num_beds;
nodes = sim.num_nodes;
n_species = sim.n_species;

% States: Ct (nodes), Ci (nodes*(n_species-1)), T (nodes), qi (nodes*n_species)
n_per_bed = nodes * (2 + n_species);

dYdt = zeros(size(Y));

% Compute boundary conditions for all beds
flow_bc = flow_network(sim, Y, current_step);

for bed_idx = 1:num_beds
    start_idx = (bed_idx - 1)*n_per_bed + 1;
    end_idx   = bed_idx * n_per_bed;
    y_bed = Y(start_idx:end_idx);

    inlet_bc = flow_bc(bed_idx).inlet;
    outlet_bc = flow_bc(bed_idx).outlet;

    dy_bed = bed_ode(t, y_bed, sim, inlet_bc, outlet_bc, bed_idx);

    dYdt(start_idx:end_idx) = dy_bed;
end
end
