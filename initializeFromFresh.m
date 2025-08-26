function [sim, bed_states, tank_states] = initializeFromFresh(sim)
%INITIALIZEFROMFRESH Initializes grid, bed, and tank states from scratch.



%% Bed State Initialization
fprintf('Initializing bed states...\n');
R = 8.314462618; % J/mol·K
bed_states = cell(sim.num_beds, 1);

for bed_id = 1:sim.num_beds
    P = sim.init_conditions{bed_id}.P0;
    T = sim.init_conditions{bed_id}.T0;
    y = sim.init_conditions{bed_id}.y0;

    bed.P = P;
    bed.T = T;
    bed.y = y;

    C_total = P ./ (R .* T);
    C_species = y .* C_total;
    bed.C = C_total;
    bed.C_species = C_species(:, 1:end-1);

    bed.q = zeros(sim.num_nodes, sim.n_species);
    for node = 1:sim.num_nodes
        layer_id = sim.layer_id(node);
        ads = sim.layers(layer_id).properties;
        bed.q(node, :) = get_equilibrium_loading(sim, ads, P(node) * y(node, :), T(node));
    end

    bed.u_z0 = 0;
    bed.u_zL = 0;
    bed.current_step_config = [];

    bed_states{bed_id} = bed;
end

%% Tank Initialization
fprintf('Initializing tank states...\n');
tank_states = initializeTankStates(sim);
sim.num_tanks = numel(sim.tanks);
fprintf('Initialization complete: %d beds, %d tanks, %d nodes.\n', ...
        sim.num_beds, sim.num_tanks, sim.num_nodes);
end
