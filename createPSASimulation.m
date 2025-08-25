function [sim, bed_states, tank_states] = createPSASimulation()
% CREATEPSASIMULATION Initializes the PSA simulation configuration and states.
% Includes optional restart from previous simulation data.

%% === General Simulation Settings ===
sim.n_cycles = 1;
sim.cycle_time = 200;
sim.num_beds = 2;

sim.species_names = {'O2', 'CO2', 'N2'};
sim.n_species = numel(sim.species_names);
sim.species_index = containers.Map(sim.species_names, 1:sim.n_species);

%% === Load Gas Properties ===
load('gas_database.mat', 'gas_database');
for i = 1:sim.n_species
    match = find(strcmp({gas_database.name}, sim.species_names{i}), 1);
    if isempty(match)
        error('Species %s not found in gas_database.mat', sim.species_names{i});
    else
        sim.gas(i) = gas_database(match);
    end
end

%% === Bed & Layer Configuration ===
sim.bed_diameter = 0.041;
sim.n_layers = 1;
sim.layers(1).length = 0.9;
sim.layers(1).num_nodes = 20;
sim.layers(1).adsorbent_name = 'ActivatedCarbon_1';
sim.bed_length = sum([sim.layers.length]);
sim.num_nodes = sum([sim.layers.num_nodes]);

% Load adsorbents
load('adsorbent_database_new.mat', 'adsorbent_database');
if iscell(adsorbent_database), adsorbent_database = [adsorbent_database{:}]; end

for i = 1:sim.n_layers
    name = sim.layers(i).adsorbent_name;
    match = find(strcmp({adsorbent_database.name}, name), 1);
    if isempty(match)
        error('Adsorbent %s not found.', name);
    else
        sim.layers(i).properties = adsorbent_database(match);
    end
end

% Assign node indices for each layer
node_idx = 1;
for i = 1:sim.n_layers
    sim.layers(i).node_start = node_idx;
    sim.layers(i).node_end = node_idx + sim.layers(i).num_nodes - 1;
    node_idx = sim.layers(i).node_end + 1;
end

%% === Tank Definitions ===
sim.tanks(1).name = 'Feed_Tank';
sim.tanks(1).type = 'infinite';
sim.tanks(1).P = 5e5;
sim.tanks(1).T = 298.15;
sim.tanks(1).y = [0.21, 0.0004, 0.7896];

sim.tanks(2).name = 'Product_Tank';
sim.tanks(2).type = 'finite';
sim.tanks(2).volume = 10.0;
sim.tanks(2).P = 1.5e5;
sim.tanks(2).T = 298.15;
sim.tanks(2).y = [1.0, 0.0, 0.0];

sim.tanks(3).name = 'Vent_Tank';
sim.tanks(3).type = 'infinite';
sim.tanks(3).P = 1e5;
sim.tanks(3).T = 298.15;
sim.tanks(3).y = [0.0, 1.0, 0.0];

%% === Restart Option ===
sim.restart_from_file = false;
sim.restart_file = 'cycle_sim.mat';

if sim.restart_from_file
    fprintf('Loading initial state from %s...\n', sim.restart_file);
    data = load(sim.restart_file);

    if isfield(data, 'bed_states')
        bed_states = data.bed_states;
    else
        error('Restart file missing bed_states.');
    end

    if isfield(data, 'tank_states')
        tank_states = data.tank_states;
    else
        warning('Restart file missing tank_states. Reinitializing tanks...');
        tank_states = initializeTankStates(sim);
    end

    fprintf('Restart loaded.\n');
else
    % Assign default initial conditions
    for i = 1:sim.num_beds
        sim.init_conditions{i}.P0 = 1.01e5 * ones(sim.num_nodes, 1);
        sim.init_conditions{i}.T0 = 298 * ones(sim.num_nodes, 1);
        sim.init_conditions{i}.y0 = repmat([0.2, 0.00, 0.80], sim.num_nodes, 1);
    end

    [sim, bed_states, tank_states] = initializeFromFresh(sim);
end
end
