function [sim, bed_states, tank_states] = createPSASimulation()
% CREATEPSASIMULATION Initializes the PSA simulation configuration and states.
% Includes optional restart from previous simulation data.

%% === General Simulation Settings ===
sim.n_cycles = 2;
sim.cycle_time = 200;
sim.num_beds = 2;

sim.species_names = {'O2', 'N2'};
sim.n_species = numel(sim.species_names);
sim.species_index = containers.Map(sim.species_names, 1:sim.n_species);
sim.R = 8.314462618;
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
sim.visc_func = @(T) (1.8e-5) * (T / 300).^0.7;  % [Pa·s], Sutherland-like law
% Extract molecular weights from gas properties
sim.MW = zeros(sim.n_species, 1);
sim.cp_g = zeros(sim.n_species, 1);

for i = 1:sim.n_species
    sim.MW(i) = sim.gas(i).MW;  % assumed to be in kg/mol
    sim.cp_g(i)=sim.gas(i).Cp; 
end

%% === Bed & Layer Configuration ===
sim.bed_diameter = 0.05;
sim.n_layers = 1;
sim.layers(1).length = 1.0;
sim.layers(1).num_nodes = 20;
sim.Ergun_A_face = ones(sim.layers(1).num_nodes+1,1);
sim.Ergun_B_face = ones(sim.layers(1).num_nodes+1,1);

sim.layers(1).adsorbent_name = 'Zeolite5A_1';
sim.bed_length = sum([sim.layers.length]);
sim.num_nodes = sum([sim.layers.num_nodes]);

% Load adsorbents
load('adsorbent_database.mat', 'adsorbent_database');
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
% Initialize sim.k_LDF as [num_nodes x n_species] matrix
sim.k_LDF = zeros(sim.num_nodes, sim.n_species);

% Assign kinetics per node
for i = 1:sim.n_layers
    nodes = sim.layers(i).node_start : sim.layers(i).node_end;
    kinetics = sim.layers(i).properties.kinetics;
    
    for j = 1:sim.n_species
        gas_name = sim.species_names{j};
        if isfield(kinetics, gas_name) && isfield(kinetics.(gas_name), 'k_LDF')
            sim.k_LDF(nodes, j) = kinetics.(gas_name).k_LDF;
        else
            error('Missing k_LDF for species %s in layer %d', gas_name, i);
        end
    end
end
%% Thermal Model
sim.heat_source=0.0;
sim.heat_sink=0.0;
%% Initialize grid
fprintf('Initializing grid...\n');

sim.z_nodes = zeros(sim.num_nodes, 1);
sim.dz = zeros(sim.num_nodes, 1);
sim.layer_id = zeros(sim.num_nodes, 1);

% Construct the grid
node_index = 1; current_z = 0;
for layer_id = 1:sim.n_layers
    L = sim.layers(layer_id).length;
    N = sim.layers(layer_id).num_nodes;
    dz = L / N;

    for n = 1:N
        sim.z_nodes(node_index) = current_z + (n - 0.5) * dz;
        sim.dz(node_index) = dz;
        sim.layer_id(node_index) = layer_id;
        node_index = node_index + 1;
    end

    current_z = current_z + L;
end

sim.z_L = current_z;
sim.A_bed = pi * (sim.bed_diameter / 2)^2;
sim.node_volume = sim.A_bed * sim.dz;

% Compute adsorbent mass
sim.adsorbent_mass = zeros(sim.num_nodes, 1);
for i = 1:sim.n_layers
    nodes = sim.layers(i).node_start : sim.layers(i).node_end;
    rho = sim.layers(i).properties.bed_density;

    sim.adsorbent_mass(nodes) = rho .* sim.node_volume(nodes);
end
sim.particle_diameter = zeros(sim.num_nodes, 1);
sim.epsilon = zeros(sim.num_nodes, 1);
sim.rho_bed = zeros(sim.num_nodes, 1);
sim.Cp_solid = zeros(sim.num_nodes, 1);  
sim.dH = zeros(sim.num_nodes,sim.n_species);

for i = 1:sim.n_layers
    nodes = sim.layers(i).node_start : sim.layers(i).node_end;
    
    % Extract the adsorbent properties for this layer
    props = sim.layers(i).properties;
    sim.epsilon(nodes) = props.void_fraction;
    sim.rho_bed(nodes) = props.bed_density;
    sim.Cp_solid(nodes) = props.heat_capacity;
    sim.particle_diameter(nodes) = props.particle_diameter;
    for s = 1:sim.n_species
        species = sim.species_names{s};
        sim.dH(nodes, s) = props.isotherms.(species).H_ads;  % [J/mol]
    end
end
 eps=sim.epsilon;
 dp=sim.particle_diameter;
 for f = 2:sim.num_nodes
    eps_face = 0.5 * (eps(f-1) + eps(f));
    dp_face  = 0.5 * (dp(f-1) + dp(f));
    A_face(f) = (150 * (1 - eps_face)^2) / (eps_face^3 * dp_face^2);
    B_face(f) = (1.75 * (1 - eps_face))   / (eps_face^3 * dp_face);
 end
 sim.Ergun_A_face = A_face;
 sim.Ergun_B_face = B_face;



%% === Tank Definitions ===
sim.tanks(1).name = 'Feed_Tank';
sim.tanks(1).type = 'infinite';
sim.tanks(1).P = 5e5;
sim.tanks(1).T = 298.15;
sim.tanks(1).y = [0.21, 0.79];

sim.tanks(2).name = 'Product_Tank';
sim.tanks(2).type = 'finite';
sim.tanks(2).volume = 1.0;
sim.tanks(2).P = 4.5e5;
sim.tanks(2).T = 298.15;
sim.tanks(2).y = [1.0, 0.0];

sim.tanks(3).name = 'Vent_Tank';
sim.tanks(3).type = 'infinite';
sim.tanks(3).P = 1e5;
sim.tanks(3).T = 298.15;
sim.tanks(3).y = [0.0,  1.0];
%% === CYCLE DEFINITION ===
% Define the duration of each step in the cycle
sim.step_times = [0, 40, 55, 60, 75, 115, 130, 135, 150]; % 8 steps

% Pre-allocate a structure array for steps
sim.step(1:8) = struct(); % For an 8-step cycle

%----------------------------------------------------------------------
% Step 1: Bed A in State 4 (Adsorption), Bed B in State 3 (Blowdown)
%----------------------------------------------------------------------
i = 1;
% --- Bed A Configuration (State 4: u_z0 > 0, u_zL > 0) ---
sim.step(i).BedA.state = 4;
% z=0: Flow IN. Source is Feed_Tank.
sim.step(i).BedA.z0.source = {"Feed_Tank"};
sim.step(i).BedA.z0.flow_law = {"molar_flow"}; %flow is constant molar flow rate
sim.step(i).BedA.z0.parameters = {struct('Fm', 0.03)}; %molar feed rate 0.016 mol/s
% z=L: Flow OUT. Destination is Product_Tank.
sim.step(i).BedA.zL.destination = {"Product_Tank"};
sim.step(i).BedA.zL.flow_law = {"valve"}; % Flow driven by P_BedA_zL - P_Product
sim.step(i).BedA.zL.parameters = {struct('Cv', 0.1)};

% --- Bed B Configuration (State 3: u_z0 < 0, u_zL = 0) ---
sim.step(i).BedB.state = 3;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedB.z0.destination = {"Vent_Tank"};
sim.step(i).BedB.z0.flow_law = {"valve"}; % Flow driven by P_BedB_z0 - P_Vent
sim.step(i).BedB.z0.parameters = {struct('Cv', 0.1)};
% z=L: No flow.
sim.step(i).BedB.zL.flow_law = {"none"};

%----------------------------------------------------------------------
% Step 2: Bed A in State 4 (Adsorption), Bed B in State 2 (Purge)
%----------------------------------------------------------------------
i = 2;
% --- Bed A Configuration (State 4: u_z0 > 0, u_zL > 0) ---
sim.step(i).BedA.state = 4;
% z=0: Flow IN. Source is Feed_Tank.
sim.step(i).BedA.z0.source = {"Feed_Tank"};
sim.step(i).BedA.z0.flow_law = {"molar_flow"}; %flow is constant molar flow rate
sim.step(i).BedA.z0.parameters = {struct('Fm', 0.03)}; %molar feed rate 0.016 mol/s
% z=L: Flow OUT. Destination is Product_Tank AND BedB_zL.
sim.step(i).BedA.zL.destination = {"Product_Tank", "BedB_zL"};
sim.step(i).BedA.zL.flow_law = {"valve_split"}; % ONE law for total flow
sim.step(i).BedA.zL.parameters = {struct('Cv', 0.1, 'split_frac', 0.1)};
% The solver will use P_BedA_zL and P_Product to get Q_total,
% then split it.

% --- Bed B Configuration (State 2: u_z0 < 0, u_zL < 0) ---
sim.step(i).BedB.state = 2;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedB.z0.destination = {"Vent_Tank"};
sim.step(i).BedB.z0.flow_law = {"valve"}; % Flow driven by P_BedB_z0 - P_Vent
sim.step(i).BedB.z0.parameters = {struct('Cv', 0.1)};
% z=L: Flow IN. Source is Bed A's zL.
sim.step(i).BedB.zL.source = {"BedA_zL"};
sim.step(i).BedB.zL.flow_law = {"linked"}; % Flow rate set by upstream split
sim.step(i).BedB.zL.parameters = {};

%----------------------------------------------------------------------
% Step 3: Pressure Equalization - Bed A to Bed B
% Bed A: State 6 (Depressurization), Bed B: State 8 (Repressurization)
%----------------------------------------------------------------------
i = 3;
% --- Bed A Configuration (State 6: u_z0 = 0, u_zL > 0) ---
sim.step(i).BedA.state = 6;
sim.step(i).BedA.z0.flow_law = {"none"};
% z=L: Flow OUT. Destination is Bed B's z=L.
sim.step(i).BedA.zL.destination = {"BedB_zL"};
sim.step(i).BedA.zL.flow_law = {"valve"}; % Flow driven by P_A_zL > P_B_zL
sim.step(i).BedA.zL.parameters = {struct('Cv', 0.1)}; % Valve coefficient for equalization line
% The solver will automatically use P at z=L for BedA and BedB.

% --- Bed B Configuration (State 8: u_z0 = 0, u_zL < 0) ---
sim.step(i).BedB.state = 8;
sim.step(i).BedB.z0.flow_law = {"none"};
% z=L: Flow IN. Source is Bed A's zL.
sim.step(i).BedB.zL.source = {"BedA_zL"};
sim.step(i).BedB.zL.flow_law = {"linked"}; % Flow rate is determined by the upstream valve
sim.step(i).BedB.zL.parameters = {};

%----------------------------------------------------------------------
% Step 4: Bed A in State 3 (Blowdown), Bed B in State 8 (Repressurization from Product)
%----------------------------------------------------------------------
i = 4;
% --- Bed A Configuration (State 3: u_z0 < 0, u_zL = 0) ---
sim.step(i).BedA.state = 3;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedA.z0.destination = {"Vent_Tank"};
sim.step(i).BedA.z0.flow_law = {"valve"}; % Flow driven by P_BedA_z0 - P_Vent
sim.step(i).BedA.z0.parameters = {struct('Cv', 0.1)};
% z=L: No flow.
sim.step(i).BedA.zL.flow_law = {"none"};

% --- Bed B Configuration (State 8: u_z0 = 0, u_zL < 0) ---
sim.step(i).BedB.state = 8;
% z=0: No flow.
sim.step(i).BedB.z0.flow_law = {"none"};
% z=L: Flow IN. Source is Product_Tank.
sim.step(i).BedB.zL.source = {"Product_Tank"};
sim.step(i).BedB.zL.flow_law = {"valve"}; % Flow driven by P_Product - P_BedB_zL
sim.step(i).BedB.zL.parameters = {struct('Cv', 0.1)};

%----------------------------------------------------------------------
% Step 5: Bed A in State 3 (Blowdown), Bed B in State 4 (Adsorption)
%----------------------------------------------------------------------
i = 5;
% --- Bed A Configuration (State 3: u_z0 < 0, u_zL = 0) ---
sim.step(i).BedA.state = 3;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedA.z0.destination = {"Vent_Tank"};
sim.step(i).BedA.z0.flow_law = {"valve"}; % Flow driven by P_BedA_z0 - P_Vent
sim.step(i).BedA.z0.parameters = {struct('Cv', 0.1)};
% z=L: No flow.
sim.step(i).BedA.zL.flow_law = {"none"};

% --- Bed B Configuration (State 4: u_z0 > 0, u_zL > 0) ---
sim.step(i).BedB.state = 4;
% z=0: Flow IN. Source is Feed_Tank.
sim.step(i).BedB.z0.source = {"Feed_Tank"};
sim.step(i).BedB.z0.flow_law = {"molar_flow"}; % Flow driven by DP
sim.step(i).BedB.z0.parameters = {struct('Fm', 0.03)};
% z=L: Flow OUT. Destination is Product_Tank.
sim.step(i).BedB.zL.destination = {"Product_Tank"};
sim.step(i).BedB.zL.flow_law = {"valve"}; % Flow driven by P_BedB_zL - P_Product
sim.step(i).BedB.zL.parameters = {struct('Cv', 0.1)};

%----------------------------------------------------------------------
% Step 6: Bed B in State 4 (Adsorption), Bed A in State 2 (Purge)
%
%----------------------------------------------------------------------
i = 6;
% --- Bed B Configuration (State 4: u_z0 > 0, u_zL > 0) ---
sim.step(i).BedB.state = 4;
% z=0: Flow IN. Source is Feed_Tank.
sim.step(i).BedB.z0.source = {"Feed_Tank"};
sim.step(i).BedB.z0.flow_law = {"molar_flow"};
sim.step(i).BedB.z0.parameters = {struct('Fm', 0.03)};
% z=L: Flow OUT. Destination is Product_Tank AND BedA_zL.
sim.step(i).BedB.zL.destination = {"Product_Tank", "BedA_zL"};
sim.step(i).BedB.zL.flow_law = {"valve_split"}; % ONE law for total flow
sim.step(i).BedB.zL.parameters = {struct('Cv', 0.1, 'split_frac', 0.1)};
% The solver will use P_BedB_zL and P_Product to get Q_total,
% then split it. 10% of flow is diverted to purge Bed A.

% --- Bed A Configuration (State 2: u_z0 < 0, u_zL < 0) ---
sim.step(i).BedA.state = 2;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedA.z0.destination = {"Vent_Tank"};
sim.step(i).BedA.z0.flow_law = {"valve"}; % Flow driven by P_BedA_z0 - P_Vent
sim.step(i).BedA.z0.parameters = {struct('Cv', 0.1)};
% z=L: Flow IN. Source is Bed B's zL.
sim.step(i).BedA.zL.source = {"BedB_zL"};
sim.step(i).BedA.zL.flow_law = {"linked"}; % Flow rate set by upstream split
sim.step(i).BedA.zL.parameters = {};

%----------------------------------------------------------------------
% Step 7: Pressure Equalization - Bed B to Bed A
% Bed B: State 6 (Depressurization), Bed A: State 8 (Repressurization)
% Mirror image of Step 3
%----------------------------------------------------------------------
i = 7;
% --- Bed B Configuration (State 6: u_z0 = 0, u_zL > 0) ---
sim.step(i).BedB.state = 6;
% z=0: No flow.
sim.step(i).BedB.z0.flow_law = {"none"};
% z=L: Flow OUT in +Z direction. Destination is Bed A's z=L.
sim.step(i).BedB.zL.destination = {"BedA_zL"};
sim.step(i).BedB.zL.flow_law = {"valve"}; % Flow driven by P_BedB_zL > P_BedA_zL
sim.step(i).BedB.zL.parameters = {struct('Cv', 0.1)}; % Valve coefficient for equalization line

% --- Bed A Configuration (State 8: u_z0 = 0, u_zL < 0) ---
sim.step(i).BedA.state = 8;
% z=0: No flow.
sim.step(i).BedA.z0.flow_law = {"none"};
% z=L: Flow IN in -Z direction. Source is Bed B's z=L.
sim.step(i).BedA.zL.source = {"BedB_zL"};
sim.step(i).BedA.zL.flow_law = {"linked"}; % Accepts the flow from Bed B
sim.step(i).BedA.zL.parameters = {};

%----------------------------------------------------------------------
% Step 8: Bed B in State 3 (Blowdown), Bed A in State 8 (Repressurization from Product)
% Mirror image of Step 4
%----------------------------------------------------------------------
i = 8;
% --- Bed B Configuration (State 3: u_z0 < 0, u_zL = 0) ---
sim.step(i).BedB.state = 3;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedB.z0.destination = {"Vent_Tank"};
sim.step(i).BedB.z0.flow_law = {"valve"}; % Flow driven by P_BedB_z0 - P_Vent
sim.step(i).BedB.z0.parameters = {struct('Cv', 0.1)};
% z=L: No flow.
sim.step(i).BedB.zL.flow_law = {"none"};

% --- Bed A Configuration (State 8: u_z0 = 0, u_zL < 0) ---
sim.step(i).BedA.state = 8;
% z=0: No flow.
sim.step(i).BedA.z0.flow_law = {"none"};
% z=L: Flow IN. Source is Product_Tank.
sim.step(i).BedA.zL.source = {"Product_Tank"};
sim.step(i).BedA.zL.flow_law = {"valve"}; % Flow driven by P_Product - P_BedA_zL
sim.step(i).BedA.zL.parameters = {struct('Cv', 0.1)};
%% ------------------SOLVER OPTIONS-----------------------------------
sim.ode_solver = 'ode15s';      % Solver choice
sim.RelTol = 1e-3;              % Relative tolerance for ODE solver
sim.AbsTol = 1e-6;              % Absolute tolerance for ODE solver
sim.dt_max = 2.0;               % Max time step [s]
sim.dt_min = 1e-5;              % Min time step [s]

%% ------------------OUTPUT AND MONITORING----------------------------
sim.save_interval = 1.0;        % Interval for saving results [s]
sim.verbose = true;             % Print progress info to console
sim.use_mass_balance_check = true;

%% === Output Control ===
sim.output.save_profiles = true;
sim.output.plot_profiles = true;
sim.output.report_recovery = true;

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
        sim.init_conditions{i}.P0 = 5.00e5 * ones(sim.num_nodes, 1);
        sim.init_conditions{i}.T0 = 298 * ones(sim.num_nodes, 1);
        sim.init_conditions{i}.y0 = repmat([1.0, 0.00], sim.num_nodes, 1);
    end

    [sim, bed_states, tank_states] = initializeFromFresh(sim);
end

end
