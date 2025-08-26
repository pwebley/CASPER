% Load existing database

Zeolite5A_1 = struct();

Zeolite5A_1.name = 'Zeolite5A_1';
Zeolite5A_1.isotherm_type = 'dsl';  % Match existing field
Zeolite5A_1.void_fraction = 0.35;   % Estimate or update
Zeolite5A_1.particle_diameter = 1.6e-3;  % Example diameter in meters

% DSL parameters moved under 'isotherms'
Zeolite5A_1.isotherms = struct();
Zeolite5A_1.isotherms.N2 = struct( ...
    'm1', 0.69, ...
    'b01', 1.2305 / 101325, ...
    'dH1', 5461.1 * 4.184, ...
    'm2', 3.34, ...
    'b02', 0.0740 / 101325, ...
    'dH2', 4265.6 * 4.184 ...
);

Zeolite5A_1.isotherms.O2 = struct( ...
    'm1', 0.69, ...
    'b01', 0.0518 / 101325, ...
    'dH1', 3716.7 * 4.184, ...
    'm2', 3.34, ...
    'b02', 0.0518 / 101325, ...
    'dH2', 3716.7 * 4.184 ...
);

% Optional physical properties (placeholders, update with real data if known)
Zeolite5A_1.heat_capacity = 1000;        % J/kg·K
Zeolite5A_1.thermal_conductivity = 0.15; % W/m·K
Zeolite5A_1.bed_density = 650;          % kg/m³

load('adsorbent_database.mat', 'adsorbent_database');

% Ensure uniform struct array
if iscell(adsorbent_database)
    adsorbent_database = [adsorbent_database{:}];
end

fields = fieldnames(adsorbent_database(1));
for f = 1:numel(fields)
    if ~isfield(Zeolite5A_1, fields{f})
        Zeolite5A_1.(fields{f}) = [];
    end
end
Zeolite5A_1 = orderfields(Zeolite5A_1, adsorbent_database(1));

adsorbent_database(end+1) = Zeolite5A_1;
save('adsorbent_database.mat', 'adsorbent_database');
fprintf('Updated adsorbent_database.mat\n');
