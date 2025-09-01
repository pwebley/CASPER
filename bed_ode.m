function dy = bed_ode(t, y, sim, inlet_bc, outlet_bc, bed_index)
% BED_ODE Right-hand side of ODE system for a single PSA bed using finite volume.
% 
% Inputs:
%   t         : current time (s)
%   y         : state vector [C_tot; C_i (1:N-1); T; q_i (1:N)]
%   sim       : simulation structure
%   inlet_bc  : structure with .u, .P, .T, .y at z = 0
%   outlet_bc : structure with .u, .P, .T, .y at z = L
%   bed_index : index of bed being computed (for logging or control)
%
% Output:
%   dy        : time derivative of state vector

%% === Unpack constants ===
NS = sim.num_nodes;         % Number of spatial nodes
N  = sim.n_species;         % Number of gas species
R  = sim.R;                 % Gas constant: 8.3144 J/mol.K
eps = sim.epsilon;          % Node-wise bed porosity
rho_b = sim.rho_bed;        % Node-wise bulk density, kg/m3
A = sim.A_bed;              % Constant cross-sectional area, m2
dz = sim.dz(:);             % Local node spacing (depends on z), m

%% === Index mapping ===
id_Ct = 1:NS;
id_Ci = reshape(NS + (1:NS*(N-1)), NS, N-1);
id_T  = NS*N + (1:NS);
id_qi = reshape(NS*(N+1) + (1:NS*N), NS, N);

%% === Unpack state variables ===
Ct = max(y(id_Ct), 1e-12);  % Total concentration, kept positive
Ci = zeros(NS, N);          % Initialize Ci
Ci(:, 1:N-1) = y(id_Ci);    % Explicit species 1 to N-1
Ci(:, N) = Ct - sum(Ci(:, 1:N-1), 2);  % Balance last species
Ci = max(Ci, 1e-12);        % Species concentration, kept positive
T = max(y(id_T), 1e-6);     % Temperature, kept positive
qi = y(id_qi);              % Adsorbed loadings

%% === Mixture properties ===
yi = Ci ./ Ct;              % Vector of mole fractions, N columns
P = Ct .* R .* T;           % Vector of pressure, Pa
mu = sim.visc_func(T);      % Vector of viscosity, Pa.s
M = yi * sim.MW(:)/1000.0;  % Vector of mixture MW, kg/mol
rho_g = Ct.* M ;            % Vector of Gas density, kg/m3
cp_g = yi * sim.cp_g(:);    % Vector of gas heat capacity, J/mol.K

%% === Compute Ergun velocity profile ===
v = zeros(NS+1, 1);
for s = 2:NS
    % Average properties at face between node s-1 and s
    mu_face     = 0.5 * (mu(s-1) + mu(s));
    rho_face    = 0.5 * (rho_g(s-1) + rho_g(s));
    dz_face     = 0.5 * (dz(s-1) + dz(s));

    dP = P(s-1) - P(s);  % Pressure drop across face: positive means flow in +ve direction
    

    % Compute Ergun coefficients
    A = sim.Ergun_A_face(s) * mu_face ;
    B = sim.Ergun_B_face(s) * rho_face ;

    RHS = abs(dP) / dz_face;

    % Solve quadratic for v: A*v + B*v^2 = RHS
    v_mag = (-A + sqrt(A^2 + 4 * B * RHS)) / (2 * B);

    % Assign sign based on pressure direction
    v(s) = sign(dP) * v_mag;
end

% Inlet/outlet velocities from boundary conditions
v(1)     = inlet_bc.u;
v(NS+1)  = outlet_bc.u;

%% === Upwind function ===
upwind = @(phi) (v >= 0).*phi(1:end-1) + (v < 0).*phi(2:end);

%% === Compute equilibrium loading (per node) and rate of adsorption ===
q_eq = zeros(NS, N);
P_partial = yi .* P;

for n = 1:NS
    layer_id = sim.layer_id(n);
    ads = sim.layers(layer_id).properties;
    q_eq(n, :) = get_equilibrium_loading(sim, ads, P_partial(n, :), T(n));
end

dqi_dt = sim.k_LDF .* (q_eq - qi);

%% === Initialize outputs ===
dCt_dt = zeros(NS,1);
dCi_dt = zeros(NS,N);
dT_dt  = zeros(NS,1);

% === Species balances ===
for i = 1:N-1
    % Boundary values
    Ci_in  = inlet_bc.y(i) * inlet_bc.P / (R * inlet_bc.T);
    Ci_out = outlet_bc.y(i) * outlet_bc.P / (R * outlet_bc.T);
    
    % Extended vector with ghost nodes
    Ci_ext = [Ci_in; Ci(:,i); Ci_out];
    flux_i = v .* upwind(Ci_ext);
    
    % Mass balance
    dCi_dt(:,i) = (1./eps) .* (-(flux_i(2:end) - flux_i(1:end-1)) ./ dz ...
                    - rho_b .* dqi_dt(:,i));
end

% === Total concentration balance ===
Ct_in  = inlet_bc.P / (R * inlet_bc.T);
Ct_out = outlet_bc.P / (R * outlet_bc.T);
Ct_ext = [Ct_in; Ct; Ct_out];
flux_Ct = v .* upwind(Ct_ext);
dCt_dt = (1./eps) .* (-(flux_Ct(2:end) - flux_Ct(1:end-1)) ./ dz ...
            - rho_b .* sum(dqi_dt, 2));

% === Energy balance ===
cp_s = sim.Cp_solid;
dH   = sim.dH;
source = sim.heat_source;
sink   = sim.heat_sink;
inlet_cp_g  = sum(inlet_bc.y(:) .* sim.cp_g(:));   % scalar [J/mol/K]
outlet_cp_g = sum(outlet_bc.y(:) .* sim.cp_g(:));  % scalar [J/mol/K]
cp_g_ext = [inlet_cp_g; cp_g; outlet_cp_g];
T_ext = [inlet_bc.T; T; outlet_bc.T];

Ct_face = upwind(Ct_ext);
T_face  = upwind(T_ext);
cp_g_face = upwind(cp_g_ext);
enthalpy_flux = v .* Ct_face .* cp_g_face .* T_face;
convective_term = -(enthalpy_flux(2:end) - enthalpy_flux(1:end-1)) ./ dz;

reaction_heat = rho_b .* sum(dqi_dt .* dH, 2);
denominator = eps .* Ct .* cp_g + (1 - eps) .* rho_b .* sim.Cp_solid;

dT_dt = (convective_term - R .* T .* dCt_dt + reaction_heat + source - sink) ./ denominator;

% === Return full state derivative ===
dy = [dCt_dt;
      reshape(dCi_dt(:,1:N-1), [], 1);
      dT_dt;
      reshape(dqi_dt, [], 1)];
end
