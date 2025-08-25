function dy = bed_ode(t, y, sim, inlet_bc, outlet_bc, bed_index)
% BED_ODE Right-hand side of ODE system for a single PSA bed using finite volume.
% 
% Inputs:
%   t         : current time (s)
%   y         : state vector [C_tot; C_i (1:N-1); T; q_i (1:N)]
%   sim       : simulation structure
%   inlet_bc  : structure with .Q, .P, .T, .y at z = 0
%   outlet_bc : structure with .Q, .P, .T, .y at z = L
%   bed_index : index of bed being computed (for logging or control)
%
% Output:
%   dy        : time derivative of state vector

% === Unpack constants ===
NS = sim.num_nodes;         % Number of spatial nodes
N  = sim.n_species;         % Number of gas species
R  = sim.R;                 % Gas constant
eps = sim.epsilon;          % Bed porosity
rho_b = sim.rho_b;          % Bulk density of bed
A = sim.A_bed;              % Cross-sectional area
dz = sim.dz(:);             % Vector of node lengths

% === Index mapping ===
id_Ct = 1:NS;
id_Ci = reshape(NS + (1:NS*(N-1)), NS, N-1);
id_T  = NS*N + (1:NS);
id_qi = reshape(NS*(N+1) + (1:NS*N), NS, N);

% === Unpack state variables ===
Ct = max(y(id_Ct), 1e-12);  % Total concentration
Ci = zeros(NS, N);
Ci(:, 1:N-1) = y(id_Ci);    % Explicit species
Ci(:, N) = Ct - sum(Ci(:, 1:N-1), 2);  % Balance last species
Ci = max(Ci, 1e-12);        % Avoid division by zero
T = max(y(id_T), 1e-6);     % Temperature
qi = y(id_qi);              % Adsorbed loadings

% === Mixture properties ===
yi = Ci ./ Ct;
P = Ct .* R .* T;
mu = sim.visc_func(T);
M = yi * sim.MW(:);                   % Mixture MW
rho_g = P .* M ./ (R * T);            % Gas density

% === Compute Ergun velocity profile ===
v = zeros(NS+1, 1);
for s = 2:NS
    dP = P(s-1) - P(s);
    v(s) = sign(dP) * sqrt(abs(dP) / (sim.Ergun_factor * rho_g(s)));
end
v(1)     = inlet_bc.Q  / (A * Ct(1));      % Inlet BC
v(NS+1)  = outlet_bc.Q / (A * Ct(end));    % Outlet BC

% === Upwind function ===
upwind = @(phi) (v >= 0).*phi(1:end-1) + (v < 0).*phi(2:end);

% === Compute equilibrium loading (per node) ===
ads_node = arrayfun(@(n) sim.layers(sim.layer_id(n)).properties, ...
                    1:NS, 'UniformOutput', false);
P_partial = yi .* P;
q_eq = get_equilibrium_loading_all(sim, vertcat(ads_node{:}), P_partial, T);
dqi_dt = sim.k_LDF .* (q_eq - qi);

% === Initialize outputs ===
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
    dCi_dt(:,i) = (1/eps) * (-(flux_i(2:end) - flux_i(1:end-1)) ./ dz ...
                    - rho_b * dqi_dt(:,i));
end

% === Total concentration balance ===
Ct_in  = inlet_bc.P / (R * inlet_bc.T);
Ct_out = outlet_bc.P / (R * outlet_bc.T);
Ct_ext = [Ct_in; Ct; Ct_out];
flux_Ct = v .* upwind(Ct_ext);
dCt_dt = (1/eps) * (-(flux_Ct(2:end) - flux_Ct(1:end-1)) ./ dz ...
            - rho_b * sum(dqi_dt, 2));

% === Energy balance ===
cp_g = sim.cp_g;
cp_s = sim.cp_s;
dH   = sim.dH;
source = sim.heat_source;
sink   = sim.heat_sink;

T_ext = [T(1); T; T(end)];
Ct_face = upwind(Ct_ext);
T_face  = upwind(T_ext);
enthalpy_flux = v .* Ct_face .* cp_g .* T_face;
convective_term = -(enthalpy_flux(2:end) - enthalpy_flux(1:end-1)) ./ dz;

reaction_heat = -rho_b * sum(dqi_dt .* dH(:)', 2);
denominator = eps .* Ct .* cp_g + (1 - eps) .* rho_b .* cp_s;

dT_dt = (convective_term - R .* T .* dCt_dt + reaction_heat + source - sink) ./ denominator;

% === Return full state derivative ===
dy = [dCt_dt;
      reshape(dCi_dt(:,1:N-1), [], 1);
      dT_dt;
      reshape(dqi_dt, [], 1)];
end
