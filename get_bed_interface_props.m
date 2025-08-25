function bc = get_bed_interface_props(Y, bed_index, sim)
% GET_BED_INTERFACE_PROPS Extracts interface conditions for a given bed.
% Returns a struct 'bc' with fields: P, T, y (all at inlet and outlet).

NS = sim.num_nodes;
N  = sim.n_species;
R  = sim.R;

% Calculate offset for the bed in the full state vector
n_per_bed = NS * (2 + N);
offset = (bed_index - 1) * n_per_bed;

% --- Unpack total concentration ---
idx_Ct = offset + (1:NS);
Ct = Y(idx_Ct);

% --- Unpack partial concentrations ---
idx_Ci = offset + NS + (1:NS*(N-1));
Ci_partial = reshape(Y(idx_Ci), NS, N-1);

% Infer the last species concentration
Ci = zeros(NS, N);
Ci(:,1:N-1) = Ci_partial;
Ci(:,N) = Ct - sum(Ci_partial, 2);

% Compute mole fractions
yi = Ci ./ Ct;

% Inlet/outlet fractions
y_in  = yi(1, :)';
y_out = yi(end, :)';

% --- Unpack temperature ---
idx_T = offset + NS*N + (1:NS);
T = Y(idx_T);
T_in  = T(1);
T_out = T(end);

% --- Compute pressures ---
P_in  = Ct(1)  * R * T_in;
P_out = Ct(end) * R * T_out;

% Return as structured BC
bc.inlet.P = P_in;
bc.inlet.T = T_in;
bc.inlet.y = y_in;

bc.outlet.P = P_out;
bc.outlet.T = T_out;
bc.outlet.y = y_out;
end