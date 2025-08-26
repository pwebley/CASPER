function bc = get_bed_interface_props(Y, bed_index, sim)
% GET_BED_INTERFACE_PROPS Extracts interface conditions at z=0 and z=L for a bed.
% Returns a struct 'bc' with fields: z0 and zL, each with P, T, y.

NS = sim.num_nodes;
N  = sim.n_species;
R  = sim.R;

% Calculate offset for the bed in the full state vector
n_per_bed = NS * (1 + 2*N);
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

% --- Unpack temperature ---
idx_T = offset + NS + NS*(N-1) + (1:NS);
T = Y(idx_T);

% Interface values at z = 0
bc.z0.P = Ct(1) * R * T(1);
bc.z0.T = T(1);
bc.z0.y = yi(1, :)';

% Interface values at z = L
bc.zL.P = Ct(end) * R * T(end);
bc.zL.T = T(end);
bc.zL.y = yi(end, :)';

end
