function [P_in, P_out, T_in, T_out, y_in, y_out] = get_bed_interface_props(y, bed_index, parm)
% GET_BED_INTERFACE_PROPS   Returns inlet/outlet P, T, y for a given bed
%
%   [P_in, P_out, T_in, T_out, y_in, y_out] = get_bed_interface_props(y, bed_index, parm)
%
%   Inputs:
%       y          : Full ODE state vector
%       bed_index  : 1 or 2
%       parm       : Parameter struct for the bed
%
%   Outputs:
%       P_in, P_out    : Pressure at z = 0 and z = L
%       T_in, T_out    : Temperature at z = 0 and z = L
%       y_in, y_out    : Mole fractions (vector) at z = 0 and z = L

NS = parm.NS;         % Number of spatial nodes
N  = parm.N_gas;      % Number of gas species
R  = parm.R;

% Determine offset for current bed
parm.varlen = parm.NS + parm.NS*(parm.N_gas - 1) + parm.NS + parm.NS * parm.N_gas;
offset = (bed_index - 1) * parm.varlen;

% === Total gas concentration ===
c_total = y(offset + 1 : offset + NS);  % mol/m³

% === Gas-phase concentrations (N-1 stored species) ===
c_partial = reshape(y(offset + NS + 1 : offset + NS + NS*(N-1)), [NS, N-1]);  % mol/m³

% === Compute concentration of species N ===
cN = c_total - sum(c_partial, 2);  % mol/m³

% === Corrected Full concentration matrix ===
c_gas = zeros(NS, N);
c_gas(:,1:N-1) = c_partial;     % stored species 2 to N
c_gas(:,N)   = cN;            % inferred species 1

% === Compute mole fractions ===
y_full = c_gas ./ c_total;

% Extract inlet and outlet mole fractions
y_in  = y_full(1, :).';
y_out = y_full(end, :).';

% === Extract temperature profile ===
T_vector = y(offset + NS + NS*(N-1) + 1 : offset + NS + NS*(N-1) + NS);
T_in  = T_vector(1);
T_out = T_vector(end);

% === Compute pressures using ideal gas law ===
P_in  = c_total(1)  * R * T_in;
P_out = c_total(end) * R * T_out;

end
