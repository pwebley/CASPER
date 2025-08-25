function dYdt = multi_bed_ode(t, Y, sim, current_step)
% MULTI_BED_ODE Compute the time derivatives for all beds simultaneously.
%
% Inputs:
%   t            – current time
%   Y            – global state vector stacking all beds
%   sim          – simulation configuration struct
%   current_step – cycle step configuration struct
%
% Output:
%   dYdt         – derivative of global state vector

num_beds = sim.num_beds;
nodes = sim.num_nodes;
n_species = sim.n_species;

% Determine state length per bed
states_per_bed = nodes*(2 + n_species);

% Prepare raw interface properties for each bed
raw = cell(num_beds, 1);
for b = 1:num_beds
    raw{b} = get_bed_interface_props(Y, b, sim);
end

% Compute inlet and outlet BCs including flow rates
[inlet_bc, outlet_bc] = flow_network_update(t, sim, raw, current_step);

% Initialize derivative vector
dYdt = zeros(size(Y));

% Loop across beds
for b = 1:num_beds
    idx_start = (b-1)*states_per_bed + 1;
    idx_end   = b*states_per_bed;
    y_bed = Y(idx_start:idx_end);

    % Compute derivatives for this bed
    dy = bed_ode(t, y_bed, sim, inlet_bc{b}, outlet_bc{b}, b);

    dYdt(idx_start:idx_end) = dy;
end

end
