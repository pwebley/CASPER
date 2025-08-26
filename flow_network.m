function [bc_z0, bc_zL] = flow_network(t, sim, raw, current_step)
% FLOW_NETWORK Computes boundary conditions for each bed at z=0 and z=L.
%
% Inputs:
%   t             – current time
%   sim           – simulation configuration struct
%   raw           – cell array of bed interface properties from get_bed_interface_props
%   current_step  – struct with current step configuration for each bed
%
% Outputs:
%   bc_z0         – cell array of structs with boundary conditions at z=0
%   bc_zL         – cell array of structs with boundary conditions at z=L

num_beds = sim.num_beds;
bc_z0 = cell(num_beds, 1);
bc_zL = cell(num_beds, 1);

for b = 1:num_beds
    bed_name = sprintf('Bed%c', 'A' + b - 1);
    step_config = current_step.(bed_name);

    bc_z0{b} = compute_boundary_condition_z('z0', raw, b, sim, step_config);
    bc_zL{b} = compute_boundary_condition_z('zL', raw, b, sim, step_config);
end
end

function bc = compute_boundary_condition_z(pos, raw, bed_id, sim, step_config)
% Helper function to compute boundary condition at z0 or zL for a given bed.

interface = step_config.(pos);
props = raw{bed_id}.(pos);

% Default boundary condition
bc.P = props.P;
bc.T = props.T;
bc.y = props.y;
bc.u = 0;

% Only apply if flow_law is defined
if isfield(interface, 'flow_law') && ~isempty(interface.flow_law)
    law = lower(interface.flow_law{1});
    params = interface.parameters{1};

    switch law
        case 'none'
            bc.u = 0;

        case 'valve'
            if strcmp(pos, 'z0') && isfield(interface, 'source')
                tank_name = interface.source{1};
                tank = getTank(sim, tank_name);
                dP = tank.P - props.P;
                bc.u = params.Cv * dP;
                bc.y = tank.y(:);
                bc.T = tank.T;
            elseif strcmp(pos, 'zL') && isfield(interface, 'destination')
                tank_name = interface.destination{1};
                tank = getTank(sim, tank_name);
                dP = props.P - tank.P;
                bc.u = params.Cv * dP;
            end

        case 'molar_flow'
            bc.u = (params.Fm * sim.R * props.T) / (props.P * sim.A_bed);  

        case 'volume_flow'
            bc.u = params.Q / sim.A_bed;

        case 'valve_split'
            tank = getTank(sim, interface.destination{1});
            dP = props.P - tank.P;
            total_flow = params.Cv * dP;
            bc.u = (1 - params.split_frac) * total_flow;

        case 'linked'
            bc.u = NaN;  % This is resolved by upstream boundary condition

        otherwise
            error('Unsupported flow law: %s', law);
    end
end
end

function tank = getTank(sim, name)
% GETTANK Retrieves tank struct from simulation config by name.
for k = 1:length(sim.tanks)
    if strcmpi(sim.tanks(k).name, name)
        tank = sim.tanks(k);
        return;
    end
end
error('Tank "%s" not found.', name);
end
