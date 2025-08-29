function [bc_z0, bc_zL] = flow_network(t, sim, raw, current_step)
% FLOW_NETWORK Computes boundary conditions for each bed at z=0 and z=L.

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

% -------------------------------------------------------------------------
function bc = compute_boundary_condition_z(pos, raw, bed_id, sim, step_config)
interface = step_config.(pos);
props = raw{bed_id}.(pos);

% Default: use local bed values
bc.P = props.P;
bc.T = props.T;
bc.y = props.y;
bc.u = 0;

if ~isfield(interface, 'flow_law') || isempty(interface.flow_law)
    return;
end

law = lower(interface.flow_law{1});
if strcmp(law, 'none')
    return;
end

params = interface.parameters{1};
bed_state = step_config.state;
sign_correction = get_flow_sign(bed_state, pos);

switch law
case 'valve'  
    if isfield(interface, 'source')
        % --- Flow from source to bed ---
        src_name = interface.source{1};
        [Tsrc, ysrc, Psrc] = resolve_source_properties(src_name, raw, sim);
        bc.T = Tsrc;
        bc.y = ysrc;

        P1 = Psrc;
        P2 = props.P;
        T = Tsrc;
        y = ysrc(:);

    elseif isfield(interface, 'destination')
        % --- Flow from bed to destination ---
        dst_name = interface.destination{1};
        [~, ~, Pdst] = resolve_source_properties(dst_name, raw, sim);

        P1 = props.P;
        P2 = Pdst;
        T = props.T;
        y = props.y(:);
        % bc.T and bc.y remain unchanged (taken from the bed)
    else
        error('Valve flow must define either a source or destination at %s', pos);
    end

    % Compute valve flow in standard L/s
    Q = valve_flow(P1, P2, T, y, sim, params.Cv);

    % Convert standard L/s to mol/s at STP (273 K, 101325 Pa)
    Q_m3_s = Q / 1000;
    n_dot = Q_m3_s * 101325 / (sim.R * 273);

    % Convert to velocity using local concentration at bed
    C = props.P / (sim.R * props.T);  % mol/m³
    bc.u = n_dot / (C * sim.A_bed);


    case 'molar_flow'
        C_local=props.P/(sim.R*props.T); 
        if isfield(interface, 'source')
            src_name = interface.source{1};
            [Tsrc, ysrc, ~] = resolve_source_properties(src_name, raw, sim);
            bc.T = Tsrc;
            bc.y = ysrc;
            bc.u = params.Fm / (C_local * sim.A_bed);
        elseif isfield(interface, 'destination')
            bc.u = -params.Fm / (C_local * sim.A_bed);
        else
            error('Molar flow requires either a source or destination at %s', pos);
        end


    case 'volume_flow'
        if isfield(interface, 'source')
            src_name = interface.source{1};
            [Tsrc, ysrc, ~] = resolve_source_properties(src_name, raw, sim);
            bc.T = Tsrc;
            bc.y = ysrc;
            bc.u = params.Q / sim.A_bed;
        elseif isfield(interface, 'destination')
            bc.u = -params.Q / sim.A_bed;
        else
            error('Volume flow requires either a source or destination at %s', pos);
        end

        
    case 'valve_split'
        % Extract split parameters
        params = interface.parameters{1};  % this has fields Cv and split_frac
        split_frac = params.split_frac;
        Cv = params.Cv;

        % Destination of split (first is tank, second is linked bed interface)
        dst_list = interface.destination;

        % Flow to tank (`destination{1}`)
        dst1 = dst_list{1};
        [~, y1, P2_1] = resolve_source_properties(dst1, raw, sim);  % ignore T
        P2_1 = P2_1;  % pressure for tank branch

        % Total flow calculation using bed properties as upstream (P1)
        P1 = props.P;
        T = props.T;
        y = props.y(:);

        Q_total = valve_flow(P1, P2_1, T, y, sim, Cv);

        % Split proportions
        Q_tank = (1 - split_frac) * Q_total;
        Q_link = split_frac * Q_total;

        % Convert Q_tank → velocity, assign to bc.u
        bc.u = convert_Ls_to_velocity(Q_tank, sim, props);

        % Temporarily store the linked branch flow in sim for the next bed to retrieve
        sim.linked_flow.(bed_name).pos = pos;  % 'zL'
        sim.linked_flow.(bed_name).flow = Q_link;

    case 'linked'
        % Retrieve upstream linked flow
        upstream = step_config.linked_from;
        Q_link = sim.linked_flow.(upstream.bed_name).flow;

        % Convert that to a velocity based on local bed properties
        bc.u = convert_Ls_to_velocity(Q_link, sim, props);

    otherwise
        error('Unsupported flow law: %s', law);
end

% Apply sign convention after velocity is computed
bc.u = sign_correction * bc.u;
end

% -------------------------------------------------------------------------
function [Tsrc, ysrc, Psrc] = resolve_source_properties(name, raw, sim)
if contains(name, 'Bed')
    bed_char = name(4);
    pos_str = name(end-1:end);  % 'z0' or 'zL'
    bed_idx = double(bed_char) - double('A') + 1;
    props = raw{bed_idx}.(pos_str);
    Tsrc = props.T;
    ysrc = props.y;
    Psrc = props.P;
else
    tank = getTank(sim, name);
    Tsrc = tank.T;
    ysrc = tank.y(:);
    Psrc = tank.P;
end
end

% -------------------------------------------------------------------------
function tank = getTank(sim, name)
for k = 1:length(sim.tanks)
    if strcmpi(sim.tanks(k).name, name)
        tank = sim.tanks(k);
        return;
    end
end
error('Tank "%s" not found.', name);
end

% -------------------------------------------------------------------------
function sgn = get_flow_sign(state, pos)
switch pos
    case 'z0'
        sgn_map = [ 1, -1, -1,  1, -1, 0, 1,  0, 0 ]; % u@z=0
    case 'zL'
        sgn_map = [-1, -1,  0,  1,  1, 1, 0, -1, 0 ]; % u@z=L
    otherwise
        error('Invalid position: %s. Must be ''z0'' or ''zL''.', pos);
end
sgn = sgn_map(state);
end

% -------------------------------------------------------------------------
function vel = convert_Ls_to_velocity(Q_Ls, sim, props)
    Q_m3_s = Q_Ls / 1000;
    n_dot = Q_m3_s * 101325 / (sim.R * 273);
    C = props.P / (sim.R * props.T);
    vel = n_dot / (C * sim.A_bed);
end
