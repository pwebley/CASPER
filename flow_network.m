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
    %% Case Valve: 
    % We model flow through a valve using the Flow Controls Institute equation. We model 
    % this as a check valve so reverse flow is not permitted 
    case 'valve'  
        if isfield(interface, 'source')
            % --- Flow from source to bed ---
            src_name = interface.source{1};
            [Tsrc, ysrc, Psrc] = resolve_source_properties(src_name, raw, sim);
            bc.T = Tsrc;
            bc.y = ysrc;

            yi = ysrc(:);  % Use source composition
            T_abs = Tsrc;  % Use source temperature
            P1 = Psrc / 101325;      % Source (upstream) pressure [atm]
            P2 = props.P / 101325;   % Bed (downstream) pressure [atm]

        elseif isfield(interface, 'destination')
            % --- Flow from bed to destination ---
            dst_name = interface.destination{1};
            [~, ~, Pdst] = resolve_source_properties(dst_name, raw, sim);

            yi = props.y(:);       % Use bed interface composition
            T_abs = props.T;       % Use bed temperature
            P1 = props.P / 101325; % Bed (upstream) pressure [atm]
            P2 = Pdst / 101325;    % Destination (downstream) pressure [atm]

            % bc.T and bc.y remain as bed properties (already assigned)
        else
            error('Valve flow must define either a source or destination at %s', pos);
        end

        % === Valve flow calculation ===
        MW_species = sim.MW(:) / 1000;      % [kg/mol]
        MW_mix = yi' * MW_species;          % Mixture MW [kg/mol]
        G = MW_mix / 0.02897;               % Specific gravity

        if P2 >= 0.53 * P1
            Q = 77.01 * params.Cv * sqrt((P1^2 - P2^2) / (G * T_abs));
        else
            Q = 77.01 * params.Cv * (P1 / sqrt(G * T_abs));
        end

        % Convert standard L/s to mol/s at STP (273 K, 101325 Pa)
        Q_m3_s = Q / 1000;
        n_dot = Q_m3_s * 101325 / (sim.R * 273);

        % Convert mol/s to superficial velocity (m/s)
        C = props.P / (sim.R * props.T);  % Use bed gas concentration
        bc.u = n_dot / (C * sim.A_bed);
%% Case Molar: This is the case where the user specifies a fixed molar flow rate (mol/s)

    case 'molar_flow'
        C_local=props.P/(sim.R*props.T); 
        if isfield(interface, 'source')
            src_name = interface.source{1};
            [Tsrc, ysrc, ~] = resolve_source_properties(src_name, raw, sim);
            C = Psrc / (sim.R * Tsrc);
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
        dst_name = interface.destination{1};
        [~, ~, Pdst] = resolve_source_properties(dst_name, raw, sim);
        dP = props.P - Pdst;
        Q_total = params.Cv * dP;
        bc.u = (1 - params.split_frac) * Q_total;

    case 'linked'
        bc.u = NaN;

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
