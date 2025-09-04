function [bc_z0, bc_zL] = flow_network(t, sim, raw, current_step)
% FLOW_NETWORK Computes boundary conditions for each bed at z=0 and z=L.

num_beds = sim.num_beds;
bc_z0 = cell(num_beds, 1);
bc_zL = cell(num_beds, 1);

% === Pass 1: process destinations first (these define outflows) ===
for b = 1:num_beds
    bed_name = sprintf('Bed%c', 'A' + b - 1);
    step_config = current_step.(bed_name);

    if isfield(step_config,'z0') && isfield(step_config.z0,'destination') && ~isempty(step_config.z0.destination)
        bc_z0{b} = compute_boundary_condition_z('z0', raw, b, sim, step_config, bed_name);
    end
    if isfield(step_config,'zL') && isfield(step_config.zL,'destination') && ~isempty(step_config.zL.destination)
        bc_zL{b} = compute_boundary_condition_z('zL', raw, b, sim, step_config, bed_name);
    end
end

% === Pass 2: now process sources (may depend on linked results) ===
for b = 1:num_beds
    bed_name = sprintf('Bed%c', 'A' + b - 1);
    step_config = current_step.(bed_name);

    if isempty(bc_z0{b})
        bc_z0{b} = compute_boundary_condition_z('z0', raw, b, sim, step_config, bed_name);
    end
    if isempty(bc_zL{b})
        bc_zL{b} = compute_boundary_condition_z('zL', raw, b, sim, step_config, bed_name);
    end
end
end

% -------------------------------------------------------------------------
function bc = compute_boundary_condition_z(pos, raw, bed_id, sim, step_config, bed_name)
% Build BC for one end (z0/zL) of one bed for the current step.

interface = step_config.(pos);
props     = raw{bed_id}.(pos);

% Default BC = local bed properties, no flow
bc.P = props.P;
bc.T = props.T;
bc.y = props.y;
bc.u = 0;

% No law or 'none' => keep defaults
if ~isfield(interface,'flow_law') || isempty(interface.flow_law)
    return;
end

% flow_law normalized to scalar, lower-case string
law = interface.flow_law; if iscell(law), law = law{1}; end
law = lower(string(law));
if law == "none", return; end

% --- Normalize parameters & endpoints
params = [];
if isfield(interface,'parameters') && ~isempty(interface.parameters)
    params = interface.parameters{1};
end

src_name = [];
dst_name = [];

if isfield(interface,'source') && ~isempty(interface.source)
    tmp = interface.source;
    if iscell(tmp)
        src_name = tmp{1};
    else
        src_name = tmp;
    end
end

if isfield(interface,'destination') && ~isempty(interface.destination)
    tmp = interface.destination;
    if iscell(tmp)
        dst_name = tmp{1};
    else
        dst_name = tmp;
    end
end

sign_correction = get_flow_sign(step_config.state, pos);

switch law

    case "valve"
        if isempty(params), error('flow_network:valve:MissingParams','Valve needs parameters (Cv) at %s of %s.', pos, bed_name); end
        if ~isempty(src_name)
            [Tsrc, ysrc, Psrc] = resolve_source_properties(src_name, raw, sim);
            bc.T = Tsrc; bc.y = ysrc;
            P1 = Psrc;  P2 = props.P;  T = Tsrc; y = ysrc(:);
        elseif ~isempty(dst_name)
            [~, ~, Pdst] = resolve_source_properties(dst_name, raw, sim);
            P1 = props.P; P2 = Pdst; T = props.T; y = props.y(:);
        else
            error('Valve flow must define source or destination at %s of %s.', pos, bed_name);
        end
        Q = valve_flow(P1, P2, T, y, sim, params.Cv);      % L/s (standard)
        bc.u = convert_Ls_to_velocity(Q, sim, props);
        % NEW: if destination is a bed endpoint, stash this flow for the linked side
        if ~isempty(dst_name) && is_bed_endpoint(dst_name)
            % store with the UPSTREAM bed+pos key, i.e., this bed end
            link_store('set', bed_name, pos, Q);
        end
    case "molar_flow"
        if isempty(params), error('flow_network:molar_flow:MissingParams','Molar flow needs parameters (Fm) at %s of %s.', pos, bed_name); end
        C_local = props.P/(sim.R*props.T);
        if ~isempty(src_name)
            [Tsrc, ysrc, ~] = resolve_source_properties(src_name, raw, sim);
            bc.T = Tsrc; bc.y = ysrc;
            bc.u =  params.Fm / (C_local * sim.A_bed);
        elseif ~isempty(dst_name)
            bc.u = -params.Fm / (C_local * sim.A_bed);
        else
            error('Molar flow requires source or destination at %s of %s.', pos, bed_name);
        end
        % NEW: stash as L/s (STP) when destination is a bed endpoint
        if ~isempty(dst_name) && is_bed_endpoint(dst_name)
            % mol/s -> m^3/s(STP) -> L/s
            Q_Ls = params.Fm * (sim.R*273/101325) * 1000;
            link_store('set', bed_name, pos, Q_Ls);
        end
    case "volume_flow"
        if isempty(params), error('flow_network:volume_flow:MissingParams','Volume flow needs parameters (Q) at %s of %s.', pos, bed_name); end
        if ~isempty(src_name)
            [Tsrc, ysrc, ~] = resolve_source_properties(src_name, raw, sim);
            bc.T = Tsrc; bc.y = ysrc;
            bc.u =  params.Q / sim.A_bed;
        elseif ~isempty(dst_name)
            bc.u = -params.Q / sim.A_bed;
        else
            error('Volume flow requires source or destination at %s of %s.', pos, bed_name);
        end
        % NEW: stash as L/s (STP) when destination is a bed endpoint
        if ~isempty(dst_name) && is_bed_endpoint(dst_name)
            % Q_actual (m^3/s) -> n_dot = C_actual*Q_actual -> Vdot_STP = n_dot*R*Tstp/Pstp
            Q_Ls = params.Q * (props.P*273/(101325*props.T)) * 1000;  % L/s at STP
            link_store('set', bed_name, pos, Q_Ls);
        end
    case "valve_split"
        if isempty(params) || ~isfield(interface,'destination') || numel(interface.destination) < 2
            error('flow_network:valve_split:BadConfig','valve_split needs two destinations and parameters (Cv, split_frac) at %s of %s.', pos, bed_name);
        end
        Cv = params.Cv; split_frac = params.split_frac;
        dst1 = interface.destination{1};   % tank branch (for P2)
        [~, ~, P2_1] = resolve_source_properties(dst1, raw, sim);
        P1 = props.P; T = props.T; y = props.y(:);

        Q_total = valve_flow(P1, P2_1, T, y, sim, Cv);  % L/s (standard)
        Q_link  = split_frac * Q_total;

        bc.u = convert_Ls_to_velocity(Q_total, sim, props);
        link_store('set', bed_name, pos, Q_link);

    case "linked"
        % Expect interface.source to be a fully-qualified bed endpoint: 'BedX_z0' or 'BedX_zL'
        Q_link = 0;

        if ~isempty(src_name) && is_bed_endpoint(src_name)
            % Parse the upstream endpoint and fetch the stored flow
            [isBed, up_idx, up_pos] = try_parse_bed_endpoint(src_name, sim);
            if isBed
                up_bed_name = sprintf('Bed%c', 'A' + up_idx - 1);
                try
                    Q_link = link_store('get', up_bed_name, up_pos);   % L/s at STP
                catch
                    Q_link = 0;   % no upstream stash found -> treat as zero
                end
                % Use upstream gas composition/temperature
                [Tsrc, ysrc, ~] = resolve_source_properties(src_name, raw, sim);
                bc.T = Tsrc; bc.y = ysrc;
            end
        else
            % If source is a tank or missing, we have no upstream bed flow to link to.
            % Keep bc.T, bc.y as local and Q_link=0 (no inflow).
            % (You can add a warning here if desired.)
        end

        % Convert linked flow (standard L/s) → interstitial velocity

        bc.u = convert_Ls_to_velocity_incoming(Q_link, sim, props.P, bc.T);

    otherwise
        error('Unsupported flow_law "%s" at %s of %s.', law, pos, bed_name);
end

% Apply sign convention after velocity is computed
bc.u = sign_correction * bc.u;
end

% ---------- per-call store for linked flows ----------
function out = link_store(op, bed_name, pos, Q)
persistent S
if isempty(S), S = struct(); end
key = sprintf('%s_%s', bed_name, pos);
switch op
    case 'set', S.(key) = Q; out = [];
    case 'get'
        if isfield(S,key), out = S.(key);
        else, error('Linked flow for %s not set.', key);
        end
    otherwise, error('Invalid link_store op: %s', op);
end
end

% -------------------------------------------------------------------------
function sgn = get_flow_sign(state, pos)
switch pos
    case 'z0', sgn_map = [ 1, -1, -1,  1, -1, 0, 1,  0, 0 ];
    case 'zL', sgn_map = [-1, -1,  0,  1,  1, 1, 0, -1, 0 ];
    otherwise, error('Invalid position: %s. Must be ''z0'' or ''zL''.', pos);
end
sgn = sgn_map(state);
end

% -------------------------------------------------------------------------
function vel = convert_Ls_to_velocity(Q_Ls, sim, props)
Q_m3_s = Q_Ls / 1000;
n_dot  = Q_m3_s * 101325 / (sim.R * 273);
C      = props.P / (sim.R * props.T);
vel    = n_dot / (C * sim.A_bed);
end

% -------------------------------------------------------------------------
function tf = is_bed_endpoint(name)
% True for strings like 'BedA_z0' or 'BedB_zL' (case-insensitive)
tf = false;
if ~(ischar(name) || isstring(name)), return; end
s = char(string(name));
tf = ~isempty(regexp(s, '^Bed[A-Za-z]_z(0|L)$', 'once', 'ignorecase'));
end

function tf = is_tank_name(name, sim)
tf = false;
if ~(ischar(name) || isstring(name)), return; end
s = string(name);
for k = 1:numel(sim.tanks)
    if strcmpi(s, sim.tanks(k).name), tf = true; return; end
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



function [Tsrc, ysrc, Psrc] = resolve_source_properties(name, raw, sim)
% RESOLVE_SOURCE_PROPERTIES
% Accepts any of the following "name" formats:
%   • struct('bed_name','BedB','pos','zL')          % preferred for beds
%   • struct('bed','BedB','pos','z0')               % alias for bed_name
%   • struct('bed_idx',2,'pos','zL')                % 1-based bed index + pos
%   • {'BedB','zL'} / {'BedA','z0'}                 % cell (bed endpoint)
%   • 'BedB_zL' / "BedA_z0"                         % string (bed endpoint)
%   • 'TankName'                                    % tank name (matches sim.tanks(k).name)
% Returns gas interface properties (Tsrc [K], ysrc [N×1], Psrc [Pa]).

    % Try to parse as a bed endpoint first
    [isBed, bed_idx, pos_str] = try_parse_bed_endpoint(name, sim);

    if isBed
        props = raw{bed_idx}.(pos_str);   % pos_str is 'z0' or 'zL'
        Tsrc  = props.T;
        ysrc  = props.y(:);
        Psrc  = props.P;
        return
    end

    % Otherwise treat as a tank reference (string/cell/struct with .name)
    tank = try_get_tank(name, sim);
    if ~isempty(tank)
        Tsrc = tank.T;
        ysrc = tank.y(:);
        Psrc = tank.P;
        return
    end

    % Nothing matched
    error('resolve_source_properties:Unresolved', ...
          'Could not resolve source "%s" as a bed endpoint or tank.', printable_name(name));
end

%==================== helpers ====================%

function [isBed, bed_idx, pos_str] = try_parse_bed_endpoint(name, sim)
% Returns true if "name" refers to a Bed endpoint; also returns bed_idx (1-based) and pos_str ('z0'/'zL').

    isBed   = false; bed_idx = []; pos_str = '';

    % ---- Struct forms ----
    if isstruct(name)
        % Accept multiple field conventions
        if all(isfield(name, {'bed_idx','pos'}))
            bed_idx = name.bed_idx;
            pos_str = lower(string(name.pos));
            validate_bed_and_pos(bed_idx, pos_str, sim);
            isBed = true; return
        end
        if all(isfield(name, {'bed_name','pos'})) || all(isfield(name, {'bed','pos'}))
            if isfield(name,'bed_name'), btoken = string(name.bed_name); else, btoken = string(name.bed); end
            bed_idx = bed_token_to_index(btoken);
            pos_str = lower(string(name.pos));
            validate_bed_and_pos(bed_idx, pos_str, sim);
            isBed = true; return
        end
        % If a struct with only .name provided, we’ll try it as a tank later
    end

    % ---- Cell forms ----
    if iscell(name)
        if numel(name) == 2
            btoken = string(name{1});
            pos_str = lower(string(name{2}));
            if startsWith(btoken, "Bed", 'IgnoreCase', true)
                bed_idx = bed_token_to_index(btoken);
                validate_bed_and_pos(bed_idx, pos_str, sim);
                isBed = true; return
            else
                % likely {'TankName', <ignored>} — let tank handler try it
            end
        elseif numel(name) == 1 && (ischar(name{1}) || isstring(name{1}))
            % Let downstream string/tank parser try it
            name = name{1};
        else
            return
        end
    end

    % ---- String forms ----
    if ischar(name) || isstring(name)
        s = string(name);

        % Pattern "BedX_z0" or "BedX_zL" (case-insensitive)
        % Examples: "BedA_z0", 'BedC_zL'
        expr = "^Bed([A-Za-z])_z([0L])$";
        tok = regexp(s, expr, 'tokens', 'once', 'ignorecase');
        if ~isempty(tok)
            bchar   = upper(string(tok{1}));           % 'A','B',...
            poscode = upper(string(tok{2}));           % '0' or 'L'
            bed_idx = double(char(bchar)) - double('A') + 1;
            pos_str = char("z"+poscode);               % 'z0' or 'zL'
            validate_bed_and_pos(bed_idx, pos_str, sim);
            isBed = true; return
        end

        % Plain "BedX" without pos is ambiguous — do NOT accept (caller must supply pos)
        if startsWith(s, "Bed", 'IgnoreCase', true)
            % Not a valid endpoint
            isBed = false; return
        end
    end
end

function validate_bed_and_pos(bed_idx, pos_str, sim)
    if ~(isscalar(bed_idx) && isnumeric(bed_idx) && bed_idx>=1 && bed_idx<=sim.num_beds)
        error('resolve_source_properties:BadBedIndex', ...
              'Invalid bed index %s. There are %d beds.', mat2str(bed_idx), sim.num_beds);
    end
    pos_str = lower(string(pos_str));
    if ~(pos_str=="z0" || pos_str=="zl")
        error('resolve_source_properties:BadPos', ...
              'Invalid position "%s". Expected "z0" or "zL".', pos_str);
    end
end

function idx = bed_token_to_index(btoken)
% btoken like "BedA" / 'BedB'
    s = char(btoken);
    if numel(s)~=4 || ~strncmpi(s,'Bed',3)
        error('resolve_source_properties:BadBedToken', ...
              'Expected "BedX" (e.g., "BedB"), got "%s".', s);
    end
    idx = double(upper(s(4))) - double('A') + 1;
end

function tank = try_get_tank(name, sim)
% Returns tank struct if name matches a tank; otherwise [].
    tank = [];
    % struct with .name
    if isstruct(name) && isfield(name,'name')
        s = string(name.name);
        for k = 1:numel(sim.tanks)
            if strcmpi(s, sim.tanks(k).name), tank = sim.tanks(k); return, end
        end
        return
    end
    % cell {'TankName'}
    if iscell(name) && numel(name)==1 && (ischar(name{1}) || isstring(name{1}))
        name = name{1};
    end
    % string / char
    if ischar(name) || isstring(name)
        s = string(name);
        for k = 1:numel(sim.tanks)
            if strcmpi(s, sim.tanks(k).name), tank = sim.tanks(k); return, end
        end
    end
end

function s = printable_name(name)
% Produce a readable token for error messages
    if ischar(name) || isstring(name)
        s = char(string(name));
    elseif iscell(name)
        try
            s = strjoin(string(name), ', ');
        catch
            s = '<cell>';
        end
    elseif isstruct(name)
        f = fieldnames(name);
        s = sprintf('struct{%s}', strjoin(f, ','));
    else
        s = class(name);
    end
end

function vel = convert_Ls_to_velocity_incoming(Q_Ls, sim, P_ref, T_in)
% Convert STP L/s to interstitial velocity using receiving boundary pressure (P_ref)
% and incoming stream temperature (T_in). Composition y does not affect total C for ideal gas.
    Q_m3_s = Q_Ls / 1000;
    n_dot  = Q_m3_s * 101325 / (sim.R * 273);      % mol/s from STP volumetric flow
    C_loc  = P_ref / (sim.R * T_in);               % mol/m^3 at boundary using incoming T
    vel    = n_dot / (C_loc * sim.A_bed);
end
