function casper_list_last_cycle_step_flows(matfile)
% List per-step inflows/outflows (mol) for the LAST cycle using boundary histories.
% Conventions:
%   z0:  +ṅ is INTO the bed;  zL: +ṅ is OUT of the bed (along +z).
% We print positive inflow and outflow amounts (mol), integrated over each step.

S   = load(matfile);
sim = S.sim;

% ---- Find last-cycle histories ----
if isfield(S,'last_cycle') && isfield(S.last_cycle,'hist')
    Hlast = S.last_cycle.hist;
    if isfield(S.last_cycle,'step_times') && ~isempty(S.last_cycle.step_times)
        t_edges = S.last_cycle.step_times(:);
    else
        % Fall back to sim step_times relative to cycle start
        t_edges = sim.step_times(:) - sim.step_times(1);
    end
else
    error('No last_cycle.hist present in %s. Re-run driver with logging enabled.', matfile);
end

nB = sim.num_beds; nS = sim.n_species;

fprintf('\n=== Per-step inflow/outflow (mol) for LAST cycle ===\n');

for s = 1:numel(t_edges)-1
    t0 = t_edges(s);
    t1 = t_edges(s+1);

    % Accumulators for this step
    StepBeds = init_step_bed_accum(nB, nS);

    for b = 1:nB
        bedlabel = sprintf('Bed%c','A'+b-1);
        hb = Hlast.(bedlabel);

        % z0 and zL time series
        z0 = hb.z0;  zL = hb.zL;

        % Integrate this step window using trapezoid on positive inflow/outflow partitions
        [in0,out0,in0_s,out0_s] = integrate_end(z0, t0, t1, sim, +1);  % +1: +ṅ is inflow at z0
        [inL,outL,inL_s,outL_s] = integrate_end(zL, t0, t1, sim, -1);  % -1: +ṅ is outflow at zL

        StepBeds(b).z0_in    = StepBeds(b).z0_in    + in0;
        StepBeds(b).z0_out   = StepBeds(b).z0_out   + out0;
        StepBeds(b).z0_in_s  = StepBeds(b).z0_in_s  + in0_s;
        StepBeds(b).z0_out_s = StepBeds(b).z0_out_s + out0_s;

        StepBeds(b).zL_in    = StepBeds(b).zL_in    + inL;
        StepBeds(b).zL_out   = StepBeds(b).zL_out   + outL;
        StepBeds(b).zL_in_s  = StepBeds(b).zL_in_s  + inL_s;
        StepBeds(b).zL_out_s = StepBeds(b).zL_out_s + outL_s;
    end

    % -------- Print step table --------
    fprintf('\nStep %d  (t = %.2f → %.2f s)\n', s, t0, t1);
    fprintf('  %-6s %-2s  %12s  %12s   %s\n', 'Bed','End','Inflow (mol)','Outflow (mol)','Per-species (first 6):  In  |  Out');

    for b = 1:nB
        bedlabel = sprintf('Bed%c','A'+b-1);
        % z0 line
        fprintf('  %-6s %-2s  %12.6e  %12.6e   ', bedlabel,'z0', StepBeds(b).z0_in, StepBeds(b).z0_out);
        print_species_pair(StepBeds(b).z0_in_s, StepBeds(b).z0_out_s, min(6,nS));
        % zL line
        fprintf('  %-6s %-2s  %12.6e  %12.6e   ', bedlabel,'zL', StepBeds(b).zL_in, StepBeds(b).zL_out);
        print_species_pair(StepBeds(b).zL_in_s, StepBeds(b).zL_out_s, min(6,nS));
        fprintf('\n');
    end
end
end

% ================= helpers =================
function SB = init_step_bed_accum(nB, nS)
SB = struct([]);
for b = 1:nB
    SB(b).z0_in    = 0; SB(b).z0_out    = 0;
    SB(b).zL_in    = 0; SB(b).zL_out    = 0;
    SB(b).z0_in_s  = zeros(nS,1); SB(b).z0_out_s  = zeros(nS,1);
    SB(b).zL_in_s  = zeros(nS,1); SB(b).zL_out_s  = zeros(nS,1);
end
end

function [inMol,outMol,inSpec,outSpec] = integrate_end(h, t0, t1, sim, inflowSign)
% h has fields t, P, T, y, u for one end.
% inflowSign: +1 means +ṅ is inflow (z0), -1 means +ṅ is outflow (zL).
% We integrate positive inflow and positive outflow (mol) separately.

% Slice the window
[tt, P, T, y, u] = slice_hist_window(h, t0, t1);

if isempty(tt)
    inMol=0; outMol=0;
    inSpec=zeros(sim.n_species,1); outSpec=zeros(sim.n_species,1);
    return
end

% Ensure y is time-by-species
y = normalize_y_shape(y, numel(tt), sim.n_species);

% Instantaneous ṅ using receiving boundary quantities
C  = P ./ (sim.R .* T);            % mol/m^3
nd = C .* u * sim.A_bed;           % mol/s, signed along +z

% Partition into inflow/outflow relative to bed
% If inflowSign = +1 (z0): inflow = max(nd,0), outflow = max(-nd,0)
% If inflowSign = -1 (zL): inflow = max(-nd,0), outflow = max(nd,0)
if inflowSign > 0
    n_in  = max(nd,0);
    n_out = max(-nd,0);
else
    n_in  = max(-nd,0);
    n_out = max(nd,0);
end

% Species rates: ṡ = y .* ṅ  (broadcast per species)
s_in  = y .* n_in;
s_out = y .* n_out;

% Integrate with trapezoid
inMol  = trapz(tt, n_in);
outMol = trapz(tt, n_out);

inSpec  = trapz(tt, s_in, 1).';   % column vector
outSpec = trapz(tt, s_out,1).';
end

function [tt,P,T,y,u] = slice_hist_window(h, t0, t1)
tt = []; P=[]; T=[]; y=[]; u=[];
if ~all(isfield(h,{'t','P','T','y','u'})), return; end

t = h.t(:);
mask = (t >= t0) & (t <= t1);
if ~any(mask)
    return
end

tt = t(mask);
P  = h.P(mask);
T  = h.T(mask);
u  = h.u(mask);

y  = h.y(mask,:);          % if stored time-by-species, this works;
if size(y,1) ~= numel(tt) && size(y,2) == numel(tt)
    y = h.y(:,mask).';     % handle species-by-time storage
end
end

function Y = normalize_y_shape(Y, nT, nS)
% Make Y size [nT x nS]
if isequal(size(Y), [nT, nS])
    return
elseif isequal(size(Y), [nS, nT])
    Y = Y.'; return
else
    % try reshape if it's a vector with nS entries replicated
    if isvector(Y) && numel(Y)==nS
        Y = repmat(Y(:).', nT, 1);
    elseif isvector(Y) && numel(Y)==nT
        Y = repmat(Y(:), 1, nS);
    else
        error('Cannot normalize y shape: got %s, expected [%d x %d] or [%d x %d].', ...
            mat2str(size(Y)), nT, nS, nS, nT);
    end
end
end

function print_species_pair(vin, vout, mmax)
vin = vin(:); vout = vout(:);
m = min([numel(vin), numel(vout), mmax]);
for j = 1:m
    fprintf('%9.3e|%9.3e ', vin(j), vout(j));
end
if numel(vin) > m
    fprintf('...');
end
fprintf('\n');
end
