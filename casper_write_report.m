function casper_write_report(sim, logger, out_filename)
% CASPER_WRITE_REPORT
% Writes a rich text report for the *last simulated cycle* plus inputs.
%
% Includes:
%   - Input summary (gases, adsorbents, bed geometry/kinetics)
%   - Step-by-step boundary moles at z=0 and z=L (total + species) per bed
%   - Step-by-step tank moles in/out (total + species)
%   - Per-step mass balances per bed (gas + adsorbed storage change)
%   - End-of-step internal node snapshots (node #1 and node #end): T, P, y
%
% Usage:
%   casper_write_report(sim, logger, 'path/to/report.txt')
%
% Expects 'logger.last_cycle.steps(k).step_flow', '...Y_end', and
% 'logger.last_cycle.Y_start_cycle' (all produced by your current driver).

if nargin < 3 || isempty(out_filename)
    ts = logger.run_info.timestamp;
    out_filename = sprintf('casper_report_%s__%s.txt', logger.run_info.run_name, ts);
end

fid = fopen(out_filename,'w');
if fid==-1, error('Could not open report file "%s"', out_filename); end
c = onCleanup(@() fclose(fid));

% ---------------- Header & run info ----------------
fprintf(fid, '====================================================\n');
fprintf(fid, 'CASPER Simulation Report (rich)\n');
fprintf(fid, '====================================================\n\n');
fprintf(fid, 'Run name        : %s\n', logger.run_info.run_name);
fprintf(fid, 'Timestamp       : %s\n', logger.run_info.timestamp);
fprintf(fid, 'Saved run file  : %s\n\n', logger.run_info.saved_as);

% ---------------- Input summary ----------------
fprintf(fid, '--- INPUT SUMMARY ---\n');
fprintf(fid, 'Simulation Name : %s\n', getf(sim,'name','(unnamed)'));
fprintf(fid, 'Species         : %s\n', strjoin(sim.species_names, ', '));
fprintf(fid, 'Beds            : %d\n', sim.num_beds);
fprintf(fid, 'Cycle time (s)  : %.6g\n', getf(sim,'cycle_time',NaN));
fprintf(fid, 'Bed diameter (m): %.4g\n', getf(sim,'bed_diameter',NaN));
fprintf(fid, 'Total length (m): %.4g\n', getf(sim,'bed_length',NaN));
fprintf(fid, 'Total nodes     : %d\n', getf(sim,'num_nodes',NaN));
fprintf(fid, 'Product tank    : %s\n', getf(sim,'product_tank_name','(unset)'));
fprintf(fid, 'Feed tank       : %s\n', getf(sim,'feed_tank_name','(unset)'));
fprintf(fid, '\n');

% Gas props
if isfield(sim,'gas')
    fprintf(fid, 'Gases (MW, Cp):\n');
    for j=1:sim.n_species
        gw = NaN; cp = NaN;
        if j<=numel(sim.gas) && isstruct(sim.gas(j))
            if isfield(sim.gas(j),'MW'), gw = sim.gas(j).MW; end
            if isfield(sim.gas(j),'Cp'), cp = sim.gas(j).Cp; end
        end
        fprintf(fid, '  %-6s  MW=%.6g kg/mol  Cp=%.6g J/mol/K\n', sim.species_names{j}, gw, cp);
    end
    fprintf(fid, '\n');
end

% Layers / adsorbents
if isfield(sim,'layers')
    fprintf(fid, 'Layers / Adsorbents:\n');
    for L=1:sim.n_layers
        props = sim.layers(L).properties;
        fprintf(fid, '  Layer %d: %s  L=%.3g m  N=%d  nodes[%d..%d]\n', ...
            L, getf(sim.layers(L),'adsorbent_name','(unknown)'), ...
            sim.layers(L).length, sim.layers(L).num_nodes, ...
            sim.layers(L).node_start, sim.layers(L).node_end);
        write_if_field(fid, '    void_fraction', props, 'void_fraction', '%.4f');
        write_if_field(fid, '    bed_density (kg/m^3)', props, 'bed_density', '%.4g');
        write_if_field(fid, '    heat_capacity (J/kg/K)', props, 'heat_capacity', '%.4g');
        write_if_field(fid, '    particle_diameter (m)', props, 'particle_diameter', '%.4g');
        fprintf(fid, '\n');
    end
end

% ---------------- Last cycle summaries ----------------
if ~isfield(logger,'last_cycle') || ~isfield(logger.last_cycle,'steps') || isempty(logger.last_cycle.steps)
    fprintf(fid, 'No last cycle step histories present.\n');
    return
end

steps = logger.last_cycle.steps;
nSteps = numel(steps);
nB = sim.num_beds; nS = sim.n_species; R = sim.R;

fprintf(fid, '--- LAST CYCLE STEP-BY-STEP ---\n');
fprintf(fid, 'Steps: %d\n\n', nSteps);

% helper lambda for per-bed storage (gas + ads) by species
function ns = bed_storage_species(bed)
    % Gas phase species moles
    % bed.C_species: [NS x nS] mol/m^3; sim.node_volume: [NS x 1] m^3
    Ngas = (bed.C_species.' * sim.node_volume(:));  % 1 x nS (mol)
    % Adsorbed species moles
    % bed.q: [NS x nS] mol/kg; sim.adsorbent_mass: [NS x 1] kg
    Nads = (bed.q.' * sim.adsorbent_mass(:));       % 1 x nS (mol)
    ns = Ngas(:) + Nads(:);                          % nS x 1
end

% Iterate steps
Y_prev = logger.last_cycle.Y_start_cycle;   % packed state at start of step 1
for s = 1:nSteps
    st = steps(s);

    fprintf(fid, 'Step %d\n', s);
    fprintf(fid, '---------\n');

    % -------- boundary flows per bed --------
    fprintf(fid, 'Boundary moles (signed integrals over step):\n');
    fprintf(fid, '  %-6s  %-2s  %14s    %s\n', 'Bed','End','Total (mol)','Species (mol)');
    for b = 1:nB
        btag = sprintf('Bed%c','A'+b-1);
        bz0T = st.step_flow.beds(b).Fz0_total;
        bzLT = st.step_flow.beds(b).FzL_total;
        bz0S = st.step_flow.beds(b).Fz0_species(:);
        bzLS = st.step_flow.beds(b).FzL_species(:);

        fprintf(fid, '  %-6s  %-2s  %14.6e    ', btag, 'z0', bz0T); print_vec(fid, bz0S);
        fprintf(fid, '  %-6s  %-2s  %14.6e    ', btag, 'zL', bzLT); print_vec(fid, bzLS);
    end
    fprintf(fid, '\n');

    % -------- tank totals this step --------
    if isfield(st.step_flow,'tanks')
        fprintf(fid, 'Tank moles (positive = into tank):\n');
        fprintf(fid, '  %-14s  %14s    %s\n', 'Tank', 'Total (mol)', 'Species (mol)');
        for k = 1:numel(st.step_flow.tanks)
            fprintf(fid, '  %-14s  %14.6e    ', st.step_flow.tanks(k).name, st.step_flow.tanks(k).n_total);
            print_vec(fid, st.step_flow.tanks(k).n_species(:));
        end
        fprintf(fid, '\n');
    end

    % -------- mass balance per bed for this step --------
    % unpack start and end bed states
    if s==1
        beds_start = unpack_bed_state_vector(Y_prev, sim);
    else
        beds_start = unpack_bed_state_vector(steps(s-1).Y_end, sim);
    end
    beds_end   = unpack_bed_state_vector(st.Y_end, sim);

    fprintf(fid, 'Per-bed mass balance (species): Net_in - ΔStorage [should ≈ 0]\n');
    for b = 1:nB
        ns_start = bed_storage_species(beds_start{b});   % nS x 1
        ns_end   = bed_storage_species(beds_end{b});     % nS x 1
        dStore   = ns_end - ns_start;                    % nS x 1

        % Net into bed using signed convention: s0 - sL
        s0 = st.step_flow.beds(b).Fz0_species(:);
        sL = st.step_flow.beds(b).FzL_species(:);
        net_in = s0 - sL;

        resid = net_in - dStore;
        fprintf(fid, '  Bed %c: ', 'A'+b-1);
        fprintf(fid, '|| ΔStore:'); print_vec(fid, dStore);
        fprintf(fid, '     NetIn:'); print_vec(fid, net_in);
        fprintf(fid, '    Resid.:'); print_vec(fid, resid);
    end
    fprintf(fid, '\n');

    % -------- end-of-step internal node snapshots (node 1 and node end) --------
    fprintf(fid, 'End-of-step internal nodes (not faces): node#1 and node#end\n');
    for b = 1:nB
        bedE = beds_end{b};
        P_profile = bedE.C .* R .* bedE.T;   % Pa
        y_profile = bedE.y;                  % NS x nS
        n1 = 1; nN = size(bedE.C,1);

        fprintf(fid, '  Bed %c:\n', 'A'+b-1);
        fprintf(fid, '    Node 1:  T=%.4g K   P=%.6g Pa   y=', bedE.T(n1), P_profile(n1)); print_vec(fid, y_profile(n1,:)');
        fprintf(fid, '    Node N:  T=%.4g K   P=%.6g Pa   y=', bedE.T(nN), P_profile(nN)); print_vec(fid, y_profile(nN,:)');
    end
    fprintf(fid, '\n');

    % prepare for next
    Y_prev = st.Y_end;
end

% ---------------- Cycle totals (from cycle_summaries) ----------------
if isfield(logger,'cycle_summaries') && ~isempty(logger.cycle_summaries)
    cs = logger.cycle_summaries(end);  % last cycle
    fprintf(fid, '--- LAST CYCLE TOTALS ---\n');
    fprintf(fid, 'Beds (signed totals):\n');
    fprintf(fid, '  %-6s  %-2s  %14s    %s\n', 'Bed','End','Total (mol)','Species (mol)');
    for b = 1:nB
        fprintf(fid, '  %-6s  %-2s  %14.6e    ', sprintf('Bed%c','A'+b-1), 'z0', cs.bed_cycle(b).Fz0_total);
        print_vec(fid, cs.bed_cycle(b).Fz0_species(:));
        fprintf(fid, '  %-6s  %-2s  %14.6e    ', sprintf('Bed%c','A'+b-1), 'zL', cs.bed_cycle(b).FzL_total);
        print_vec(fid, cs.bed_cycle(b).FzL_species(:));
    end
    fprintf(fid, '\n  Tanks (positive = into tank):\n');
    fprintf(fid, '  %-14s  %14s    %s\n', 'Tank','Total (mol)','Species (mol)');
    for k=1:numel(cs.tank_cycle)
        fprintf(fid, '  %-14s  %14.6e    ', cs.tank_cycle(k).name, cs.tank_cycle(k).n_total);
        print_vec(fid, cs.tank_cycle(k).n_species(:));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '====================================================\n');
fprintf(fid, 'End of report\n');
fprintf(fid, '====================================================\n');

fprintf('✅ Wrote report to %s\n', out_filename);

% ----- helpers -----
function v = getf(S,fld,default)
    if isstruct(S) && isfield(S,fld), v = S.(fld); else, v = default; end
end
end % casper_write_report

function write_if_field(fid, label, S, fld, fmt)
if isfield(S,fld)
    fprintf(fid, '%s: ', label);
    fprintf(fid, [fmt '\n'], S.(fld));
end
end

function print_vec(fid, v)
v = v(:);
fprintf(fid, '[');
for j=1:numel(v)
    if j>1, fprintf(fid, ' '); end
    fprintf(fid, '%.6e', v(j));
end
fprintf(fid, ']\n');
end
