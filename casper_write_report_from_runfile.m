function casper_write_report_from_runfile(run_mat_file, out_filename, K_last_cycles)
% Generate a CASPER report from a saved logger .mat WITHOUT re-running sims.
% Usage:
%   casper_write_report_from_runfile('casper_run_run__20250908_135717.mat');
%   casper_write_report_from_runfile(file, 'my_report.txt', 3);

if nargin < 1 || isempty(run_mat_file)
    run_mat_file = casper_choose_runfile('runs','casper_run_*.mat');
end
if nargin < 2, out_filename = ''; end
if nargin < 3 || isempty(K_last_cycles), K_last_cycles = 1; end

S = load(run_mat_file);  % expects fields from casper_logger_save: run_info, sim, cycle_summaries, last_cycle, etc.

% --- sanity
req = {'sim','cycle_summaries','run_info'};
for r = 1:numel(req)
    if ~isfield(S, req{r}), error('Run file missing "%s".', req{r}); end
end
sim   = S.sim;
runs  = S.run_info;
summ  = S.cycle_summaries;

if isempty(summ), error('No cycle_summaries found in run file.'); end

% --- choose cycles to average (last K)
all_idx = [summ.cycle_index];
[~,ord] = sort(all_idx,'ascend');
summ = summ(ord);
take = max(1, numel(summ)-K_last_cycles+1) : numel(summ);
pick = summ(take);

% --- resolve tank and species indices
prod_idx = find(strcmpi({sim.tanks.name}, sim.product_tank_name), 1);
if isempty(prod_idx), error('Product tank "%s" not found.', sim.product_tank_name); end
feed_idx = find(strcmpi({sim.tanks.name}, sim.feed_tank_name), 1);
if isempty(feed_idx), error('Feed tank "%s" not found.', sim.feed_tank_name); end

nS  = sim.n_species;
spn = sim.species_names(:);

% --- accumulate over chosen cycles
prod_total  = 0;                 % total mol to Product tank
prod_sp     = zeros(nS,1);       % species mol to Product tank
feed_total  = 0;                 % signed (will be negative = out of tank)
feed_sp     = zeros(nS,1);

for i = 1:numel(pick)
    tc = pick(i).tank_cycle;                 % struct array over all tanks
    % product tank
    prod_total = prod_total + tc(prod_idx).n_total;
    prod_sp    = prod_sp    + tc(prod_idx).n_species(:);
    % feed tank (NOTE: out of feed tank is NEGATIVE in your accounting)
    feed_total = feed_total + tc(feed_idx).n_total;
    feed_sp    = feed_sp    + tc(feed_idx).n_species(:);
end

K = numel(pick);
cycle_time = sim.cycle_time;  % seconds per cycle (from your sim setup)

% --- metrics
prod_total_pos = prod_total;                    % already positive into Product tank by construction
prod_sp_pos    = max(prod_sp, 0);               % guard (should already be positive)

feed_in_sp     = -feed_sp;                      % make positive “into the process” from Feed tank
feed_in_total  = -feed_total;                   % same

avg_prod_flow_mol_s = prod_total_pos / (K * cycle_time);    % mol/s averaged over K cycles
purity_vec           = prod_sp_pos / max(prod_total_pos, eps);

% species headline (optional)
headline_idx = [];
if isfield(sim,'primary_product_species') && ~isempty(sim.primary_product_species)
    headline_idx = find(strcmpi(spn, sim.primary_product_species), 1);
end
if isempty(headline_idx), headline_idx = 1; end  % default to species 1 if not specified

% component recovery (%): product of species i / feed of species i
recovery_vec = zeros(nS,1);
for j = 1:nS
    denom = max(feed_in_sp(j), eps);
    recovery_vec(j) = 100 * (prod_sp_pos(j) / denom);
end

% --- output file name
if isempty(out_filename)
    out_filename = sprintf('%s__report.txt', erase(runs.saved_as, '.mat'));
end
fid = fopen(out_filename,'w');
if fid == -1, error('Could not open report file "%s"', out_filename); end

% === WRITE ===
fprintf(fid, '====================================================\n');
fprintf(fid, 'CASPER Simulation Report\n');
fprintf(fid, '====================================================\n\n');
fprintf(fid, 'Run file        : %s\n', run_mat_file);
fprintf(fid, 'Saved as        : %s\n', runs.saved_as);
fprintf(fid, 'Run name        : %s\n', runs.run_name);
fprintf(fid, 'Timestamp       : %s\n', runs.timestamp);
fprintf(fid, '\n');

fprintf(fid, '--- Simulation Summary ---\n');
fprintf(fid, 'Simulation Name : %s\n', getfield_def(sim,'name','(unnamed)'));
fprintf(fid, 'Species         : %s\n', strjoin(spn.', ', '));
fprintf(fid, 'Beds            : %d\n', sim.num_beds);
fprintf(fid, 'Cycle time (s)  : %.4g\n', cycle_time);
fprintf(fid, 'Averaged cycles : last %d cycle(s)\n', K);
fprintf(fid, '\n');

fprintf(fid, '--- Tanks ---\n');
fprintf(fid, 'Feed tank       : %s\n', sim.feed_tank_name);
fprintf(fid, 'Product tank    : %s\n', sim.product_tank_name);
fprintf(fid, '\n');

fprintf(fid, '--- Product Totals (over selected cycles) ---\n');
fprintf(fid, 'Total to Product (mol) : %.6e\n', prod_total_pos);
fprintf(fid, 'Avg Product Flow (mol/s): %.6e\n', avg_prod_flow_mol_s);
fprintf(fid, 'Species to Product (mol):\n');
for j = 1:nS
    fprintf(fid, '  %-6s : %.6e\n', spn{j}, prod_sp_pos(j));
end
fprintf(fid, '\n');

fprintf(fid, '--- Product Purity (%%) ---\n');
for j = 1:nS
    fprintf(fid, '  %-6s : %7.3f %%\n', spn{j}, 100*purity_vec(j));
end
fprintf(fid, '\n');

fprintf(fid, '--- Component Recoveries (%% of feed-in) ---\n');
for j = 1:nS
    fprintf(fid, '  %-6s : %7.3f %%\n', spn{j}, recovery_vec(j));
end
fprintf(fid, '\n');

fprintf(fid, 'Headline species: %s\n', spn{headline_idx});
fprintf(fid, '  Purity   : %7.3f %%\n', 100*purity_vec(headline_idx));
fprintf(fid, '  Recovery : %7.3f %%\n', recovery_vec(headline_idx));
fprintf(fid, '\n');

% optional: bed/step diagnostic snippet if present
if isfield(S, 'last_cycle') && isfield(S.last_cycle, 'steps') && ~isempty(S.last_cycle.steps)
    fprintf(fid, '--- Last Cycle: Steps present = %d (decimated histories)\n', numel(S.last_cycle.steps));
end

fclose(fid);
fprintf('✅ Report written: %s\n', out_filename);
end

function val = getfield_def(S, field, default)
if isstruct(S) && isfield(S, field), val = S.(field); else, val = default; end
end
