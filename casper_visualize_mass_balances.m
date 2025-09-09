function casper_visualize_mass_balances(matfile)
% Compute per-cycle mass-balance errors (%) for each bed, overall and per species,
% and plot convergence toward cyclic steady state. Works with:
%   • S.cycle_summaries.beds   (preferred)
%   • S.cycles(c).steps        (recompute)
%   • S.last_cycle.steps       (recompute last cycle only)

S   = load(matfile);
sim = S.sim;
species = sim.species_names; nS = sim.n_species; nB = sim.num_beds;

% ---- Source of bed totals per cycle ----
if isfield(S, 'cycle_summaries') && isstruct(S.cycle_summaries) ...
        && isfield(S.cycle_summaries, 'beds')
    % Already summarized
    C = S.cycle_summaries;
elseif isfield(S, 'cycles') && ~isempty(S.cycles)
    % Recompute for every cycle
    C = recompute_all_cycles(sim, S.cycles);
elseif isfield(S, 'last_cycle') && isfield(S.last_cycle, 'steps')
    % Recompute only for the last cycle; warn
    warning('Only last_cycle found in %s. Computing errors for the last cycle only.', matfile);
    C = struct('beds', recompute_cycle_bed_totals(sim, S.last_cycle.steps));
else
    error('No suitable data found in %s (need cycle_summaries.beds or cycles/last_cycle with steps).', matfile);
end

% Some files store a scalar struct (single cycle). Normalize to array.
if ~isstruct(C) || isempty(C)
    error('Unable to extract cycle summaries from %s.', matfile);
end
if ~isfield(C, 'beds')
    % Defensive: if cycle_summaries exists but lacks .beds, fall back to recompute (if possible)
    if isfield(S, 'cycles') && ~isempty(S.cycles)
        C = recompute_all_cycles(sim, S.cycles);
    elseif isfield(S, 'last_cycle') && isfield(S.last_cycle,'steps')
        warning('cycle_summaries missing .beds; using last_cycle only.');
        C = struct('beds', recompute_cycle_bed_totals(sim, S.last_cycle.steps));
    else
        error('Cycle summaries exist but no .beds field and no steps to recompute.');
    end
end

% ---- Compute errors ----
nCycles = numel(C);
err_total = nan(nCycles, nB);
err_spec  = nan(nCycles, nB, nS);

fprintf('\n=== Mass-balance errors by cycle (%% of inflow) ===\n');
for cyc = 1:nCycles
    B = C(cyc).beds; % array over beds with fields: Fz0_total, FzL_total, Fz0_species, FzL_species
    fprintf('\nCycle %d\n', cyc);
    fprintf('  %-6s %-6s %12s   %s\n', 'Bed', 'Type', 'Error %', 'Per-species Error % (first 6)');

    for b = 1:nB
        [e_tot, e_spec] = bed_errors_from_signed(B(b));
        err_total(cyc,b)  = e_tot;
        err_spec(cyc,b,:) = e_spec;

        bedlabel = sprintf('Bed%c','A'+b-1);
        fprintf('  %-6s %-6s %12.3f   ', bedlabel, 'total', e_tot);
        print_vec(e_spec, min(6,nS));
    end
end

% ---- Plots: convergence across cycles ----
for b = 1:nB
    bedlabel = sprintf('Bed%c','A'+b-1);

    % Overall
    figure('Color','w','Name',[bedlabel ' — Overall balance error (% vs inflow)']);
    plot(1:nCycles, err_total(:,b), 'o-','LineWidth',1.6); grid on
    xlabel('Cycle'); ylabel('Error (%)');
    title([bedlabel ' — Overall mass-balance error (% of inflow)']);

    % Per-species
    figure('Color','w','Name',[bedlabel ' — Per-species balance error (% vs inflow)']);
    hold on; grid on
    for j = 1:nS
        plot(1:nCycles, squeeze(err_spec(:,b,j)), 'o-','LineWidth',1.4, 'DisplayName', species{j});
    end
    xlabel('Cycle'); ylabel('Error (%)');
    title([bedlabel ' — Per-species mass-balance error (% of inflow)']);
    legend('show','Location','best');
end
end

%================== helpers ==================%

function [err_total_pct, err_spec_pct] = bed_errors_from_signed(Bb)
% Signed convention:
%   z0:  + is INTO the bed
%   zL:  + is OUT of the bed (i.e., along +z)
Fz0  = Bb.Fz0_total;
FzL  = Bb.FzL_total;
S0   = Bb.Fz0_species(:);
SL   = Bb.FzL_species(:);

% Inflow to bed = into at z0 (positive) + into at zL (negative along +z)
in_total  = max(Fz0, 0) + max(-FzL, 0);
out_total = max(-Fz0, 0) + max( FzL, 0);

in_spec   = max(S0,0) + max(-SL,0);
out_spec  = max(-S0,0) + max( SL,0);

err_total_pct = 100 * (in_total - out_total) / max(in_total, eps);
den = in_spec; den(den==0) = eps;
err_spec_pct  = 100 * (in_spec - out_spec) ./ den;
end

function print_vec(v, mmax)
v = v(:); m = min(numel(v), mmax);
for k = 1:m, fprintf('%9.3f', v(k)); end
if numel(v) > m, fprintf('  ...'); end
fprintf('\n');
end

% --------- recompute from stored cycles/steps when summaries missing ---------

function C = recompute_all_cycles(sim, cycles)
C = struct([]);
for c = 1:numel(cycles)
    C(c).beds = recompute_cycle_bed_totals(sim, cycles(c).steps);
end
end

function beds_acc = recompute_cycle_bed_totals(sim, steps)
% Recompute signed bed totals over a cycle (sum of steps), like your integrator.
nB = sim.num_beds; nS = sim.n_species;
beds_acc = init_bed_accum(sim);

for s = 1:numel(steps)
    st   = steps(s);
    if ~isfield(st,'Tseg') || ~isfield(st,'Yseg') || isempty(st.Tseg)
        % If a step has no stored segment (e.g., logging off), skip it
        continue;
    end
    Tseg = st.Tseg(:);
    Yseg = st.Yseg;

    prev.has = false;
    for k = 1:numel(Tseg)
        t  = Tseg(k);
        Yk = Yseg(k,:).';

        % Interface bed properties
        raw = cell(nB,1);
        for b = 1:nB
            raw{b} = get_bed_interface_props(Yk, b, sim);
        end
        % Boundary conditions / velocities
        [bc0, bcL] = flow_network(t, sim, raw, st.step);

        % Instantaneous signed molar rates at each end (using receiving boundary P,T)
        for b = 1:nB
            C0 = raw{b}.z0.P/(sim.R*raw{b}.z0.T);
            CL = raw{b}.zL.P/(sim.R*raw{b}.zL.T);
            n0 = C0 * bc0{b}.u * sim.A_bed;          % mol/s (+ along +z)
            nL = CL * bcL{b}.u * sim.A_bed;          % mol/s (+ along +z)
            s0 = bc0{b}.y(:) * n0;                   % species mol/s
            sL = bcL{b}.y(:) * nL;

            inst(b).n0 = n0; inst(b).nL = nL;
            inst(b).s0 = s0; inst(b).sL = sL;
        end

        % Trapezoid accumulate
        if prev.has
            dt = Tseg(k) - Tseg(k-1);
            w  = 0.5*dt;
            for b = 1:nB
                beds_acc(b).Fz0_total   = beds_acc(b).Fz0_total   + w*(prev.inst(b).n0 + inst(b).n0);
                beds_acc(b).FzL_total   = beds_acc(b).FzL_total   + w*(prev.inst(b).nL + inst(b).nL);
                beds_acc(b).Fz0_species = beds_acc(b).Fz0_species + w*(prev.inst(b).s0 + inst(b).s0);
                beds_acc(b).FzL_species = beds_acc(b).FzL_species + w*(prev.inst(b).sL + inst(b).sL);
            end
        end
        prev.inst = inst; prev.has = true;
    end
end
end

function S = init_bed_accum(sim)
S = struct();
for b = 1:sim.num_beds
    S(b).Fz0_total   = 0;
    S(b).FzL_total   = 0;
    S(b).Fz0_species = zeros(sim.n_species,1);
    S(b).FzL_species = zeros(sim.n_species,1);
end
end
