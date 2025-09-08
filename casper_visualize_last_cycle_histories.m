function casper_visualize_last_cycle_histories(matfile)
% CASPER_VISUALIZE_LAST_CYCLE_HISTORIES
% Show boundary histories over the LAST cycle and mark step boundaries.

S       = load(matfile);
sim     = S.sim;
last    = S.last_cycle;
steps   = last.steps;
species = sim.species_names;
nB      = sim.num_beds;
nS      = sim.n_species;

% ---- Gather per-bed concatenated histories ----
concat = struct();
for b = 1:nB
    bedlabel = sprintf('Bed%c','A'+b-1);
    concat.(bedlabel).z0 = init_series(nS);
    concat.(bedlabel).zL = init_series(nS);
end

% Also track adjusted step boundaries for xline() annotations
stepStart = zeros(1, numel(steps));  % adjusted absolute start times
stepEnd   = zeros(1, numel(steps));  % adjusted absolute end times

t_prev_end = -inf;  % for offsetting non-absolute step times
for s = 1:numel(steps)
    st = steps(s);
    H  = st.hist;

    % time vector: prefer hist.t, else Tseg
    if isfield(H,'t') && ~isempty(H.t)
        t = H.t(:);
    elseif isfield(st,'Tseg') && ~isempty(st.Tseg)
        t = st.Tseg(:);
    else
        error('Step %d has no time vector (neither hist.t nor Tseg).', s);
    end

    % offset if this step’s time isn’t strictly increasing vs previous
    if ~isempty(t) && t(1) <= t_prev_end
        t = t + (t_prev_end - t(1) + eps);
    end
    if ~isempty(t)
        stepStart(s) = t(1);
        stepEnd(s)   = t(end);
        t_prev_end   = t(end);
    else
        stepStart(s) = t_prev_end;
        stepEnd(s)   = t_prev_end;
    end

    % append per-bed, per-end signals (assign back)
    for b = 1:nB
        bedlabel = sprintf('Bed%c','A'+b-1);
        Hb = H.(bedlabel);
        concat.(bedlabel).z0 = append_series(concat.(bedlabel).z0, t, Hb.z0);
        concat.(bedlabel).zL = append_series(concat.(bedlabel).zL, t, Hb.zL);
    end
end

% ---- Plot per bed ----
for b = 1:nB
    bedlabel = sprintf('Bed%c','A'+b-1);

    % ---------- z=0 ----------
    H0 = concat.(bedlabel).z0;
    f = figure('Name',sprintf('%s — z0 histories (last cycle)', bedlabel),'Color','w'); %#ok<NASGU>
    tl = tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

    nexttile; plot(H0.t, H0.P,'LineWidth',1.4); grid on
    title(sprintf('%s z0 — Pressure'), 'Interpreter','none'); xlabel('t (s)'); ylabel('P (Pa)');
    add_step_markers(stepStart, tl);

    nexttile; plot(H0.t, H0.T,'LineWidth',1.4); grid on
    title(sprintf('%s z0 — Temperature'), 'Interpreter','none'); xlabel('t (s)'); ylabel('T (K)');
    add_step_markers(stepStart, tl);

    nexttile; plot(H0.t, H0.u,'LineWidth',1.4); grid on
    title(sprintf('%s z0 — Velocity'), 'Interpreter','none'); xlabel('t (s)'); ylabel('u (m/s)');
    add_step_markers(stepStart, tl);

    nexttile; plot(H0.t, H0.n,'LineWidth',1.4); grid on
    title(sprintf('%s z0 — Molar flow n'), 'Interpreter','none'); xlabel('t (s)'); ylabel('mol/s');
    add_step_markers(stepStart, tl);

    ax = nexttile([1 2]); hold on; grid on
    for j = 1:nS
        plot(H0.t, H0.y(:,j), 'LineWidth',1.4, 'DisplayName', species{j});
    end
    title(sprintf('%s z0 — Mole fractions'), 'Interpreter','none');
    xlabel('t (s)'); ylabel('y_i'); ylim([0 1]); legend('show','Location','bestoutside');
    add_step_markers(stepStart, tl, ax);

    % ---------- z=L ----------
    HL = concat.(bedlabel).zL;
    f = figure('Name',sprintf('%s — zL histories (last cycle)', bedlabel),'Color','w'); %#ok<NASGU>
    tl = tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

    nexttile; plot(HL.t, HL.P,'LineWidth',1.4); grid on
    title(sprintf('%s zL — Pressure'), 'Interpreter','none'); xlabel('t (s)'); ylabel('P (Pa)');
    add_step_markers(stepStart, tl);

    nexttile; plot(HL.t, HL.T,'LineWidth',1.4); grid on
    title(sprintf('%s zL — Temperature'), 'Interpreter','none'); xlabel('t (s)'); ylabel('T (K)');
    add_step_markers(stepStart, tl);

    nexttile; plot(HL.t, HL.u,'LineWidth',1.4); grid on
    title(sprintf('%s zL — Velocity'), 'Interpreter','none'); xlabel('t (s)'); ylabel('u (m/s)');
    add_step_markers(stepStart, tl);

    nexttile; plot(HL.t, HL.n,'LineWidth',1.4); grid on
    title(sprintf('%s zL — Molar flow n'), 'Interpreter','none'); xlabel('t (s)'); ylabel('mol/s');
    add_step_markers(stepStart, tl);

    ax = nexttile([1 2]); hold on; grid on
    for j = 1:nS
        plot(HL.t, HL.y(:,j), 'LineWidth',1.4, 'DisplayName', species{j});
    end
    title(sprintf('%s zL — Mole fractions'), 'Interpreter','none');
    xlabel('t (s)'); ylabel('y_i'); ylim([0 1]); legend('show','Location','bestoutside');
    add_step_markers(stepStart, tl, ax);
end
end

% ===== helpers =====
function S = init_series(nS)
S.t = [];
S.P = []; S.T = []; S.u = []; S.n = [];
S.y = zeros(0,nS);
end

function Sdst = append_series(Sdst, t, Ssrc_end)
Sdst.t = [Sdst.t; t(:)];
Sdst.P = [Sdst.P; Ssrc_end.P(:)];
Sdst.T = [Sdst.T; Ssrc_end.T(:)];
Sdst.u = [Sdst.u; Ssrc_end.u(:)];
Sdst.n = [Sdst.n; Ssrc_end.n(:)];
ysrc = Ssrc_end.y;
if isvector(ysrc), ysrc = ysrc(:).'; end
if size(ysrc,1) ~= numel(t)
    if size(ysrc,1)==1, ysrc = repmat(ysrc, numel(t), 1);
    else, error('Mole fraction array has incompatible size.'); end
end
Sdst.y = [Sdst.y; ysrc];
end

function add_step_markers(stepStart, tl, ax)
% Draw red dashed vertical lines at each step boundary with labels.
% Works with either current axes or provided ax.
if nargin < 3 || isempty(ax), ax = gca; end
axes(ax); %#ok<LAXES>
for s = 1:numel(stepStart)
    xline(stepStart(s), 'r--', sprintf('Step %d', s), ...
        'HandleVisibility','off', ...
        'LabelVerticalAlignment','bottom', ...
        'LabelHorizontalAlignment','center');
end
end
