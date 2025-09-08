function casper_visualize_last_cycle(matfile)
% Plot spatial profiles at the end of each step of the LAST cycle.
% - Pressure & Temperature: one figure per bed (overlaid by step)
% - Mole fractions: one figure per (bed, species) using species names
% - Loadings:      one figure per (bed, species) using species names
% - Total loading: one figure per bed

S = load(matfile);
sim     = S.sim;
steps   = S.last_cycle.steps;
nS      = sim.n_species;
nSteps  = numel(steps);
species = sim.species_names;    % cell array of names, e.g., {'O2','N2',...}

% styling
cmap       = parula(max(nSteps, 8));
linestyles = {'-','--',':','-.'};
markers    = {'o','s','^','d','v','>','<','p'};
markerEvery = max(1, round(sim.num_nodes/10));

for b = 1:sim.num_beds
    bedlabel = sprintf('Bed%c','A'+b-1);

    % ---------- Pressure & Temperature ----------
    f = figure('Name',sprintf('%s — P/T (end of step)',bedlabel),'Color','w');
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % Pressure
    nexttile; hold on; grid on; title([bedlabel ' — Pressure (end of step)']);
    for s = 1:nSteps
        pr = cv_last_unpack(sim, b, steps(s).Y_end);
        plot(sim.z_nodes, pr.P, ...
            'Color', cmap(s,:), ...
            'LineStyle', linestyles{1+mod(s-1,numel(linestyles))}, ...
            'LineWidth', 1.5, ...
            'Marker', markers{1+mod(s-1,numel(markers))}, ...
            'MarkerIndices', 1:markerEvery:numel(sim.z_nodes), ...
            'DisplayName', sprintf('Step %d', s));
    end
    xlabel('z (m)'); ylabel('P (Pa)'); legend('show','Location','bestoutside');

    % Temperature
    nexttile; hold on; grid on; title([bedlabel ' — Temperature (end of step)']);
    for s = 1:nSteps
        pr = cv_last_unpack(sim, b, steps(s).Y_end);
        plot(sim.z_nodes, pr.T, ...
            'Color', cmap(s,:), ...
            'LineStyle', linestyles{1+mod(s-1,numel(linestyles))}, ...
            'LineWidth', 1.5, ...
            'Marker', markers{1+mod(s-1,numel(markers))}, ...
            'MarkerIndices', 1:markerEvery:numel(sim.z_nodes), ...
            'DisplayName', sprintf('Step %d', s));
    end
    xlabel('z (m)'); ylabel('T (K)'); legend('show','Location','bestoutside');

    % ---------- Mole fractions: per species ----------
    for j = 1:nS
        spname = species{j};
        f = figure('Name',sprintf('%s — y(%s) (end of step)',bedlabel, spname),'Color','w'); %#ok<NASGU>
        hold on; grid on; title(sprintf('%s — Mole fraction %s (end of step)', bedlabel, spname));
        for s = 1:nSteps
            pr = cv_last_unpack(sim, b, steps(s).Y_end);
            plot(sim.z_nodes, pr.y(:,j), ...
                'Color', cmap(s,:), ...
                'LineStyle', linestyles{1+mod(s-1,numel(linestyles))}, ...
                'LineWidth', 1.5, ...
                'Marker', markers{1+mod(s-1,numel(markers))}, ...
                'MarkerIndices', 1:markerEvery:numel(sim.z_nodes), ...
                'DisplayName', sprintf('Step %d', s));
        end
        xlabel('z (m)'); ylabel(sprintf('y(%s)', spname));
        legend('show','Location','bestoutside');
    end

    % ---------- Loadings: per species ----------
    for j = 1:nS
        spname = species{j};
        f = figure('Name',sprintf('%s — Loading %s (end of step)',bedlabel, spname),'Color','w'); %#ok<NASGU>
        hold on; grid on; title(sprintf('%s — Loading %s (end of step)', bedlabel, spname));
        for s = 1:nSteps
            pr = cv_last_unpack(sim, b, steps(s).Y_end);
            plot(sim.z_nodes, pr.q(:,j), ...
                'Color', cmap(s,:), ...
                'LineStyle', linestyles{1+mod(s-1,numel(linestyles))}, ...
                'LineWidth', 1.5, ...
                'Marker', markers{1+mod(s-1,numel(markers))}, ...
                'MarkerIndices', 1:markerEvery:numel(sim.z_nodes), ...
                'DisplayName', sprintf('Step %d', s));
        end
        xlabel('z (m)'); ylabel(sprintf('q(%s)', spname));
        legend('show','Location','bestoutside');
    end

    % ---------- Total loading ----------
    f = figure('Name',sprintf('%s — q_{tot} (end of step)',bedlabel),'Color','w'); %#ok<NASGU>
    hold on; grid on; title([bedlabel ' — Total loading q_{tot} (end of step)']);
    for s = 1:nSteps
        pr = cv_last_unpack(sim, b, steps(s).Y_end);
        plot(sim.z_nodes, sum(pr.q,2), ...
            'Color', cmap(s,:), ...
            'LineStyle', linestyles{1+mod(s-1,numel(linestyles))}, ...
            'LineWidth', 1.8, ...
            'Marker', markers{1+mod(s-1,numel(markers))}, ...
            'MarkerIndices', 1:markerEvery:numel(sim.z_nodes), ...
            'DisplayName', sprintf('Step %d', s));
    end
    xlabel('z (m)'); ylabel('q_{tot}');
    legend('show','Location','bestoutside');
end
end

% ===== helper =====
function prof = cv_last_unpack(sim, b, Y)
beds = unpack_bed_state_vector(Y, sim);
bed  = beds{b};
prof.P = bed.C .* sim.R .* bed.T;
prof.T = bed.T;
prof.y = bed.y;     % [N x n_species]
prof.q = bed.q;     % [N x n_species]
end
