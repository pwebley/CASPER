function plot_diagnostics(sim, Y)
% Plot diagnostics directly from the current state vector Y
% Shows per-bed: P, T, y_j, and q_j (plus q_total).

    beds     = unpack_bed_state_vector(Y, sim);
    nBeds    = sim.num_beds;
    nSpecies = sim.n_species;
    N        = sim.num_nodes;

    figure;
    for b = 1:nBeds
        bed = beds{b};
        P   = bed.C .* sim.R .* bed.T;   % ideal-gas pressure [Pa]
        T   = bed.T;                     % temperature [K]
        y   = bed.y;                     % mole fractions [N x nSpecies]
        q   = bed.q;                     % loadings [N x nSpecies] (units per your sim)
        qtot = sum(q, 2);

        % --- Pressure ---
        subplot(nBeds,4,(b-1)*4+1);
        plot(1:N, P, 'LineWidth', 1.5);
        ylabel('P (Pa)'); xlabel('node'); grid on
        ylim([1e3, 5e5]);                       % fixed axis so tiny ripples don’t dominate
        title(sprintf('Pressure (Bed %c)', 'A'+b-1));

        % --- Temperature ---
        subplot(nBeds,4,(b-1)*4+2);
        plot(1:N, T, 'LineWidth', 1.5);
        ylabel('T (K)'); xlabel('node'); grid on
        title(sprintf('Temperature (Bed %c)', 'A'+b-1));

        % --- Mole fractions ---
        subplot(nBeds,4,(b-1)*4+3); hold on
        for j = 1:nSpecies
            plot(1:N, y(:,j), 'LineWidth', 1.5, 'DisplayName', sprintf('y_%d', j));
        end
        hold off
        ylabel('Mole fraction'); xlabel('node'); grid on
        ylim([0, 1]);
        title(sprintf('Mole fractions (Bed %c)', 'A'+b-1));
        if nSpecies <= 6, legend('show'); end

        % --- Loadings ---
        subplot(nBeds,4,(b-1)*4+4); hold on
        % total loading
        plot(1:N, qtot, 'LineWidth', 1.8, 'DisplayName', 'q_{tot}');
        % per-species loadings (only show legend if reasonable)
        for j = 1:nSpecies
            plot(1:N, q(:,j), 'LineWidth', 1.1, 'DisplayName', sprintf('q_%d', j));
        end
        hold off
        ylabel('q (loading)'); xlabel('node'); grid on
        title(sprintf('Loadings (Bed %c)', 'A'+b-1));
        if nSpecies <= 6, legend('show'); end
    end
end
