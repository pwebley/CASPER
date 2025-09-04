
function plot_diagnostics(sim)
% PLOT_DIAGNOSTICS — Minimal placeholder (extend as needed)
% Looks at the latest state saved in 'cycle_sim.mat' if present; otherwise does nothing.

    if exist('cycle_sim.mat','file')
        S = load('cycle_sim.mat');
        if isfield(S,'Y_end') && isfield(S,'sim')
            beds = unpack_bed_state_vector(S.Y_end, S.sim);
            figure; 
            for b = 1:sim.num_beds
                subplot(sim.num_beds,1,b);
                P = beds{b}.C .* S.sim.R .* beds{b}.T;
                plot(P,'LineWidth',1.5); ylabel('P (Pa)'); xlabel('node'); grid on
                title(sprintf('Pressure profile at end (Bed %c)', 'A'+b-1));
            end
        end
    end
end
