
function print_step_summary(cycle_idx, step_idx, step_name, sim, Y_end, tspan, bc_z0, bc_zL)
% PRINT_STEP_SUMMARY — Compact sanity check after each step
%   Prints inlet/outlet P, y(species 1), q(species 1) for each bed.

    beds = unpack_bed_state_vector(Y_end, sim);
    N = sim.n_species;
    fprintf('    ── Summary @ t=%.2f s (cycle %d, step %d: %s) ──\n', ...
            tspan(2), cycle_idx, step_idx, step_name);

    for b = 1:sim.num_beds
        bed = beds{b};
        P = bed.C .* sim.R .* bed.T;        % ideal-gas P profile
        y = bed.y;                           % NS x N
        qi = bed.q;                          % NS x N

        P_in = P(1);   P_out = P(end);
        y1_in = y(1,1); y1_out = y(end,1);
        q1_in = qi(1,1); q1_out = qi(end,1);

        % Pull simple BC velocity signs if available
        u_in = NaN;  u_out = NaN;
        if b <= numel(bc_z0) && isstruct(bc_z0{b}) && isfield(bc_z0{b},'u')
            u_in = bc_z0{b}.u;
        end
        if b <= numel(bc_zL) && isstruct(bc_zL{b}) && isfield(bc_zL{b},'u')
            u_out = bc_zL{b}.u;
        end

        in_arrow  = velocity_arrow(u_in, 'z0');
        out_arrow = velocity_arrow(u_out, 'zL');

        fprintf('      Bed %c | Pin %.3f bar %s  Pout %.3f bar %s  | y1_in %.3f  y1_out %.3f  | q1_in %.3f  q1_out %.3f\n', ...
                'A'+b-1, 1e-5*P_in, in_arrow, 1e-5*P_out, out_arrow, y1_in, y1_out, q1_in, q1_out);
    end
end

function s = velocity_arrow(u, pos)
    if isnan(u), s = ''; return; end
    if strcmp(pos,'z0')
        if u > 0, s = '→'; elseif u < 0, s = '←'; else, s = '·'; end
    else
        if u > 0, s = '→'; elseif u < 0, s = '←'; else, s = '·'; end
    end
end
