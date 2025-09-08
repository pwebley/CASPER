function casper_visualize_mass_balances(matfile)
% Print per-step mass balances for last cycle and the cycle totals.

S = load(matfile);
sim = S.sim;
steps = S.last_cycle.steps;

fprintf('\n=== Mass balances (last cycle) ===\n');
for s = 1:numel(steps)
    fprintf('Step %d:\n', s);
    local_print_step_flow_report_signed(sim, s, steps(s).step_flow);
end

fprintf('\n=== Cycle totals ===\n');
cs = S.cycle_summaries(end); % last cycle summary
local_print_cycle_totals(sim, cs.cycle_index, cs.bed_cycle, cs.tank_cycle);

end

% --- local print helpers (standalone so you don't need the driver) ---

function local_print_step_flow_report_signed(sim, step_idx, step_flow)
    fprintf('    — Step %d signed boundary flows (mol) —\n', step_idx);
    fprintf('      %-6s  %-2s  %14s    %s\n', 'Bed', 'End', 'Total (mol)', 'Species (mol) [first 6]');

    for b = 1:sim.num_beds
        bedlabel = sprintf('Bed%c','A'+b-1);
        sp = step_flow.beds(b).Fz0_species;
        fprintf('      %-6s  %-2s  %14.6e    ', bedlabel, 'z0', step_flow.beds(b).Fz0_total);
        local_print_species_vec(sp, 6);

        sp = step_flow.beds(b).FzL_species;
        fprintf('      %-6s  %-2s  %14.6e    ', bedlabel, 'zL', step_flow.beds(b).FzL_total);
        local_print_species_vec(sp, 6);
        fprintf('\n');
    end

    fprintf('      \n');
    fprintf('      %-14s  %14s    %s\n', 'Tank', 'Total (mol)', 'Species (mol) [first 6]');
    for k = 1:numel(step_flow.tanks)
        fprintf('      %-14s  %14.6e    ', step_flow.tanks(k).name, step_flow.tanks(k).n_total);
        local_print_species_vec(step_flow.tanks(k).n_species, 6);
    end
end

function local_print_cycle_totals(sim, cycle_idx, bed_cycle, tank_cycle)
    fprintf('\n=== Cycle %d totals (signed) ===\n', cycle_idx);
    fprintf('  %-6s  %-2s  %14s    %s\n', 'Bed', 'End', 'Total (mol)', 'Species (mol) [first 6]');

    for b = 1:sim.num_beds
        bedlabel = sprintf('Bed%c','A'+b-1);
        sp = bed_cycle(b).Fz0_species;
        fprintf('  %-6s  %-2s  %14.6e    ', bedlabel, 'z0', bed_cycle(b).Fz0_total);
        local_print_species_vec(sp, 6);
        sp = bed_cycle(b).FzL_species;
        fprintf('  %-6s  %-2s  %14.6e    ', bedlabel, 'zL', bed_cycle(b).FzL_total);
        local_print_species_vec(sp, 6);
        fprintf('\n');
    end

    fprintf('  \n');
    fprintf('  %-14s  %14s    %s\n', 'Tank', 'Total (mol)', 'Species (mol) [first 6]');
    for k = 1:numel(tank_cycle)
        fprintf('  %-14s  %14.6e    ', tank_cycle(k).name, tank_cycle(k).n_total);
        local_print_species_vec(tank_cycle(k).n_species, 6);
    end
end

function local_print_species_vec(v, maxshow)
    v = v(:);
    m = min(numel(v), maxshow);
    for j = 1:m
        fprintf('%10.3e', v(j));
    end
    if numel(v) > m, fprintf('  ...'); end
    fprintf('\n');
end
