function run_casper_simulation
% RUN_CASPER Run PSA simulation with multiple cycles and progress display.
%

[sim, bed_states, tank_states] = createPSASimulation();
for b = 1:sim.num_beds
    sim.diagnostics(b).time = [];
    sim.diagnostics(b).Ct_1 = [];
    sim.diagnostics(b).Ct_N = [];
    sim.diagnostics(b).y1_1 = [];
    sim.diagnostics(b).y1_N = [];
    sim.diagnostics(b).P_1 = [];
    sim.diagnostics(b).q1_1 = [];
    sim.diagnostics(b).q1_N = [];
    sim.diagnostics(b).q2_1 = [];
    sim.diagnostics(b).q2_N = [];
    sim.diagnostics(b).T_1 = [];
    sim.diagnostics(b).T_N = [];
end

num_cycles=sim.n_cycles; 

% Build initial state  Y0
Y0 = [];
for b = 1:sim.num_beds
    bed = bed_states{b};
    Y0 = [Y0; bed.C; reshape(bed.C_species, [],1); bed.T; reshape(bed.q, [],1)];
end

% opts = odeset('RelTol', sim.RelTol, ...
%               'AbsTol', sim.AbsTol, ...
%               'Stats', 'on', ...
%               'MaxStep', sim.dt_max);
opts = odeset('RelTol', sim.RelTol, ...
              'AbsTol', sim.AbsTol, ...
              'Stats', 'on', ...
              'MaxStep', sim.dt_max, ...
              'OutputFcn', @odeProgress);
% 
T_all = [];
Y_all = [];
bed_states_all = cell(num_cycles, 1);  % Store final bed states per cycle
tank_states_all = cell(num_cycles, 1); % (optional) for tanks too

for cycle = 1:num_cycles
    fprintf('\n=== Beginning cycle %d of %d ===\n', cycle, num_cycles);
    
    for idx = 1:length(sim.step_times)-1
        t0 = sim.step_times(idx);
        tf = sim.step_times(idx+1);
        current_step = sim.step(idx);
        
        fprintf('  Cycle %d — Step %d: t = %.2f to %.2f s\n', cycle, idx, t0, tf);
        
        % Integration across this step
        [Tseg, Yseg] = ode23t(@(t,Y) multi_bed_ode(t, Y, sim, current_step), ...
                              [t0, tf], Y0, opts);
        
        T_all = [T_all; Tseg];
        Y_all = [Y_all; Yseg];
        Y0 = Yseg(end,:)';  % Update initial state for next step
        plot_diagnostics(sim);
        pause;

    end
    bed_states_all{cycle} = unpack_bed_state_vector(Y0, sim);
    tank_states_all{cycle} = tank_states;  % If tanks are updated per step

end

save('cycle_sim.mat', 'sim', 'bed_states_all', 'tank_states_all');

fprintf('\n✅ Simulation complete. Final state saved to cycle_sim.mat\n');


% Example output: final concentration profile of Bed A
Ct_end = Y0(1:sim.num_nodes);
figure; plot(sim.z_nodes, Ct_end, 'LineWidth', 2);
xlabel('z (m)'); ylabel('Total Concentration (mol/m³)');
title('Final Bed A Concentration Profile'); grid on;
end


function status = combinedOutputFcn(t, y, flag)
    status1 = odeProgress(t, y, flag);
    status2 = odeplot(t, y, flag);
    status = status1;
    status = status1 || status2;
end
