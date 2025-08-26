function run_casper_simulation
% RUN_CASPER Run PSA simulation with multiple cycles and progress display.
%
% Usage: run_casper(5)  % Runs 5 cycles
if nargin < 1
    num_cycles = 1;
end

[sim, bed_states, tank_states] = createPSASimulation();
num_cycles=sim.n_cycles; 

% Flatten initial state into Y0
Y0 = [];
for b = 1:sim.num_beds
    bed = bed_states{b};
    Y0 = [Y0; bed.C; reshape(bed.C_species, [],1); bed.T; reshape(bed.q, [],1)];
end

opts = odeset('RelTol', sim.RelTol, 'AbsTol', sim.AbsTol, ...
    'MaxStep', sim.dt_max);

% Optional: use OutputFcn for progress updates
opts = odeset(opts, 'OutputFcn', @(t,y,flag) odeProgress(t,y,flag));

T_all = [];
Y_all = [];

for cycle = 1:num_cycles
    fprintf('\n=== Beginning cycle %d of %d ===\n', cycle, num_cycles);
    
    for idx = 1:length(sim.step_times)-1
        t0 = sim.step_times(idx);
        tf = sim.step_times(idx+1);
        current_step = sim.step(idx);
        
        fprintf('  Cycle %d — Step %d: t = %.2f to %.2f s\n', cycle, idx, t0, tf);
        
        % Integration across this step
        [Tseg, Yseg] = ode15s(@(t,Y) multi_bed_ode(t, Y, sim, current_step), ...
                              [t0, tf], Y0, opts);
        
        T_all = [T_all; Tseg];
        Y_all = [Y_all; Yseg];
        Y0 = Yseg(end,:)';  % Update initial state for next step
    end
end

% Save final state for restart
bed_states_final = reshape(Y0, [], sim.num_beds);
save('cycle_sim.mat', 'sim', 'bed_states_final', 'tank_states');

fprintf('\nSimulation complete! Final state saved to cycle_sim.mat\n');

% Example output: final concentration profile of Bed A
Ct_end = Y0(1:sim.num_nodes);
figure; plot(sim.z_nodes, Ct_end, 'LineWidth', 2);
xlabel('z (m)'); ylabel('Total Concentration (mol/m³)');
title('Final Bed A Concentration Profile'); grid on;
end
