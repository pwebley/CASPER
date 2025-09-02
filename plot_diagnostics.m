function plot_diagnostics(sim)
% PLOT_DIAGNOSTICS Visualize key PSA variables after a cycle step.

num_beds = sim.num_beds;
N = sim.n_species;
bed_labels = arrayfun(@(b) sprintf('Bed%c', 'A' + b - 1), 1:num_beds, 'UniformOutput', false);

% === Plot A: Total concentration at node 1 and N ===
figure('Name', 'Total Concentration');
hold on;
for b = 1:num_beds
    t = sim.diagnostics(b).time;
    plot(t, sim.diagnostics(b).Ct_1, '-', 'DisplayName', [bed_labels{b}, ' Ct@1']);
    plot(t, sim.diagnostics(b).Ct_N, '--', 'DisplayName', [bed_labels{b}, ' Ct@N']);
end
xlabel('Time (s)'); ylabel('Total Concentration (mol/m³)');
legend; title('Total Concentration at Nodes 1 and N');

% === Plot B: Mole fraction of species 1 at nodes 1 and N ===
figure('Name', 'Mole Fraction of Species 1');
hold on;
for b = 1:num_beds
    t = sim.diagnostics(b).time;
    plot(t, sim.diagnostics(b).y1_1, '-', 'DisplayName', [bed_labels{b}, ' y1@1']);
    plot(t, sim.diagnostics(b).y1_N, '--', 'DisplayName', [bed_labels{b}, ' y1@N']);
end
xlabel('Time (s)'); ylabel('Mole Fraction of Species 1');
legend; title('Mole Fraction at Nodes 1 and N');

% === Plot C: Pressure at node 1 ===
figure('Name', 'Pressure at Node 1');
hold on;
for b = 1:num_beds
    t = sim.diagnostics(b).time;
    plot(t, sim.diagnostics(b).P_1, '-', 'DisplayName', [bed_labels{b}, ' P@1']);
end
xlabel('Time (s)'); ylabel('Pressure (Pa)');
legend; title('Pressure at Node 1');

% === Plot D: Loadings of Species 1 and 2 at nodes 1 and N ===
figure('Name', 'Adsorbed Loadings');
hold on;
for b = 1:num_beds
    t = sim.diagnostics(b).time;
    plot(t, sim.diagnostics(b).q1_1, '-', 'DisplayName', [bed_labels{b}, ' q1@1']);
    plot(t, sim.diagnostics(b).q1_N, '--', 'DisplayName', [bed_labels{b}, ' q1@N']);
    plot(t, sim.diagnostics(b).q2_1, '-.', 'DisplayName', [bed_labels{b}, ' q2@1']);
    plot(t, sim.diagnostics(b).q2_N, ':', 'DisplayName', [bed_labels{b}, ' q2@N']);
end
xlabel('Time (s)'); ylabel('Loading (mol/kg)');
legend; title('Adsorbed Loadings at Nodes 1 and N');

% === Plot E: Temperature at nodes 1 and N ===
figure('Name', 'Temperature');
hold on;
for b = 1:num_beds
    t = sim.diagnostics(b).time;
    plot(t, sim.diagnostics(b).T_1, '-', 'DisplayName', [bed_labels{b}, ' T@1']);
    plot(t, sim.diagnostics(b).T_N, '--', 'DisplayName', [bed_labels{b}, ' T@N']);
end
xlabel('Time (s)'); ylabel('Temperature (K)');
legend; title('Temperature at Nodes 1 and N');
end
