function [k, cfg] = casper_step_for_time(sim, steps, t_now)
% Returns step index k and its config for global time t_now within a cycle.

% Assumes sim.step_times defines edges in [0, sim.cycle_time]
edges = sim.step_times(:);
if t_now <= edges(1), k = 1;
elseif t_now >= edges(end), k = numel(steps);
else
    k = find(edges <= t_now, 1, 'last');
    if k == numel(edges), k = k-1; end
end
cfg = steps(k).config;
end
