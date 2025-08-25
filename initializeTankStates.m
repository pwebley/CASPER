function tank_states = initializeTankStates(sim)
%INITIALIZETANKSTATES Returns initialized tank state structs.

R = 8.314;
tank_states = struct('name', {}, 'P', {}, 'T', {}, 'n', {}, 'y', {}, 'molar_hold_up', {});

for k = 1:numel(sim.tanks)
    tank = sim.tanks(k);
    ts.name = tank.name;
    ts.P = tank.P;
    ts.T = tank.T;
    ts.y = tank.y(:)';

    if strcmpi(tank.type, 'finite')
        ts.n = (ts.P * tank.volume) / (R * ts.T);
    else
        ts.n = 1e6;
    end

    ts.molar_hold_up = ts.n * ts.y;
    ts.initial_hold_up = ts.molar_hold_up;
    tank_states(k) = ts;
end
end
