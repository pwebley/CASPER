% File: odeProgress.m
function stop = odeProgress(t, y, flag)
    stop = false;
    persistent lastT
    switch flag
        case 'init'
            fprintf('Solver started...\n');
            lastT = -inf;
        case ''
            tnow = t(end);
            if tnow - lastT >= 1  % Update every 1s
                fprintf('  Solver time: %.2f s\n', tnow);
                lastT = tnow;
            end
        case 'done'
            fprintf('Solver finished step.\n');
    end
end
