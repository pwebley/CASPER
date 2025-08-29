function status = customOdePlot(t, y, flag)
    persistent h indices ax

    if isempty(indices)
        indices = [1, 10, 20];  % Only plot these indices
    end

    switch flag
        case 'init'
            ax = gca;
            hold(ax, 'on');
            h = gobjects(size(indices));
            for i = 1:length(indices)
                h(i) = plot(ax, t, y(indices(i)), '.-', 'DisplayName', sprintf('y(%d)', indices(i)));
            end
            legend(ax, 'show');

        case ''
            for i = 1:length(indices)
                set(h(i), 'XData', [get(h(i), 'XData'), t], ...
                          'YData', [get(h(i), 'YData'), y(indices(i))]);
            end

            % === Diagnostic: identify the max y-element at this time ===
            [maxVal, maxIdx] = max(y);
            fprintf('t = %.3f s | Max y: %.4f at index %d\n', t, maxVal, maxIdx);

            % Optional: visually mark it
            plot(ax, t, maxVal, 'ko', 'MarkerSize', 6, ...
                 'DisplayName', sprintf('Max: y(%d)', maxIdx));

        case 'done'
            % Nothing to clean up
    end

    status = 0;
end
