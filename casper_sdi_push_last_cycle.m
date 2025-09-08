function casper_sdi_push_last_cycle(matfile)
% Push last-cycle boundary histories to Simulink Data Inspector (SDI).
% SDI is great for interactive time exploration: zoom, cursors, compare.

S = load(matfile);
sim = S.sim; steps = S.last_cycle.steps;

signals = struct('name',{},'ts',{});

for b = 1:sim.num_beds
  bedlabel = sprintf('Bed%c','A'+b-1);
  for endlabel = ["z0","zL"]
    % concat histories
    t=[]; P=[]; T=[]; u=[]; n=[];
    ymat = [];  % will build ns columns

    for s=1:numel(steps)
        H = steps(s).hist.(bedlabel).(endlabel);
        t = [t; H.t(:)];
        P = [P; H.P(:)];
        T = [T; H.T(:)];
        u = [u; H.u(:)];
        n = [n; H.n(:)];
        ymat = [ymat; H.y]; %#ok<AGROW>
    end

    % timeseries
    signals(end+1).name = sprintf('%s_%s_P',bedlabel,endlabel); %#ok<AGROW>
    signals(end).ts = timeseries(P, t);
    signals(end+1).name = sprintf('%s_%s_T',bedlabel,endlabel);
    signals(end).ts = timeseries(T, t);
    signals(end+1).name = sprintf('%s_%s_u',bedlabel,endlabel);
    signals(end).ts = timeseries(u, t);
    signals(end+1).name = sprintf('%s_%s_n',bedlabel,endlabel);
    signals(end).ts = timeseries(n, t);

    for j = 1:sim.n_species
      signals(end+1).name = sprintf('%s_%s_y%d',bedlabel,endlabel,j);
      signals(end).ts = timeseries(ymat(:,j), t);
    end
  end
end

run_label = S.run_info.run_name;
runID = Simulink.sdi.createRun(run_label, 'vars', signals); %#ok<NASGU>
Simulink.sdi.view;
end
