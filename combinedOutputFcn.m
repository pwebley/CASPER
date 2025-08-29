function status = combinedOutputFcn(t, y, flag)
    % Call your custom plot function
     status1 = customOdePlot(t, y, flag);

    % Call the odeProgress function
    status2 = odeProgress(t, y, flag);

    % Return true if either signals a stop
     status = status1 || status2;
end
