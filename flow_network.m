function [inlet_bc, outlet_bc] = flow_network(t, sim, raw, current_step)
% FLOW_NETWORK_UPDATE computes the inlet and outlet boundary conditions 
% for each bed based on the current step's configuration and flow laws.

num_beds = sim.num_beds;
R = sim.R;

inlet_bc = cell(num_beds, 1);
outlet_bc = cell(num_beds, 1);

for b = 1:num_beds
    % Get interface conditions for this bed
    bed_raw = raw{b};

    % Determine which config applies (BedA or BedB)
    if b == 1
        config = current_step.BedA;
    else
        config = current_step.BedB;
    end

    %% --- z = 0 Interface ---
    bc0 = struct();
    bc0.P = bed_raw.inlet.P;
    bc0.T = bed_raw.inlet.T;
    bc0.y = bed_raw.inlet.y;
    bc0.F = 0;  % Default flow rate

    if isfield(config, 'z0') && ~isempty(config.z0.flow_law)
        flow_law = lower(config.z0.flow_law{1});
        params = config.z0.parameters{1};

        switch flow_law
            case 'valve'
                % Assume external pressure is known from sim or tank
                if isfield(config.z0, 'source')
                    src = config.z0.source{1};
                    P_ext = get_external_pressure(sim, src);
                else
                    P_ext = bed_raw.inlet.P;
                end
                dP = P_ext - bed_raw.inlet.P;
                bc0.F = params.Cv * dP;

            case 'molar_flow'
                bc0.F = params.Fm;

            case 'volume_flow'
                % F_vol in m³/s -> convert to mol/s using ideal gas law
                F_vol = params.Fv;
                bc0.F = (F_vol * bed_raw.inlet.P) / (R * bed_raw.inlet.T);

            case 'none'
                bc0.F = 0;

            otherwise
                error('Unsupported flow law at z=0: %s', flow_law);
        end
    end

    %% --- z = L Interface ---
    bcL = struct();
    bcL.P = bed_raw.outlet.P;
    bcL.T = bed_raw.outlet.T;
    bcL.y = bed_raw.outlet.y;
    bcL.F = 0;

    if isfield(config, 'zL') && ~isempty(config.zL.flow_law)
        flow_law = lower(config.zL.flow_law{1});
        params = config.zL.parameters{1};

        switch flow_law
            case 'valve'
                if isfield(config.zL, 'destination')
                    dst = config.zL.destination{1};
                    P_ext = get_external_pressure(sim, dst);
                else
                    P_ext = bed_raw.outlet.P;
                end
                dP = bed_raw.outlet.P - P_ext;
                bcL.F = params.Cv * dP;

            case 'molar_flow'
                bcL.F = params.Fm;

            case 'volume_flow'
                F_vol = params.Fv;
                bcL.F = (F_vol * bed_raw.outlet.P) / (R * bed_raw.outlet.T);

            case 'none'
                bcL.F = 0;

            otherwise
                error('Unsupported flow law at z=L: %s', flow_law);
        end
    end

    inlet_bc{b} = bc0;
    outlet_bc{b} = bcL;
end
end
