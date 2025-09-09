function q = get_equilibrium_loading(sim, adsorbent, P_partial, T)
% GET_EQUILIBRIUM_LOADING
% Supports your existing models and adds 'coads_wadst' for DAC:
%  - H2O by GAB with temperature-dependent c,k
%  - CO2 by WADST = weighted average of dry/wet Toth
%  - N2 linear (very small K)

    n = sim.n_species;
    q = zeros(1, n);
    R = 8.314462618; % J/mol-K
    
    model  = adsorbent.isotherm_type;
    params = adsorbent.isotherms;

    % ---------- SPECIAL CASE: WADST (CO2/H2O co-adsorption) ----------
    if strcmpi(model, 'coads_wadst')
        % Map species indices (robust to order)
        idx = containers.Map(sim.species_names, num2cell(1:n));
        iCO2 = idx('CO2');  iH2O = idx('H2O');  iN2 = idx('N2');
        pCO2 = P_partial(iCO2);  pH2O = P_partial(iH2O);

        % ---- 1) H2O by GAB (Anderson form) ----
        gab = params.H2O.gab;    % fields: qm, C, D, F, G
        % saturation vapor pressure (Tetens-like; good up to ~100C)
        psat = 610.94 * exp(17.625*(T-273.15)/(T-30.11));  % Pa
        RH   = max(0, min(0.9999, pH2O / psat));           % in [0,1)
        % Energies:
        E10 = -44.38*T + 57220;                    % J/mol
        E1  = gab.C - exp(gab.D*T);                % J/mol
        E29 = gab.F + gab.G*T;                     % J/mol
        c   = exp((E1  - E10)/(R*T));
        k   = exp((E29 - E10)/(R*T));
        qm  = gab.qm;                               % mol/kg
        qH2O = qm * (k*c*RH) / ((1 - k*RH) * (1 + (c-1)*k*RH));
        q(iH2O) = qH2O;

        % ---- 2) CO2 by Toth (dry & wet) ----
        CO2d = params.CO2.dry_toth;   % fields: T0,q1_0,w,b0,DH0,t0,a
        CO2w = params.CO2.wet_toth;   % same fields

        toth_eval = @(tp, p, T) ...
            (tp.q1_0 * exp(tp.w*(1 - T/tp.T0))) .* ...
            ( (tp.b0 * exp(-tp.DH0/(R*T))) .* p ) ./ ...
            ( 1 + ( (tp.b0 * exp(-tp.DH0/(R*T))) * p ).^(tp.t0 + tp.a*(1 - tp.T0/T)) ).^(1./(tp.t0 + tp.a*(1 - tp.T0/T)));

        qCO2_dry = toth_eval(CO2d, pCO2, T);
        qCO2_wet = toth_eval(CO2w, pCO2, T);

        % ---- 3) WADST blend with critical water loading A ----
        A   = params.CO2.coads.A;     % mol/kg
        phi = exp(-A * max(qH2O, 0)); % probability kernel
        q(iCO2) = (1 - phi) * qCO2_dry + phi * qCO2_wet;

        % ---- 4) N2 (linear, near-zero) ----
        if isfield(params,'N2') && isfield(params.N2,'linear')
            q(iN2) = params.N2.linear.K * P_partial(iN2);
        end

        return;  % *** WADST handled completely; skip the generic path ***
    end

    % ---------- GENERIC PRECOMPUTE (your original pattern) ----------
    switch lower(model)
        case 'langmuir'
            sum_bP = 0;
            for i = 1:n
                gas_name   = sim.species_names{i};
                gas_params = params.(gas_name);
                Pi = P_partial(i);
                % NOTE: many forms use b = b0*exp(-dH/RT). Keep your sign if desired.
                bi = gas_params.b0 * exp(-gas_params.dH / (R * T));
                sum_bP = sum_bP + bi * Pi;
            end

        case 'dsl'
            sum_b1P = 0;  sum_b2P = 0;
            for i = 1:n
                gas_name   = sim.species_names{i};
                gas_params = params.(gas_name);
                Pi = P_partial(i);
                b1_i = gas_params.b01 * exp(-gas_params.dH1 / (R * T));
                b2_i = gas_params.b02 * exp(-gas_params.dH2 / (R * T));
                sum_b1P = sum_b1P + b1_i * Pi;
                sum_b2P = sum_b2P + b2_i * Pi;
            end

        case {'langmuir_freundlich','sips','lang_freund'}
            sum_bP_n = 0;
            for i = 1:n
                gas_name   = sim.species_names{i};
                gas_params = params.(gas_name);
                Pi = P_partial(i);
                bi = gas_params.b0 * exp(-gas_params.dH / (R * T));
                sum_bP_n = sum_bP_n + (bi * Pi)^gas_params.n;
            end

        otherwise
            % linear and others don't need precompute
    end

    % ---------- PER-SPECIES LOADINGS (your original loop) ----------
    for i = 1:n
        gas_name   = sim.species_names{i};
        gas_params = params.(gas_name);
        Pi = P_partial(i);

        switch lower(model)
            case 'langmuir'
                bi   = gas_params.b0 * exp(-gas_params.dH / (R * T));
                q(i) = gas_params.qmax * (bi * Pi) / (1 + sum_bP);

            case 'linear'
                q(i) = gas_params.K * Pi;

            case 'dsl'
                b1_i = gas_params.b01 * exp(-gas_params.dH1 / (R * T));
                b2_i = gas_params.b02 * exp(-gas_params.dH2 / (R * T));
                q(i) = (gas_params.m1 * b1_i * Pi) / (1 + sum_b1P) + ...
                       (gas_params.m2 * b2_i * Pi) / (1 + sum_b2P);

            case {'langmuir_freundlich','sips','lang_freund'}
                bi    = gas_params.b0 * exp(-gas_params.dH / (R * T));
                bP_n  = (bi * Pi)^gas_params.n;
                q(i)  = gas_params.qmax * bP_n / (1 + sum_bP_n);

            otherwise
                error('Unsupported isotherm model: %s', model);
        end
    end
end
