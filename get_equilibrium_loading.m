function q = get_equilibrium_loading(sim, adsorbent, P_partial, T)
    n = sim.n_species;
    q = zeros(1, n);
    R = 8.314462618; % J/mol·K
    
    model = adsorbent.isotherm_type;
    params = adsorbent.isotherms;
    
    % Pre-calculate competitive terms if needed
    switch lower(model)
        case 'langmuir'
            sum_bP = 0;
            for i = 1:n
                gas_name = sim.species_names{i};
                gas_params = params.(gas_name);
                Pi = P_partial(i);
                bi = gas_params.b0 * exp(-gas_params.dH / (R * T));
                sum_bP = sum_bP + bi * Pi;
            end
            
        case 'dsl'
            sum_b1P = 0;
            sum_b2P = 0;
            for i = 1:n
                gas_name = sim.species_names{i};
                gas_params = params.(gas_name);
                Pi = P_partial(i);
                b1_i = gas_params.b01 * exp(-gas_params.dH1 / (R * T));
                b2_i = gas_params.b02 * exp(-gas_params.dH2 / (R * T));
                sum_b1P = sum_b1P + b1_i * Pi;
                sum_b2P = sum_b2P + b2_i * Pi;
            end
            
        case {'langmuir_freundlich', 'sips', 'lang_freund'}
            sum_bP_n = 0;
            for i = 1:n
                gas_name = sim.species_names{i};
                gas_params = params.(gas_name);
                Pi = P_partial(i);
                bi = gas_params.b0 * exp(-gas_params.dH / (R * T));
                % Calculate (b_i * P_i)^n_i
                bP_n = (bi * Pi)^gas_params.n;
                sum_bP_n = sum_bP_n + bP_n;
            end
    end
    
    % Calculate loadings
    for i = 1:n
        gas_name = sim.species_names{i};
        gas_params = params.(gas_name);
        Pi = P_partial(i);
        
        switch lower(model)
            case 'langmuir'
                bi = gas_params.b0 * exp(-gas_params.dH / (R * T));
                q(i) = gas_params.qmax * (bi * Pi) / (1 + sum_bP);
                
            case 'linear'
                q(i) = gas_params.K * Pi;
                
            case 'dsl'
                b1_i = gas_params.b01 * exp(-gas_params.dH1 / (R * T));
                b2_i = gas_params.b02 * exp(-gas_params.dH2 / (R * T));
                q(i) = (gas_params.m1 * b1_i * Pi) / (1 + sum_b1P) + ...
                       (gas_params.m2 * b2_i * Pi) / (1 + sum_b2P);
                
            case {'langmuir_freundlich', 'sips', 'lang_freund'}
                bi = gas_params.b0 * exp(-gas_params.dH / (R * T));
                % Calculate (b_i * P_i)^n_i
                bP_n = (bi * Pi)^gas_params.n;
                q(i) = gas_params.qmax * bP_n / (1 + sum_bP_n);
                
            otherwise
                error('Unsupported isotherm model: %s', model);
        end
    end
end