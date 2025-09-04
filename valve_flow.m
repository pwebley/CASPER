function Q = valve_flow(P1, P2, T, y, sim, Cv)
% VALVE_FLOW  Regularized gas valve equation -> STP volumetric flow [L/s]
% Keeps your original physics but blends smoothly between choked & subsonic.

    % --- Constants & guards
    Rcrit   = 0.53;      % critical pressure ratio
    band    = 0.04;      % smoothing half-width in ratio units (±4%)  <-- tune 0.02–0.07
    Tmin    = 200;       % K, floor to avoid 1/sqrt(T) blowups
    Pminatm = 1e-4;      % atm, to avoid divide-by-zero in ratio

    % Convert to atm for your existing coefficient (matches your original)
    P1_atm = max(P1/101325, Pminatm);
    P2_atm = max(P2/101325, 0);

    % No forward flow if P1<=P2 (soft-clip with tiny epsilon if you prefer)
    if P1_atm <= P2_atm
        Q = 0;
        return;
    end

    % Mixture MW and specific gravity
    MW_species = sim.MW(:) / 1000;       % kg/mol
    MW_mix     = max(y(:)'*MW_species, 1e-6);
    G          = MW_mix / 0.02897;       % relative to air

    T_eff = max(T, Tmin);
    r     = min(P2_atm / P1_atm, 0.999999);  % ratio, capped below 1

    % Your two modes (same formulas as you had)
    Q_sub    = 77.01 * Cv * sqrt( max(P1_atm^2 - P2_atm^2, 0) / (G * T_eff) );
    Q_choked = 77.01 * Cv * ( P1_atm / sqrt(G * T_eff) );

    % Smooth blend across r ~ Rcrit using a tanh window
    % s ≈ 0 => choked, s ≈ 1 => subsonic
    s   = 0.5 * (1 + tanh( (r - Rcrit) / band ));
    Q   = s * Q_sub + (1 - s) * Q_choked;

    % Numerical hygiene (exact zeros when very close)
    if Q < 1e-12
        Q = 0;
    end
end
