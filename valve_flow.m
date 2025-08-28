function Q = valve_flow(P1, P2, T, y, sim, Cv)
% VALVE_FLOW Computes flow rate [L/s] through a valve using standard valve equation
% Inputs:
%   P1 – upstream pressure [Pa]
%   P2 – downstream pressure [Pa]
%   T  – temperature [K]
%   y  – mole fractions (column vector)
%   sim – simulation struct (must include sim.MW)
%   Cv – valve flow coefficient
% Output:
%   Q  – volumetric flow rate at STP [L/s]

    % Convert pressures to atm
    P1 = P1 / 101325;
    P2 = P2 / 101325;

    % Check valve behavior: flow only if P1 > P2
    if P1 <= P2
        Q = 0;
        return;
    end

    % Compute gas mixture molecular weight and specific gravity
    MW_species = sim.MW(:) / 1000;     % kg/mol
    MW_mix = y' * MW_species;          % mixture MW
    G = MW_mix / 0.02897;              % specific gravity rel. to air

    if P2 >= 0.53 * P1
        Q = 77.01 * Cv * sqrt((P1^2 - P2^2) / (G * T));
    else
        Q = 77.01 * Cv * (P1 / sqrt(G * T));
    end
end
