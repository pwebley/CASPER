function bed_states = unpack_bed_state_vector(Y, sim)
% UNPACK_BED_STATE_VECTOR Unpacks flattened ODE vector into bed states.
%
% Inputs:
%   Y   - ODE state vector at a given time [column vector]
%   sim - simulation structure with configuration
%
% Output:
%   bed_states - cell array of bed state structs {1:sim.num_beds}

    NS = sim.num_nodes;
    N  = sim.n_species;
    bed_states = cell(sim.num_beds, 1);
    
    % Total length of the state vector per bed
    varlen = NS + NS * (N-1) + NS + NS * N;
    
    for b = 1:sim.num_beds
        offset = (b-1)*varlen;

        Ct = Y(offset + 1 : offset + NS);
        Ci = reshape(Y(offset + NS + 1 : offset + NS + NS*(N-1)), NS, N-1);
        T  = Y(offset + NS + NS*(N-1) + 1 : offset + NS + NS*(N-1) + NS);
        qi = reshape(Y(offset + NS + NS*(N-1) + NS + 1 : offset + varlen), NS, N);

        % Compute full species concentrations
        Ci_full = zeros(NS, N);
        Ci_full(:, 1:N-1) = Ci;
        Ci_full(:, N) = Ct - sum(Ci, 2);
        yi = Ci_full ./ Ct;

        % Store in bed struct
        bed.P = Ct .* sim.R .* T;
        bed.C = Ct;
        bed.C_species = Ci;
        bed.T = T;
        bed.q = qi;
        bed.y = yi;

        % Placeholder velocity fields
        bed.u_z0 = 0;
        bed.u_zL = 0;
        bed.current_step_config = [];

        bed_states{b} = bed;
    end
end
