function [outputArg1,outputArg2] = untitled11(inputArg1,inputArg2)
% Parametres for N2
b0_paper=1.2305;  %1/atm
d0_paper=0.0740;  %1/atm
m1_paper=0.69;    %mmol/g
m2_paper=3.34;    %mmol/g
dH1_paper=5461.1; %cal/mol
dH1_paper=dH1_paper*4.184; %convert to J/mol
dH2_paper=4265.6; %cal/mol
dH2_paper=dH2_paper*4.184; %convert to J/mol
T0 = 296.15;      %Reference T, 296.15K
R = 8.3144;       %Ideal Gas Constant, J/mol.K
b0_model = (b0_paper / 101325) * exp(-dH1_paper / (R * T0)); %b0 in terms of 1/Pa
d0_model = (d0_paper / 101325) * exp(-dH2_paper / (R * T0)); %d0 in terms of 1/Pa


% --- Pressure range (0 to 20 atm converted to Pa) ---
P_atm = linspace(0, 20, 200);     % atm
P_Pa = P_atm * 101325;            % convert to Pa

% --- Evaluation Temperatures (in K) ---
T_list = [296, 318];
colors = ['b', 'r'];
legend_entries = cell(1, numel(T_list));
q_all = zeros(numel(T_list), numel(P_atm));

% --- Loop Over Temperatures and Compute Isotherms ---
figure;
hold on;
for idx = 1:numel(T_list)
    T = T_list(idx);
% --- Calculate b(T) and d(T) at this T ---
b1 = b0_model * exp(dH1_paper / (R * T));
b2 = d0_model * exp(dH2_paper / (R * T));

% --- DSL isotherm equation (q in mmol/g) ---
q = m1_paper * (b1 .* P_Pa) ./ (1 + b1 .* P_Pa) + ...
    m2_paper * (b2 .* P_Pa) ./ (1 + b2 .* P_Pa);
    % Store and Plot
    q_all(idx, :) = q;
    plot(P_atm, q, [colors(idx) '-'], 'LineWidth', 2);
    legend_entries{idx} = [num2str(T) ' K'];
end

% --- Final Plot Formatting ---
xlabel('Pressure [atm]');
ylabel('q [mmol/g]');
title('N_2 DSL Isotherm on Zeolite 5A');
legend(legend_entries, 'Location', 'southeast');
grid on;
hold off;


% Parameters for O2
b0_paper=0.0518;  %1/atm
d0_paper=0.0518;  %1/atm
m1_paper=0.69;    %mmol/g
m2_paper=3.34;    %mmol/g
dH1_paper=3716.7; %cal/mol
dH1_paper=dH1_paper*4.184; %convert to J/mol
dH2_paper=3716.7; %cal/mol
dH2_paper=dH2_paper*4.184; %convert to J/mol
T0 = 296.15;      %Reference T, 296.15K
R = 8.3144;       %Ideal Gas Constant, J/mol.K
b0_model = (b0_paper / 101325) * exp(-dH1_paper / (R * T0)); %b0 in terms of 1/Pa
d0_model = (d0_paper / 101325) * exp(-dH2_paper / (R * T0)); %d0 in terms of 1/Pa


% --- Pressure range (0 to 20 atm converted to Pa) ---
P_atm = linspace(0, 20, 200);     % atm
P_Pa = P_atm * 101325;            % convert to Pa

% --- Evaluation Temperatures (in K) ---
T_list = [296, 318];
colors = ['b', 'r'];
legend_entries = cell(1, numel(T_list));
q_all = zeros(numel(T_list), numel(P_atm));

% --- Loop Over Temperatures and Compute Isotherms ---
figure;
hold on;
for idx = 1:numel(T_list)
    T = T_list(idx);
% --- Calculate b(T) and d(T) at this T ---
b1 = b0_model * exp(dH1_paper / (R * T));
b2 = d0_model * exp(dH2_paper / (R * T));

% --- DSL isotherm equation (q in mmol/g) ---
q = m1_paper * (b1 .* P_Pa) ./ (1 + b1 .* P_Pa) + ...
    m2_paper * (b2 .* P_Pa) ./ (1 + b2 .* P_Pa);
    % Store and Plot
    q_all(idx, :) = q;
    plot(P_atm, q, [colors(idx) '-'], 'LineWidth', 2);
    legend_entries{idx} = [num2str(T) ' K'];
end

% --- Final Plot Formatting ---
xlabel('Pressure [atm]');
ylabel('q [mmol/g]');
title('O_2 DSL Isotherm on Zeolite 5A');
legend(legend_entries, 'Location', 'southeast');
grid on;
hold off;

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end