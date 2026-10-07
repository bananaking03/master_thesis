clear;
clc;
close all;

filenames = ["comp_-075mV.matlab", "comp_-05mV.matlab", "comp_-025mV.matlab", ...
             "comp_-0125mV.matlab", "comp_0mV.matlab", "comp_0125mV.matlab", ...
             "comp_025mV.matlab", "comp_05mV.matlab", "comp_075mV.matlab"];

Voltages = [-0.75e-3 -0.5e-3 -0.25e-3 -0.125e-3 0 0.125e-3 0.25e-3 ...
             0.5e-3 0.75e-3];

one_ratios = zeros(1,length(filenames));
N_ones = zeros(1,length(filenames));
N_total = zeros(1,length(filenames));

for i = 1:length(filenames)

    %% Load exported Cadence data
    filename = filenames(i);

    data = readtable(filename, ...
        'FileType', 'text', ...
        'Delimiter', ',', ...
        'VariableNamingRule', 'preserve');

    time = data{:,1};
    Vout = data{:,2};

    time = time(:);
    Vout = Vout(:);

    %% Convert comparator output to 0/1
    bin_out = Vout > 0;

    %% Statistics
    N_ones(i) = sum(bin_out);
    N_total(i) = length(bin_out);

    one_ratios(i) = N_ones(i) / N_total(i);
end


%% ============================================================
%  Fit entire probability curve
%  P(1) = normcdf((Vdiff - Vos)/sigma)
% =============================================================

% Initial guesses
Vos0 = 0;
sigma0 = 2e-3;

x0 = [Vos0 sigma0];

% Objective function
% x(1) = offset
% x(2) = input-referred noise
objective = @(x) -sum( ...
    N_ones .* log(normcdf((Voltages-x(1))/x(2))) + ...
    (N_total-N_ones) .* log(1-normcdf((Voltages-x(1))/x(2))) );

% Fit
options = optimset('Display','iter');

fit_params = fminsearch(objective, x0, options);

Vos = fit_params(1);
sigma_in = abs(fit_params(2));


%% Fitted curve
Vfit = linspace(min(Voltages), max(Voltages), 1000);

Pfit = normcdf((Vfit-Vos)/sigma_in);


%% Display results
fprintf('\nInput-referred offset = %.3f mV\n', Vos*1e3);
fprintf('Input-referred noise  = %.3f mV RMS\n', sigma_in*1e3);


%% Plot
figure;

plot(Voltages*1e3, one_ratios, 'o', 'LineWidth', 1.5);
hold on;

plot(Vfit*1e3, Pfit, 'LineWidth', 1.5);

xlabel('Differential input voltage [mV]');
ylabel('Probability of output = 1');
legend('Simulation data', 'Gaussian CDF fit', 'Location', 'best');
grid on;