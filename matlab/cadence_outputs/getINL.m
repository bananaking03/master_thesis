clear;
clc;
close all;

%% Load exported Cadence data
filename = '9f_DAC_sweep.matlab';

% Read the comma-separated data
data = readtable(filename, ...
    'FileType', 'text', ...
    'Delimiter', ',', ...
    'VariableNamingRule', 'preserve');

% Extract time and output voltage
time = data{:,1};
Vout = data{:,2};

% Convert to column vectors
time = time(:);
Vout = Vout(:);

%% DAC parameters
N = length(Vout);       % Number of samples
N_steps = N - 1;        % Number of DAC transitions

%% Calculate endpoint INL
Vmin = Vout(1);
Vmax = Vout(end);

% Endpoint-fit LSB
LSB = (Vmax - Vmin) / N_steps;

% Ideal output for each code
code = (0:N_steps).';
Videal = Vmin + code * LSB;

% INL in volts and LSB
INL_V = Vout - Videal;
INL_LSB = INL_V / LSB;

%% Display results
fprintf('Number of samples: %d\n', N);
fprintf('Number of DAC steps: %d\n', N_steps);
fprintf('Endpoint LSB: %.10f V\n', LSB);
fprintf('Maximum INL: %.6f LSB\n', max(INL_LSB));
fprintf('Minimum INL: %.6f LSB\n', min(INL_LSB));
fprintf('Peak-to-peak INL: %.6f LSB\n', ...
    max(INL_LSB) - min(INL_LSB));

%% Plot INL
figure;
plot(code, INL_LSB, 'LineWidth', 1.5);
grid on;
xlabel('DAC code');
ylabel('INL (LSB)');
title('DAC Integral Nonlinearity (Endpoint Method)');

%% Plot transfer characteristic
figure;
plot(code, Vout, 'LineWidth', 1.5);
hold on;
plot(code, Videal, '--', 'LineWidth', 1.5);
grid on;
xlabel('DAC code');
ylabel('Output voltage (V)');
title('DAC Transfer Characteristic');
legend('Measured output', 'Ideal endpoint line', ...
    'Location', 'best');