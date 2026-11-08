clear;
clc;
close all;

%% Parameters

N_bits = 13;
N_caps = 2^N_bits;

sigma_mismatch = 0.01;  % 1% standard deviation
C_unit = 1e-15;         % Nominal unit capacitance: 1 fF

rng(42);                % Fixed seed for reproducibility

%% Generate mismatch for two independent DACs

N_caps = 8192;
sigma_mismatch = 0.01;
C_unit = 1e-15;

rng(42);

% DAC 1
C1 = C_unit * (1 + sigma_mismatch * randn(N_caps, 1));

% DAC 2
C2 = C_unit * (1 + sigma_mismatch * randn(N_caps, 1));

% Normalized output levels
dac_levels_1 = [0; cumsum(C1)] / sum(C1);
dac_levels_2 = [0; cumsum(C2)] / sum(C2);

% Save DAC 1
fid = fopen('dac1_mismatch_lut.txt', 'w');
fprintf(fid, '%.16e\n', dac_levels_1);
fclose(fid);

% Save DAC 2
fid = fopen('dac2_mismatch_lut.txt', 'w');
fprintf(fid, '%.16e\n', dac_levels_2);
fclose(fid);

%% Plot results

figure('Color', 'w');

plot(0:N_caps, INL, 'LineWidth', 1);

xlabel('DAC code');
ylabel('INL [LSB]');
title('DAC Integral Nonlinearity');
grid on;

set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');

figure('Color', 'w');

plot(1:N_caps, DNL, 'LineWidth', 1);

xlabel('DAC code');
ylabel('DNL [LSB]');
title('DAC Differential Nonlinearity');
grid on;

set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');