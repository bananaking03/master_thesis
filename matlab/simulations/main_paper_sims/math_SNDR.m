% clc; clear; close all;

%% Parameters
N_bits = 9;                 % ADC resolution
VFS = 2;                    % Full-scale voltage (e.g., 1V)
mu = 0.25;                  % Calibration constant (IMPORTANT: must be small now!)
N_vals = logspace(2, 6, 200);   % cal_len range

%% Ideal SNDR
SNDR_ideal = 6.02*N_bits + 1.76;

%% Compute SNDR
scaling = (6 * mu * 2^(2*N_bits)) / (VFS^2);
SNDR = SNDR_ideal - 10*log10(1 + scaling ./ N_vals);

%% Plot
figure;
semilogx(N_vals, SNDR, 'LineWidth', 2);
grid on;
xlabel('cal\_len (N)');
ylabel('SNDR (dB)');
title(['SNDR vs cal\_len (Physical Scaling, ', num2str(N_bits), '-bit ADC)']);

hold on;
yline(SNDR_ideal, '--r', 'Ideal SNDR');

legend('With calibration noise', 'Ideal SNDR');