clear; clc;

%% Parameters
k = 1.38064852e-23;   % Boltzmann constant
T = 300;              % Temperature (K)

Nbits = 10;           % ADC resolution
Vfs = 2;              % Full-scale peak-to-peak voltage (V)

%% Derived quantities
Delta = Vfs / (2^Nbits);         % LSB size
Pq = (Delta^2) / 12;             % Quantization noise power

Ps = (Vfs^2) / 8;                % Signal power (full-scale sine)

%% Capacitance sweep
C = logspace(-15, -9, 200);      % 1 fF to 1 nF

%% kT/C noise
Pktc = k*T ./ C;

%% Total noise
Ptotal = Pq + Pktc;

%% SNDR calculation
SNDR = Ps ./ Ptotal;
SNDR_dB = 10*log10(SNDR);

%% Ideal SNDR (no kT/C)
SNDR_ideal = Ps / Pq;
SNDR_ideal_dB = 10*log10(SNDR_ideal);

%% Plot
figure;
semilogx(C, SNDR_dB, 'LineWidth', 2); hold on;
yline(SNDR_ideal_dB, '--r', 'Ideal SNDR');

grid on;
xlabel('Sampling Capacitance (F)');
ylabel('SNDR (dB)');
title('SNDR vs Sampling Capacitance (including kT/C noise)');
legend('With kT/C noise', 'Ideal (quantization-limited)', 'Location', 'best');