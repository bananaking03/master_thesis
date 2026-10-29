clear;
clc;
close all;

% Read exported SimVision data H_matched
data = readtable('results_callen3000.csv', ...
    'VariableNamingRule', 'preserve');

%% Parameters
N_bits = 10;

f0 = 5435051;       % Input frequency (Hz)
Ts = 11e-9;         % Sampling period (s)
Fs = 1/Ts;          % Sampling frequency (Hz)

% Coherent sampling
N = 184;            % Number of FFT samples
k = 11;             % Number of input periods

f0_coherent = k*Fs/N;

fprintf('Sampling frequency: %.6f MHz\n', Fs/1e6);
fprintf('Input frequency: %.6f MHz\n', f0/1e6);
fprintf('Coherent frequency: %.6f MHz\n', f0_coherent/1e6);
fprintf('Frequency error: %.3f Hz\n', f0-f0_coherent);

%% Extract ADC output
% Column 1: Simulation time (ns)
% Column 2: ADC output (hexadecimal)

t = data.SimTime;
adc_hex = string(data{:,2});

% Convert hexadecimal to decimal
adc = double(hex2dec(adc_hex));

% Convert to signed 10-bit two's complement
adc(adc >= 2^(N_bits-1)) = ...
    adc(adc >= 2^(N_bits-1)) - 2^N_bits;

%% Plot ADC output versus time
figure('Color', 'w');

plot(t, adc, 'b-', 'LineWidth', 1);
grid on;
box on;

xlabel('Time (ns)', 'Color', 'k');
ylabel('ADC Output Code', 'Color', 'k');
title('ADC Output versus Time', 'Color', 'k');

ylim([-512 511]);

set(gca, 'Color', 'w', ...
         'XColor', 'k', ...
         'YColor', 'k', ...
         'FontSize', 12);

%% Select samples for coherent FFT
% Skip the first sample (initialization)
start_idx = 2;

if length(adc) < start_idx + N - 1
    error('Not enough ADC samples for N = %d.', N);
end

adc_fft = adc(start_idx:start_idx+N-1);

% Remove DC component
adc_fft = adc_fft - mean(adc_fft);

%% Calculate FFT
% Rectangular window (no windowing for coherent sampling)

Y = fft(adc_fft);

% Single-sided amplitude spectrum
P2 = abs(Y)/N;
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Convert to dBFS
% Full-scale sine-wave peak amplitude for 10-bit ADC
FS_amplitude = 2^(N_bits-1);

P1_dBFS = 20*log10(max(P1/FS_amplitude, eps));

% Frequency axis in MHz
f = (0:N/2)*Fs/N/1e6;

%% Plot frequency spectrum
figure('Color', 'w');

plot(f, P1_dBFS, 'b-', 'LineWidth', 1);
grid on;
box on;

xlabel('Frequency (MHz)', 'Color', 'k');
ylabel('Amplitude (dBFS)', 'Color', 'k');
title('ADC Output Frequency Spectrum', 'Color', 'k');

xlim([0 Fs/2/1e6]);
ylim([-120 10]);

set(gca, 'Color', 'w', ...
         'XColor', 'k', ...
         'YColor', 'k', ...
         'FontSize', 12);

%% Print FFT information
fprintf('\nFFT information:\n');
fprintf('FFT samples: %d\n', N);
fprintf('FFT bin spacing: %.3f kHz\n', Fs/N/1e3);
fprintf('Fundamental bin: %d\n', k);
fprintf('Fundamental frequency: %.6f MHz\n', f(k+1));

%% Calculate SNDR

% Exclude only DC and fundamental
exclude_sndr = false(N,1);

exclude_sndr(1) = true;
exclude_sndr([k+1, N-k+1]) = true;

noise_distortion_power = sum(power_spectrum(~exclude_sndr));

SNDR = 10*log10(signal_power/noise_distortion_power);

% Effective Number of Bits
ENOB = (SNDR - 1.76)/6.02;

fprintf('SNDR = %.2f dB\n', SNDR);
fprintf('ENOB = %.2f bits\n', ENOB);

%% Calculate SNR

% Calculate two-sided FFT power spectrum
power_spectrum = abs(Y/N).^2;

% Fundamental frequency bin (MATLAB indexing)
fund_bin = k + 1;

% Fundamental signal power (positive and negative frequencies)
signal_power = power_spectrum(fund_bin) + ...
               power_spectrum(N-k+1);

% Identify DC, fundamental, and harmonics
exclude_bins = false(N,1);

% Exclude DC
exclude_bins(1) = true;

% Exclude fundamental
exclude_bins([k+1, N-k+1]) = true;

% Exclude harmonics (2nd to 5th)
num_harmonics = 5;

for h = 2:num_harmonics
    harmonic_bin = mod(h*k, N);

    % Exclude positive and negative frequency components
    exclude_bins(harmonic_bin+1) = true;
    exclude_bins(mod(N-harmonic_bin,N)+1) = true;
end

% Calculate noise power
noise_power = sum(power_spectrum(~exclude_bins));

% Calculate SNR
SNR = 10*log10(signal_power/noise_power);

fprintf('\nADC Performance:\n');
fprintf('SNR = %.2f dB\n', SNR);