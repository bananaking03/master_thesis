
%% ADC output analysis
clearvars -except results_file

Nfft = 1000;
N_bits = 10;
f0 = 10e6;

% results_file = 'results_done_matdeltareg_dis_va.csv';
% results_file = 'results_done_dis_va.csv';
results_file = 'results_done_ideal_va.csv';

[results, timestamps] = read_histogram(results_file);

% flip collumns
results(:,:) = results(:,[2 1]);

%% Extract output

out = double(results(1:2:end,1));
time = timestamps(1:2:end)*1e-9;
out = out(1:1000);
time = time(1:1000);
timestamps = timestamps(1:2000);

% Convert to signed 10-bit two's complement
out(out >= 2^(N_bits-1)) = ...
    out(out >= 2^(N_bits-1)) - 2^N_bits;

% Remove DC offset for spectral analysis
out_ac = out - mean(out);

%% Sampling information

dt = mean(diff(time));
Fs = 1/dt;

N = length(out);

% Use the actual number of samples if Nfft is too small
if Nfft < N
    warning('Nfft is smaller than the number of samples. Increasing Nfft.');
    Nfft = 2^nextpow2(N);
end

fprintf('\n========== ADC INFORMATION ==========\n');
fprintf('Number of samples: %d\n', N);
fprintf('Sampling frequency: %.6f MHz\n', Fs/1e6);
fprintf('Sampling period: %.6f ns\n', dt*1e9);
fprintf('Observation time: %.6f us\n', (time(end)-time(1))*1e6);
fprintf('Number of input cycles: %.8f\n', N*f0/Fs);

%% Time-domain measurements

DC_offset = mean(out);
RMS_value = rms(out_ac);
Peak_value = max(abs(out_ac));
Peak_to_peak = max(out) - min(out);

% Check for ADC clipping
max_code = 2^(N_bits-1)-1;
min_code = -2^(N_bits-1);

clipped_samples = sum(out >= max_code | out <= min_code);

fprintf('\n========== TIME DOMAIN ==========\n');
fprintf('DC offset: %.4f LSB\n', DC_offset);
fprintf('AC RMS: %.4f LSB\n', RMS_value);
fprintf('Peak amplitude: %.4f LSB\n', Peak_value);
fprintf('Peak-to-peak: %.4f LSB\n', Peak_to_peak);
fprintf('Clipped samples: %d (%.4f%%)\n', ...
    clipped_samples, 100*clipped_samples/N);

%% FFT
N = length(out);

Fs = 1/mean(diff(time));
f0 = 110/1000 * Fs;

Y = fft(out - mean(out), N);

frequency = (0:N-1)*Fs/N;
magnitudeSpectrum = abs(Y)/N;

%% Single-sided power spectrum
P2 = abs(Y).^2/N^2;
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

frequency_single = (0:N/2)*Fs/N;

%% SNDR calculation
fundamental_bin = 111;
signal_power = P1(fundamental_bin);

harmonic_bins = [];

for h = 2:floor((Fs/2)/f0)
    bin = round(h*f0/(Fs/N))+1;

    if bin <= length(P1)
        harmonic_bins(end+1) = bin;
    end
end

all_bins = 1:length(P1);
noise_distortion_bins = setdiff(all_bins, [1 fundamental_bin]);
noise_bins = setdiff(noise_distortion_bins, harmonic_bins);

noise_power = sum(P1(noise_bins));
distortion_power = sum(P1(harmonic_bins));

SNDR = 10*log10(signal_power/(noise_power+distortion_power));
SNR = 10*log10(signal_power/noise_power);
THD = 10*log10(distortion_power/signal_power);

fprintf('\nSNDR: %.2f dB\n', SNDR);
fprintf('SNR:  %.2f dB\n', SNR);
fprintf('THD:  %.2f dB\n', THD);
fprintf('ENOB: %.2f bits\n', (SNDR-1.76)/6.02);

%% Plot time-domain signal
figure;
plot(timestamps(1:2:end), out, 'o-');
grid on;
xlabel('Time');
ylabel('ADC output');
title('ADC output');
% xlim([12 112])

%% Plot full FFT
figure;
plot(frequency, magnitudeSpectrum);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('ADC output spectrum');
xlim([0 Fs/2]);

%% Plot single-sided spectrum in dB
figure;
plot(frequency_single, 10*log10(P1 + eps));
grid on;
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Single-sided power spectrum');
xlim([0 Fs/2]);


%% IDEAL ADC COMPARISON
% Generate a perfect 10-bit ADC output using the same sampling conditions

fprintf('\n\n========== IDEAL ADC COMPARISON ==========\n');

%% Generate ideal input

ideal_time = (0:N-1)/Fs;

% Coherent input frequency: 110 cycles over N samples
ideal_f0 = 110*Fs/N;

% Full-scale sine wave for a signed 10-bit ADC
ideal_input = (2^(N_bits-1)-1)*sin(2*pi*ideal_f0*ideal_time);

% Ideal quantization
ideal_out = round(ideal_input);

% Saturate to the signed 10-bit range
ideal_out = max(min(ideal_out, 2^(N_bits-1)-1), ...
    -2^(N_bits-1));

%% Ideal ADC FFT

ideal_Y = fft(ideal_out - mean(ideal_out), N);

ideal_frequency = (0:N-1)*Fs/N;
ideal_magnitudeSpectrum = abs(ideal_Y)/N;

% Single-sided power spectrum
ideal_P2 = abs(ideal_Y).^2/N^2;
ideal_P1 = ideal_P2(1:N/2+1);
ideal_P1(2:end-1) = 2*ideal_P1(2:end-1);

ideal_frequency_single = (0:N/2)*Fs/N;

%% Ideal SNDR calculation

ideal_fundamental_bin = round(ideal_f0/(Fs/N))+1;

ideal_signal_power = ideal_P1(ideal_fundamental_bin);

ideal_all_bins = 1:length(ideal_P1);

ideal_noise_distortion_bins = setdiff(ideal_all_bins, ...
    [1 ideal_fundamental_bin]);

ideal_noise_distortion_power = ...
    sum(ideal_P1(ideal_noise_distortion_bins));

ideal_SNDR = 10*log10(ideal_signal_power / ...
    ideal_noise_distortion_power);

ideal_ENOB = (ideal_SNDR-1.76)/6.02;

%% Ideal ADC results

fprintf('Number of samples: %d\n', N);
fprintf('Input frequency: %.6f MHz\n', ideal_f0/1e6);
fprintf('Number of input cycles: %.0f\n', 110);
fprintf('Fundamental bin: %d\n', ideal_fundamental_bin);

fprintf('\nIdeal ADC SNDR: %.2f dB\n', ideal_SNDR);
fprintf('Ideal ADC ENOB: %.2f bits\n', ideal_ENOB);

%% Comparison

fprintf('\n========== SNDR COMPARISON ==========\n');
fprintf('Simulated ADC SNDR: %.2f dB\n', SNDR);
fprintf('Ideal ADC SNDR:     %.2f dB\n', ideal_SNDR);
fprintf('Difference:         %.2f dB\n', ideal_SNDR-SNDR);

%% Plot ideal time-domain signal

figure;
plot(ideal_time*1e6, ideal_out, 'o-');
grid on;
xlabel('Time (us)');
ylabel('ADC output (LSB)');
title('Ideal 10-bit ADC output');
set(gca, 'Color', 'w');
set(gcf, 'Color', 'w');
% xlim([0 0.1])

%% Plot ideal FFT

figure;
plot(ideal_frequency, ideal_magnitudeSpectrum);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Ideal ADC output spectrum');
xlim([0 Fs/2]);
set(gca, 'Color', 'w');
set(gcf, 'Color', 'w');

%% Plot ideal single-sided spectrum

figure;
plot(ideal_frequency_single, 10*log10(ideal_P1 + eps));
grid on;
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('Ideal ADC spectrum - SNDR = %.2f dB', ideal_SNDR));
xlim([0 Fs/2]);
set(gca, 'Color', 'w');
set(gcf, 'Color', 'w');

%% different way 
SNDR = calculate_SNDR(out_ac,ideal_input,N)