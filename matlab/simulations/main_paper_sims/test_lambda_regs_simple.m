clear; close all; clc;
% Lambda regularization sweep for flash ADC dither calibration
Vhigh = 1;
Vlow = -1;
N_bits = 8;
cal_cycles = 1000000;
cal_len = 10000;
N = (2048*2^-3) - 1; % fft size
fs = 48000;  % coherent sampling
f0 = (13/N)*fs;
f1 = (15.24532/N)*fs;
non_lin_parameters = [0 1 0];
cal_cutoff = 0;
cal_constant = 0.01;

L = 2^N_bits;
LSB = (Vhigh - Vlow) / L;
Vinc = 1 * LSB;

lambda_regs = logspace(-3, 2, 6); % regularization parameter sweep

% Create ideal thresholds
ideal_thresholds = linspace(Vlow, Vhigh, L+1)';
ideal_thresholds = ideal_thresholds(2:end-1);

% Initialize noisy thresholds
init_thresholds = linspace(Vlow, Vhigh, L+1)';
init_thresholds(6) = init_thresholds(6) + LSB*0.9;
init_thresholds = init_thresholds(2:end-1);
noise_amp = 10 * LSB;
noisy = ideal_thresholds + noise_amp*(2*rand(size(ideal_thresholds))-1);
noisy = sort(noisy);
noisy = interp1(ideal_thresholds, noisy, ideal_thresholds, 'linear', 'extrap');
init_thresholds = noisy;

num_cases = length(lambda_regs);
SNDRs_cases = zeros(num_cases,1);
SNDRs_cases_max = zeros(num_cases,1);
converge_times = zeros(num_cases,1);

% Test and calibration inputs
t = (0:1/fs:(cal_len-1)/fs)';
analog_in = sin(2*pi*f1*t);
analog_in2 = sin(2*pi*f0*t);

[ideal_digi_out] = flash_adc(analog_in2, N_bits, Vhigh, Vlow, ideal_thresholds);
[initial_digi_out] = flash_adc(analog_in2, N_bits, Vhigh, Vlow, init_thresholds);
[post_calib_digi_out] = post_calib_flash(analog_in2, N_bits, Vhigh, Vlow, init_thresholds);

for i = 1:num_cases
    lambda_reg = lambda_regs(i);
    
    [digi_out, SNDRs, last_thresholds] = flash_adc_dither_sim_simple(analog_in, cal_len, cal_cycles, ...
        cal_constant, cal_cutoff, init_thresholds, Vhigh, Vlow, Vinc, N_bits, non_lin_parameters, lambda_reg, N, analog_in2);

    SNDRs_cases(i) = mean(SNDRs(end:end-500000));
    SNDRs_cases_max(i) = max(SNDRs);

    if max(SNDRs) <= SNDRs(1)
        converge_times(i) = 0;
    else
        SNDR_intercept = (max(SNDRs) - min(SNDRs))*0.95 + min(SNDRs);
        j = 1;
        done = false;
        while ~done && j <= cal_cycles
            if SNDRs(j) > SNDR_intercept
                converge_times(i) = j * cal_len;
                done = true;
            end
            j = j + 1;
        end
    end

    figure;
    plot(SNDRs);
    xlabel('Calibration cycle');
    ylabel('SNDR (dB)');
    title(sprintf('SNDR evolution for lambda\\_reg = %g', lambda_reg));
    grid on;

    figure;
    final_INL = (last_thresholds(1:end-1) - ideal_thresholds) / LSB;
    plot(final_INL);
    xlabel('Threshold index');
    ylabel('INL error (LSB)');
    title(sprintf('Final INL for lambda\\_reg = %g', lambda_reg));
    grid on;
end

figure;
semilogx(lambda_regs, SNDRs_cases, '*-', 'MarkerSize', 10);
hold on;
semilogx(lambda_regs, SNDRs_cases_max, 'o--', 'MarkerSize', 8);
yline(calculate_SNDR(initial_digi_out, analog_in2, N), '--r', 'Initial SNDR');
yline(calculate_SNDR(ideal_digi_out, analog_in2, N), '--g', 'Ideal SNDR');
yline(calculate_SNDR(post_calib_digi_out, analog_in2, N), '--m', 'Post-calibration SNDR');
xlabel('lambda\\_reg');
ylabel('SNDR (dB)');
title('SNDR vs lambda\\_reg');
grid on;
legend('Final SNDR', 'Max SNDR', 'Initial SNDR', 'Ideal SNDR', 'Post-calibration SNDR', 'Location', 'best');

figure;
loglog(lambda_regs, max(converge_times, 1), 's-', 'MarkerSize', 10);
xlabel('lambda\\_reg');
ylabel('Convergence time (samples)');
title('Convergence time vs lambda\\_reg');
grid on;
