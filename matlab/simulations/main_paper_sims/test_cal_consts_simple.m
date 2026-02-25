clear; close all; clc;
% change to differantial system!!!!! Check distortie
Vhigh = 1;
Vlow = -1;
N_bits =10; % probeer meer bits
cal_cycles = 50000;
N = (2048*2^-3) - 1; % fft size
fs = 48000;  % coherent sampling
f0 = (13/N)*fs;
non_lin_parameters = [0 1 0];
cal_len = 400000;  % need increase for more bits
% N = 2048*2^4; % fft size
analyze_spesific = 1;
cal_cutoff = 0;

save_video = true;

L = 2^N_bits;
Vinc = (Vhigh - Vlow) / L; % 1 LSB

% cal_consts = [8*10^(-5) 9*10^(-5) 10^(-4) 2*10^(-4) 3*10^(-4) 4*10^(-4)];
% cal_consts = [10^(-3) 10^(-4) 10^(-5) 10^(-6) 10^(-7) 10^(-8)];
% cal_consts = logspace(-8, -4, 10)';
cal_consts = [10^-6];

to_dB = @(x) 20*log10(abs(x) + eps);  % Safe dB conversion

% % % --- Create ideal thresholds ---
ideal_thresholds = linspace(Vlow, Vhigh, L+1)';
ideal_thresholds = ideal_thresholds(2:end-1);

% Initialize thresholds
init_thresholds = linspace(Vlow, Vhigh, L+1)'; % L+1 edges
init_thresholds(4) = init_thresholds(4) + Vinc*0.9;
init_thresholds = init_thresholds(2:end-1); % remove 0 and Vref
% noise added around each threshold independently
noise_amp = 40 * Vinc;
noisy = ideal_thresholds + noise_amp*(2*rand(size(ideal_thresholds))-1);

% but ensure thresholds remain sorted and equally spaced on average
noisy = sort(noisy);

% now RE-MAP them to preserve spacing
noisy = interp1(ideal_thresholds, noisy, ideal_thresholds, 'linear', 'extrap');

init_thresholds = noisy;

% Preallocate storage
num_cases = length(cal_consts);
thresholds = zeros(L-1, (cal_cycles+1)*num_cases);
SNDRs_cases = zeros(num_cases,1);
SNDRs = zeros(cal_cycles,1);

% analog_in = 0.5 + 0.5*sin(2*pi*20*t);
t = (0:1/fs:(cal_len-1)/fs)';
analog_in =sin(2*pi*f0*t);
% analog_in = mod(t,1);

[ideal_digi_out] = flash_adc(analog_in,N_bits,Vhigh,Vlow,ideal_thresholds);

[initial_digi_out] = flash_adc(analog_in,N_bits,Vhigh,Vlow,init_thresholds);

% Run calibration for each calibration length
for i = 1:num_cases
    cal_constant = cal_consts(i);
    
    % Call your calibration function
    [digi_out, SNDRs] = flash_adc_dither_sim_simple(analog_in, cal_len, cal_cycles, ...
        cal_constant, cal_cutoff, init_thresholds, Vhigh, Vlow, Vinc, N_bits, non_lin_parameters,N);

    figure;
    plot(SNDRs);

    SNDRs_cases(i) = SNDRs(cal_cycles); 
end

figure;
semilogx(cal_consts, SNDRs_cases, '*');
hold on;
yline(calculate_SNDR(initial_digi_out,analog_in,N), '--r', 'initial SNDR');
yline(calculate_SNDR(ideal_digi_out,analog_in,N), '--r', 'ideal SNDR');
xlabel('cal constant');
ylabel('SNDR achieved in last iteration (dB)');
title('SNDR for each cal constant');
grid on;
hold off;


