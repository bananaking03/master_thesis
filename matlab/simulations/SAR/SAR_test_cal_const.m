clear; close all; clc;
% change to differantial system!!!!! Check distortie
Vhigh = 1;
Vlow = -1;
N_bits =10; % probeer meer bits
N_bits_caps = 15; 
cal_cycles = 1;
N = (2048*2^-3) - 1; % fft size
fs = 48000;  % coherent sampling
f0 = (13/N)*fs;
f1 = (15.24532/N)*fs;
LSC = 10^-15; % least signifigant capacitance
mismatch_LSC = 0.1;
non_lin_parameters = [0 1 0];
% cal_len = 5000000;  % need increase for more bits
cal_len = 1000;
% cal_len = 10001
% N = 2048*2^4; % fft size
analyze_spesific = 1;
cal_cutoff = 0;

save_video = true;

L = 2^N_bits;
Vinc = (Vhigh - Vlow) / L; % 1 LSB

% cal_consts = [8*10^(-5) 9*10^(-5) 10^(-4) 2*10^(-4) 3*10^(-4) 4*10^(-4)];
% cal_consts = [10^(-3) 10^(-4) 10^(-5) 10^(-6) 10^(-7) 10^(-8)];
% cal_consts = logspace(-12, -3, 10)';
cal_consts = [1*10^-8];

to_dB = @(x) 20*log10(abs(x) + eps);  % Safe dB conversion

% % % --- Create ideal thresholds ---
ideal_thresholds = linspace(Vlow, Vhigh, L+1)';
ideal_thresholds = ideal_thresholds(2:end-1);

% Create the caps
caps = derrive_caps(N_bits_caps,mismatch_LSC,LSC);

% thresholds from caps
init_thresholds = cap_to_thr(caps(1:N_bits+1),Vlow,Vhigh);
pos_thresholds = cap_to_thr(caps,Vlow,Vhigh); % possible thresholds
[~, idx] = min(abs(pos_thresholds(:) - ideal_thresholds(:).'), [], 1);
best_thresholds = pos_thresholds(idx);

% Preallocate storage
num_cases = length(cal_consts);
thresholds = zeros(L-1, (cal_cycles+1)*num_cases);
SNDRs_cases = zeros(num_cases,1);
SNDRs = zeros(cal_cycles,1);

% analog_in = 0.5 + 0.5*sin(2*pi*20*t);
t = (0:1/fs:(cal_len-1)/fs)';
analog_in =sin(2*pi*f1*t);
% analog_in = mod(t,1);

t = (0:1/fs:(cal_len-1)/fs)';
analog_in2 =sin(2*pi*f0*t);
%% 

[ideal_digi_out] = flash_adc(analog_in2,N_bits,Vhigh,Vlow,ideal_thresholds);

[best_digi_out] = flash_adc(analog_in2,N_bits,Vhigh,Vlow,best_thresholds);

[initial_digi_out] = flash_adc(analog_in2,N_bits,Vhigh,Vlow,init_thresholds);
%% 

% Run calibration for each calibration constant
for i = 1:num_cases
    cal_constant = cal_consts(i);
    
    % Call your calibration function
    [digi_out, SNDRs] = SAR_adc_dither_sim(analog_in, cal_len, cal_cycles, ...
        cal_constant, cal_cutoff, init_thresholds(2:end), pos_thresholds, Vhigh, Vlow, Vinc, N_bits, non_lin_parameters,N, analog_in2);

    figure;
    plot(SNDRs);

    SNDRs_cases(i) = SNDRs(cal_cycles); 
end
%% 

% Plot SNDRs
figure;
semilogx(cal_consts, SNDRs_cases, '*');
hold on;
yline(calculate_SNDR(initial_digi_out,analog_in2,N), '--r', 'initial SNDR');
yline(calculate_SNDR(ideal_digi_out,analog_in2,N), '--r', 'ideal SNDR');
yline(calculate_SNDR(best_digi_out,analog_in2,N), '--r', 'best SNDR');
xlabel('cal constant');
ylabel('SNDR achieved in last iteration (dB)');
title('SNDR for each cal constant');
grid on;
hold off;