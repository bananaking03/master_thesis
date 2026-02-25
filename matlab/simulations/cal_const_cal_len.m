clear; close all; clc;
% change to differantial system!!!!! Check distortie
Vhigh = 1;
Vlow = -1;
N_bits =10; % probeer meer bits
% cal_cycles = 5000;
cal_cycles = 50000;
N = (2048*2^-3) - 1; % fft size
fs = 48000;  % coherent sampling
f0 = (13/N)*fs;
% f0 = 990;
non_lin_parameters = [0 1 0];
cal_len = 40000;  % need increase for more bits
% N = 2048*2^4; % fft size
analyze_spesific = 1;
cal_cutoff = 0;

L = 2^N_bits;
Vinc = (Vhigh - Vlow) / L; % 1 LSB

% % % --- Create ideal thresholds ---
ideal_thresholds = linspace(Vlow, Vhigh, L+1)';
ideal_thresholds = ideal_thresholds(2:end-1);

% Initialize thresholds
init_thresholds = linspace(Vlow, Vhigh, L+1)'; % L+1 edges
init_thresholds(4) = init_thresholds(4) + Vinc*0.9;
init_thresholds(7) = init_thresholds(7) - Vinc*0.9;
init_thresholds = init_thresholds(2:end-1); % remove 0 and Vref

t = (0:1/fs:(cal_cycles-1)/fs)';
analog_in =sin(2*pi*f0*t);

delta_thr1 = zeros(cal_cycles,1);
delta_thr2 = zeros(cal_cycles,1);

thresholds = init_thresholds;

D = randi([0 1], cal_cycles,1);

% add dither
analog_in_inc = analog_in + D*Vinc;

% adc
adc_out = flash_adc(analog_in_inc,N_bits,Vhigh,Vlow,thresholds(1:end));

for i = 1:cal_cycles
    % create digital output
    %     digi_out((i-1)*cal_len+1:i*cal_len) = adc_out - D;
    out = adc_out(1:i) - D(1:i);
    
    % seperate incremented values and non-incremented values
    D_plus = out(D(1:i) == 1);
    D_min = out(D(1:i) == 0);
    
    % Compute histogram counts for each integer in range, ²one extra bin at
    % the end which gets discarded due to overflow
    edges = -0.5:1:(L + 0.5);   % L = 2^N_bits, so for 5-bit => -0.5:1:32.5

    H_plus = histcounts(D_plus, edges);
    H_min  = histcounts(D_min , edges);        % Hmin has bin 33????????????????????????????????
    
    binCenters = edges(1:end-1) + diff(edges)/2;
    
    H_delta = H_plus(1:L) - H_min(1:L);   % use bins 0..31, discard overflow (bin 32)

    delta_thr1(i) = H_delta(3);
    delta_thr2(i) = H_delta(6);
end
%% plots

figure;
plot(delta_thr1);
figure;
plot(delta_thr2);
%% calculate slopes
x = (1:cal_cycles)';

p1 = polyfit(x, delta_thr1, 1);
p2 = polyfit(x, delta_thr2, 1);

slope1 = p1(1);
slope2 = p2(1);

disp(['Slope delta\_thr1 = ', num2str(slope1)])
disp(['Slope delta\_thr2 = ', num2str(slope2)])