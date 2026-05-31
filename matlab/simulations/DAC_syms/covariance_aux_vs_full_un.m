%% Threshold Error Analysis for SAR ADC Calibration
% This script compares:
%   1. Ideal thresholds
%   2. Calibrated thresholds
%   3. Threshold error covariance
%   4. Correlation structure
%
% e_k = V_th,k - V_ideal,k

clear;
clc;
close all;

%% ADC PARAMETERS

N = 8;                      % ADC resolution
Vref = 1.0;

num_thresholds = 2^N - 1;

%% IDEAL THRESHOLDS

% Ideal transition levels centered between codes
Videal = ((1:num_thresholds) / 2^N) * Vref;

%% EXAMPLE 1:
% "Unary-like" calibration
% Small mostly uncorrelated errors

sigma_unary = 0.0002 * Vref;

Vth_unary = Videal + sigma_unary * randn(size(Videal));

%% EXAMPLE 2:
% "Auxiliary DAC-like" calibration
% Correlated low-rank residual errors

x = linspace(0,1,num_thresholds);

% Smooth correlated distortion shape
corr_shape = ...
      0.0015*sin(2*pi*x) ...
    + 0.0008*cos(4*pi*x);

% Small random component
noise = 0.0001 * randn(size(x));

Vth_aux = Videal + corr_shape + noise;

%% THRESHOLD ERRORS

e_unary = Vth_unary - Videal;
e_aux   = Vth_aux   - Videal;

%% RMS ERROR

rms_unary = rms(e_unary);
rms_aux   = rms(e_aux);

fprintf('Unary RMS threshold error: %e\n', rms_unary);
fprintf('Aux DAC RMS threshold error: %e\n', rms_aux);

%% PLOTS

figure;

subplot(2,1,1)
plot(e_unary,'LineWidth',1.5)
grid on
title('Unary CDAC Threshold Errors')
xlabel('Threshold Index')
ylabel('Error (V)')

subplot(2,1,2)
plot(e_aux,'LineWidth',1.5)
grid on
title('Auxiliary DAC Threshold Errors')
xlabel('Threshold Index')
ylabel('Error (V)')

%% ERROR AUTOCORRELATION
% Correlated errors generate harmonics more efficiently

figure;

[acf_unary,lags1] = xcorr(e_unary,'coeff');
[acf_aux,lags2]   = xcorr(e_aux,'coeff');

subplot(2,1,1)
plot(lags1,acf_unary,'LineWidth',1.5)
grid on
title('Unary CDAC Error Autocorrelation')
xlabel('Lag')
ylabel('Correlation')

subplot(2,1,2)
plot(lags2,acf_aux,'LineWidth',1.5)
grid on
title('Aux DAC Error Autocorrelation')
xlabel('Lag')
ylabel('Correlation')

%% FFT OF THRESHOLD ERRORS
% Shows harmonic structure in threshold placement

figure;

E1 = abs(fft(e_unary));
E2 = abs(fft(e_aux));

f = 0:length(E1)-1;

subplot(2,1,1)
stem(f,E1,'filled')
xlim([0 40])
grid on
title('Unary CDAC Error Spectrum')
xlabel('Spatial Frequency')
ylabel('|FFT|')

subplot(2,1,2)
stem(f,E2,'filled')
xlim([0 40])
grid on
title('Aux DAC Error Spectrum')
xlabel('Spatial Frequency')
ylabel('|FFT|')

%% OPTIONAL:
% Covariance matrices from Monte Carlo runs

num_mc = 500;

E_unary_mc = zeros(num_mc,num_thresholds);
E_aux_mc   = zeros(num_mc,num_thresholds);

for k = 1:num_mc

    % Unary-like
    E_unary_mc(k,:) = ...
        sigma_unary * randn(1,num_thresholds);

    % Aux-like correlated structure
    phase = 2*pi*rand;

    corr_shape = ...
        0.0015*sin(2*pi*x + phase) + ...
        0.0008*cos(4*pi*x + phase);

    noise = 0.0001 * randn(size(x));

    E_aux_mc(k,:) = corr_shape + noise;
end

%% Covariance matrices

C_unary = cov(E_unary_mc);
C_aux   = cov(E_aux_mc);

%% Plot covariance matrices

figure;

subplot(1,2,1)
imagesc(C_unary)
axis image
colorbar
title('Unary CDAC Error Covariance')

subplot(1,2,2)
imagesc(C_aux)
axis image
colorbar
title('Aux DAC Error Covariance')

%% Singular value decomposition
% Measures effective dimensionality of correction space

s_unary = svd(C_unary);
s_aux   = svd(C_aux);

figure;

semilogy(s_unary,'o-','LineWidth',1.5)
hold on
semilogy(s_aux,'o-','LineWidth',1.5)

grid on
xlabel('Singular Value Index')
ylabel('Singular Value')

legend('Unary','Aux DAC')
title('Effective Rank of Threshold Error Space')