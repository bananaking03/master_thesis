%% Distortion Consequences Visualization
clear; close all; clc;

fs = 48000;              % Sampling rate
t = 0:1/fs:0.01;         % Short time window for display

to_dB = @(x) 20*log10(abs(x) + eps);  % Safe dB conversion

%% ------------------------------------------------------------
% 1. REDUCED DYNAMIC RANGE (Clipping)
%% ------------------------------------------------------------
x = 0.8*sin(2*pi*1000*t) + 0.3*sin(2*pi*300*t);   % Original signal
clipLevel = 0.5;
x_clip = max(min(x,clipLevel), -clipLevel);        % Hard clipping

% FFT
N = length(x);
f = (0:N-1)*(fs/N);
X  = abs(fft(x));
Xc = abs(fft(x_clip));

figure;
subplot(2,1,1)
plot(t,x, 'b'); hold on;
plot(t,x_clip,'r');
title('Reduced Dynamic Range (Clipping)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original','Clipped');

subplot(2,1,2)
plot(f(1:N/2), to_dB(X(1:N/2)), 'b', 'LineStyle', 'none', 'Marker','.'); hold on;
plot(f(1:N/2), to_dB(Xc(1:N/2)), 'r', 'LineStyle', 'none','Marker','.');
title('Frequency Domain – Clipping creates spectral spreading');
xlabel('Frequency (Hz)');
xlim([0 5000]);
ylim([-25 50]);
ylabel('Magnitude (dB)');
legend('Original','Clipped');


%% ------------------------------------------------------------
% 2. HARMONIC DISTORTION
%% ------------------------------------------------------------
f0 = 1000;                         % Fundamental tone
x = sin(2*pi*f0*t);                % Clean tone
nonlinear = @(u) u + 0.3*u.^3;     % Cubic nonlinearity
x_harm = nonlinear(x);             % Harmonic distortion

% FFT
X  = abs(fft(x));
Xh = abs(fft(x_harm));

figure;
subplot(2,1,1)
plot(t, x, 'b'); hold on;
plot(t, x_harm,'r');
title('Harmonic Distortion (Nonlinear System)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original','Distorted');

subplot(2,1,2)
plot(f(1:N/2), to_dB(X(1:N/2)), 'b', 'LineStyle', 'none', 'Marker','.'); hold on;
plot(f(1:N/2), to_dB(Xh(1:N/2)), 'r', 'LineStyle', 'none', 'Marker','.');
title('Frequency Domain – Harmonics appear at 2f, 3f, etc');
xlabel('Frequency (Hz)');
xlim([0 5000]);
ylim([-25 50]);
ylabel('Magnitude (dB)');
legend('Original','Distorted');


%% ------------------------------------------------------------
% 3. INTERMODULATION DISTORTION
%% ------------------------------------------------------------
f1 = 1000; 
f2 = 1200;
x = 0.8*sin(2*pi*f1*t) + 0.8*sin(2*pi*f2*t);   % Two-tone test signal
x_IMD = nonlinear(x);                          % Apply nonlinearity

% FFT
X  = abs(fft(x));
Xi = abs(fft(x_IMD));

figure;
subplot(2,1,1)
plot(t, x, 'b'); hold on;
plot(t, x_IMD,'r');
title('Intermodulation Distortion');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original','Distorted');

subplot(2,1,2)
plot(f(1:N/2), to_dB(X(1:N/2)), 'b','LineStyle', 'none','Marker','.'); hold on;
plot(f(1:N/2), to_dB(Xi(1:N/2)), 'r','LineStyle', 'none','Marker','.');
title('Frequency Domain – IMD products appear at |mf1 ± nf2|');
xlabel('Frequency (Hz)');
xlim([0 5000]);
ylim([-25 50]);
ylabel('Magnitude (dB)');
legend('Original','Distorted');
