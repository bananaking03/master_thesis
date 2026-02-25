%% Distortion Consequences Visualization (Full Upgraded Version)
clear; close all; clc;

fs = 48000;              % Sampling rate
t = 0:1/fs:0.01;         % Short time window (10 ms)

%% Common FFT settings
Nfft   = 2048;               % Zero-padding (8×) for smoother FFT
t = (0:Nfft-1)/fs;         % Short time window (10 ms)
win    = hann(Nfft)';            % Hann window
Wnorm  = sum(win) / Nfft;        % Normalize amplitude for window loss
f      = (0:Nfft-1)*(fs/Nfft);
dbFloor = -80;                % Baseline level for all spectra
to_dB = @(x) 20*log10(abs(x) + eps);

% Universal FFT helper (window + padding + dB + baseline)
% fft_db = @(x) max( to_dB( fft(x .* win, Nfft) / Wnorm ), dbFloor );
fft_db = @(x) max( to_dB( fft(x , Nfft) ), dbFloor );


%% ========================================================================
% 1. REDUCED DYNAMIC RANGE (Clipping)
% ========================================================================
x = 0.8*sin(2*pi*1000*t) + 0.3*sin(2*pi*300*t);  % Original signal
clipLevel = 0.5;
x_clip = max(min(x,clipLevel), -clipLevel);     % Hard clipping

% FFTs
X  = fft_db(x);
Xc = fft_db(x_clip);

figure;
subplot(2,1,1)
plot(t,x,'b'); hold on;
plot(t,x_clip,'r');
title('Reduced Dynamic Range (Clipping)');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Original','Clipped');
grid on;

subplot(2,1,2)
bar(f(1:Nfft/2), X(1:Nfft/2), 0.05, 'FaceColor','b','EdgeColor','b', 'LineWidth', 0.8,'BaseValue',-80); hold on;
bar(f(1:Nfft/2), Xc(1:Nfft/2),0.05,'FaceColor','r','EdgeColor','r', 'LineWidth', 0.1,'FaceAlpha',0.7,'BaseValue',-80);
title('Frequency Domain – Clipping Creates Spectral Spreading');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
xlim([0 5000]); ylim([dbFloor 40]);
legend('Original','Clipped');
grid on;


%% ========================================================================
% 2. HARMONIC DISTORTION
% ========================================================================
f0 = (128/2048)*fs;
x = sin(2*pi*f0*t);               % Fundamental tone
nonlinear = @(u) u + 0.3*u.^3;    % Cubic nonlinearity
x_harm = nonlinear(x);            % Distorted signal

% FFTs
X  = fft_db(x);
Xh = fft_db(x_harm);

figure;
subplot(2,1,1)
plot(t,x,'b'); hold on;
plot(t,x_harm,'r');
title('Harmonic Distortion');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Original','Distorted');
grid on;

subplot(2,1,2)
bar(f(1:Nfft/2), X(1:Nfft/2),0.05,'FaceColor','b','EdgeColor','b', 'LineWidth', 0.8,'BaseValue',-80); hold on;
bar(f(1:Nfft/2), Xh(1:Nfft/2),0.05,'FaceColor','r','EdgeColor','r', 'LineWidth', 0.1,'FaceAlpha',0.7,'BaseValue',-80);
title('Frequency Domain – Harmonics at 2f, 3f, …');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
% xlim([0 5000]); 
ylim([dbFloor max([40 X Xh])]);
legend('Original','Distorted');
grid on;


%% ========================================================================
% 3. INTERMODULATION DISTORTION (Two-tone test)
% ========================================================================
f1 = 1000; 
f2 = 1200;

x = 0.8*sin(2*pi*f1*t) + 0.8*sin(2*pi*f2*t);   % Two sinewaves
x_IMD = nonlinear(x);                          % Apply nonlinearity

% FFTs
X  = fft_db(x);
Xi = fft_db(x_IMD);

figure;
subplot(2,1,1)
plot(t, x,'b'); hold on;
plot(t, x_IMD,'r');
title('Intermodulation Distortion');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Original','Distorted');
grid on;

subplot(2,1,2)
bar(f(1:Nfft/2), X(1:Nfft/2),0.05,'FaceColor','b','EdgeColor','b', 'LineWidth', 0.8,'BaseValue',-80); hold on;
bar(f(1:Nfft/2), Xi(1:Nfft/2),0.05,'FaceColor','r','EdgeColor','r', 'LineWidth', 0.1,'FaceAlpha',0.7,'BaseValue',-80);
title('Frequency Domain – IMD at |m f1 ± n f2|');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
xlim([0 5000]); ylim([dbFloor 40]);
legend('Original','Distorted');
grid on;
