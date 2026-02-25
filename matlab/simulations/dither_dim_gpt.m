% ---------------- PARAMETERS (keep from your script) ----------------
L = 2^3;
Vref = 1;
Vinc = Vref / L; % 1 LSB
fs = 48000;
x_max = 0.01;
f0 = 1000;

% ... (threshold generation stays the same) ...
% init_thresholds = noisy; (your noisy thresholds)
% ideal_thresholds = ideal_thresholds;

% time & clean analog input
x = (0:1/fs:x_max)';                      % column
analog_in = 0.5 + 0.5*sin(2*pi*f0*x);     % between 0 and 1

N = length(analog_in);

% ---------------- PROPER DITHER (triangular, zero-mean) ----------------
% Triangular dither ~ sum of two independent uniform(-0.5,0.5) * 1LSB
u1 = rand(N,1) - 0.5;
u2 = rand(N,1) - 0.5;
dither_volts = (u1 + u2) * Vinc;   % triangular, range roughly [-1*Vinc, +1*Vinc], zero-mean

% Apply dither to input amplitude (NOT to time)
analog_in_dithered = analog_in + dither_volts;

% Make sure we don't exceed ADC range [0 Vref]
analog_in_dithered = min(max(analog_in_dithered, 0), Vref);

% Send through ADC (use your flash_adc function)
output_analog = flash_adc(analog_in, 3, Vref, 0, init_thresholds);            % no dither
output_analog_dith = flash_adc(analog_in_dithered, 3, Vref, 0, init_thresholds); % with dither

% ---------------- OPTIONAL: exact integer-LSB dithers and compensation ----------
% If you instead add deterministic integer-LSB dither (e.g. add 0 or 1 LSB),
% then you can remove it by subtracting that integer value:
% Example integer dither bits:
% Dbits = randi([0 1],N,1);                   % random 0/1 LSB add
% analog_inc = analog_in + Dbits * Vinc;      % add exactly 1 LSB when D=1
% output_inc = flash_adc(analog_inc,3,Vref,0,init_thresholds);
% output_comp = output_inc - Dbits;           % subtract the integer LSB in code domain
% But note: this integer approach is not triangular and may leave harmonics
% unless Dbits are truly random and uncorrelated with input.

% ---------------- PLOT & FFT CHECK ----------------
% center outputs before FFT (subtract midscale 4)
out_no_dith_centered = double(output_analog) - 4;
out_dith_centered = double(output_analog_dith) - 4;

% window
w = hann(N);
X_no_dith = fft(out_no_dith_centered .* w);
X_dith = fft(out_dith_centered .* w);

% single-sided
f = (0:N-1)*(fs/N);
half = 1:floor(N/2);
f_half = f(half);

% amplitude correction using sum(w)
X_no_mag = abs(X_no_dith(half))/sum(w);
X_dith_mag   = abs(X_dith(half))/sum(w);

figure;
semilogy(f_half, X_no_mag); hold on;
semilogy(f_half, X_dith_mag, '--');
xlim([0 8000]);
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Single-sided Spectra (proper triangular dither)');
legend('No dither','Triangular dither','Location','best');
grid on;

% mark fundamentals/harmonics for visual inspection
harmonics = [1 2 3 5] * f0;
for k=1:length(harmonics)
    [~, idx]=min(abs(f_half - harmonics(k)));
    plot(f_half(idx), X_no_mag(idx), 'ro');  % mark on no-dith curve
    plot(f_half(idx), X_dith_mag(idx), 'kx');% mark on dithered curve
end
