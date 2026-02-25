N_bits = 3;
L = 2^N_bits;
Vref = 1;
Vinc = Vref / L; % 1 LSB
fs = 48000;
x_max = 0.01;
f0 = 1000;

to_dB = @(x) 20*log10(abs(x) + eps);  % Safe dB conversion

init_thresholds = linspace(0, Vref, L+1)'; % L+1 edges
init_thresholds(4) = init_thresholds(4) + Vinc*0.9;
init_thresholds = init_thresholds(2:end-1); % remove 0 and Vref

% ---------------------------------noisy thresholds---------------------------------------

% % % --- Create ideal thresholds ---
ideal_thresholds = linspace(0, Vref, L+1)';
ideal_thresholds = ideal_thresholds(2:end-1);
% 
% % noise added around each threshold independently
% noise_amp = 100 * Vinc;
% noisy = ideal_thresholds + noise_amp*(2*rand(size(ideal_thresholds))-1);
% 
% % but ensure thresholds remain sorted and equally spaced on average
% noisy = sort(noisy);
% 
% % now RE-MAP them to preserve spacing
% noisy = interp1(ideal_thresholds, noisy, ideal_thresholds, 'linear', 'extrap');
% 
% init_thresholds = noisy;
% --------------------------------------------------------------------------
    
t = 0:1/fs:1;
% add dither
t_dith = t + Vinc;

x = (0:1/fs:x_max)';
analog_in = 0.5 + 0.5*sin(2*pi*f0*x);
D = randi([0 1], x_max*fs+1,1);
analog_in_inc = analog_in + D*Vinc;

output = flash_adc(t,N_bits,Vref,0,init_thresholds);
output_id = flash_adc(t,N_bits,Vref,0,ideal_thresholds);
output_dith = flash_adc(t_dith,N_bits,Vref,0,init_thresholds);
output_dith_comp = output_dith - 1;

% test signal without dither
output_analog = flash_adc(analog_in,N_bits,Vref,0,init_thresholds);

% test signal with dither
output_analog_dith = flash_adc(analog_in_inc,N_bits,Vref,0,init_thresholds); % Currently gives wrong result because adc returns [0 7]
output_analog_comp_dith = output_analog_dith - D;

%  ideal output
output_ideal = flash_adc(analog_in,N_bits,Vref,0,ideal_thresholds);

figure;
plot(t, output, 'LineWidth', 2)
hold on
plot(t,output_dith_comp, 'g', 'LineWidth', 2, 'LineStyle','-.')
plot(t,t*8, 'r', 'LineWidth', 2, 'LineStyle','-.')

% Add legend
legend('Original output', 'Dithered output', 'Ideal output')
xlabel('input')
ylabel('ADC output code')
title('Flash ADC Output with Dithering')
grid on

%% Plot frequency response
figure;
N = length(output_analog_comp_dith);            % number of samples
X_out_dith = fft(output_analog_comp_dith-L/2);               % compute FFT

% Frequency vector
f = (0:N-1)*(fs/N);

% Plot magnitude spectrum
plot(f(1:N/2), to_dB(X_out_dith(1:N/2)), 'Marker','.');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Spectrum');

hold on

% plot original frequency of input signal
X_ideal = fft(output_ideal-L/2);
plot(f(1:N/2), to_dB(X_ideal(1:N/2)), 'Marker','.');

% plot output without dither
X_out = fft(output_analog-L/2);
plot(f(1:N/2), to_dB(X_out(1:N/2)), 'Marker','.');
xlim([0 0.5e4]);

% Add legend
legend('Output with dither (compensated)', ...
       'Ideal output', ...
       'Output without dither');