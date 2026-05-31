clear; close all; clc;
% change to differantial system!!!!! Check distortie
Vhigh = 1;
Vlow = -1;
N_bits =8; % probeer meer bits
cal_cycles =  1000000;
N = (2048*2^-3) - 1; % fft size
fs = 48000;  % coherent sampling
f0 = (13/N)*fs;
f1 = (15.24532/N)*fs;
non_lin_parameters = [0 1 0];
% cal_len = 10000000/(2^5);  % need increase for more bits
% cal_len = 10000
% cal_len = 10001
% N = 2048*2^4; % fft size
analyze_spesific = 1;
cal_cutoff = 0;
cal_constant = 0.1;
cal_len = 5000;

save_video = true;

L = 2^N_bits;
LSB = (Vhigh - Vlow) / L;
Vincs = linspace(0.5,1.5,9) * LSB;

to_dB = @(x) 20*log10(abs(x) + eps);  % Safe dB conversion

% % % --- Create ideal thresholds ---
ideal_thresholds = linspace(Vlow, Vhigh, L+1)';
ideal_thresholds = ideal_thresholds(2:end-1);

% Initialize thresholds
init_thresholds = linspace(Vlow, Vhigh, L+1)'; % L+1 edges
init_thresholds(6) = init_thresholds(6) + LSB*0.9;
init_thresholds = init_thresholds(2:end-1); % remove 0 and Vref
% noise added around each threshold independently
noise_amp = 10 * LSB;
noisy = ideal_thresholds + noise_amp*(2*rand(size(ideal_thresholds))-1);

% but ensure thresholds remain sorted and equally spaced on average
noisy = sort(noisy);

% now RE-MAP them to preserve spacing
noisy = interp1(ideal_thresholds, noisy, ideal_thresholds, 'linear', 'extrap');

init_thresholds = noisy;

% Preallocate storage
num_cases = length(Vincs);
% thresholds = zeros(L-1, (cal_cycles+1)*num_cases);
SNDRs_cases = zeros(num_cases,1);
SNDRs_cases_max = zeros(num_cases,1);
converge_times = zeros(num_cases,1);
% final_INL = zeros(num_cases,1);

% analog_in = 0.5 + 0.5*sin(2*pi*20*t);
t = (0:1/fs:(cal_len-1)/fs)';
analog_in =sin(2*pi*f1*t);
% analog_in = mod(t,1);

t = (0:1/fs:(cal_len-1)/fs)';
analog_in2 =sin(2*pi*f0*t);

% Run calibration for each calibration length
%% 
for i = 1:num_cases
    % cal_cycles = cal_cycless(i);
%%
    Vinc = Vincs(i);
    SNDRs = zeros(cal_cycles,1);
    
    [ideal_digi_out] = flash_adc(analog_in2,N_bits,Vhigh,Vlow,ideal_thresholds);
    
    [initial_digi_out] = flash_adc(analog_in2,N_bits,Vhigh,Vlow,init_thresholds);
    
    [post_calib_digi_out] = post_calib_flash(analog_in2,N_bits,Vhigh,Vlow,init_thresholds);
    
    % Call your calibration function
    [digi_out, SNDRs,last_thresholds] = flash_adc_dither_sim_simple(analog_in, cal_len, cal_cycles, ...
        cal_constant, cal_cutoff, init_thresholds, Vhigh, Vlow, Vinc, N_bits, non_lin_parameters,N, analog_in2);

    figure;
    plot(SNDRs);
    hold on
    xlabel('calibration cycle')
    ylabel('SNDR (dB)')

    SNDRs_cases(i) = mean(SNDRs(end-30000:end)); 
    SNDRs_cases_max(i) = max(SNDRs); 

    figure;
    final_INL = (last_thresholds(1:end-1)-ideal_thresholds)/LSB;
    plot(final_INL);
    xlabel('threshold')
    ylabel('INL error')

    % Calculate convergence time of SNDR
    if max(SNDRs) <= SNDRs(1)
        converge_times(i) = 0;
    else
        SNDR_intercept = (max(SNDRs) - min(SNDRs))*0.95 + min(SNDRs);
        j = 1;
        done = false;
        while (done == false)
            if (SNDRs(j) > SNDR_intercept)
                converge_times(i) = j*cal_len;
                done = true;
            end
            j = j+1;
        end
    end
    %%
end

%% Ideal SNDR
SNDR_ideal = 6.02*N_bits + 1.76;

% %% Compute SNDR
% N_vals = logspace(2, 6, 200);   % cal_len range
% scaling = (12 * cal_constant * 2^(N_bits)) / (2*(Vhigh-Vlow));
% SNDR = SNDR_ideal - 10*log10(1 + scaling ./ N_vals);
% SNDR_fit = SNDR_ideal - 10*log10(1 + 16.*scaling ./ N_vals);
% 
% %% Compute Convergence time
% converge_times_math = 1.5*N_vals*L/(cal_constant);
% % converge_times_math_fit = 10*N_vals*L/(cal_constant);
% % converge_times_math_fit = 100*N_vals*L;
% converge_times_math_fit = N_vals*log(0.05)/(2*log(1-cal_constant/L));

%% SNDR plot
figure;
yyaxis left
plot(linspace(-0.5,0.5,9), SNDRs_cases, '*','MarkerSize',12, 'Color',[0 0 1]);
hold on
plot(linspace(-0.5,0.5,9), SNDRs_cases_max, '*','MarkerSize',12,'Color',[0.7 0 1 0.3]);
% semilogx(N_vals, SNDR, 'LineWidth', 2,'Color',[0 0 1 0.3]);
% semilogx(N_vals, SNDR_fit, 'LineWidth', 2,'Color',[0 0 1]);
yline(calculate_SNDR(initial_digi_out,analog_in2,N), '--r', 'initial SNDR');
yline(calculate_SNDR(ideal_digi_out,analog_in2,N), '--r', 'ideal SNDR');
yline(calculate_SNDR(post_calib_digi_out,analog_in2,N), '--r', 'post callibration maximum SNDR');
xlabel('Dither mismatch');
ylabel('SNDR achieved in last iteration (dB)');
title('SNDR and convergency time for each cal length');
grid on;

yyaxis right
plot(linspace(-0.5,0.5,9), converge_times, '.','MarkerSize',12);
% semilogx(N_vals, converge_times_math, 'Linewidth',2,'Color',[255/255 165/255 0 0.3]);
% semilogx(N_vals, converge_times_math_fit, 'Linewidth',2,'Color',[255/255 165/255 0]);
ylabel('convege time (95% of final SNDR)');
set(gca, 'YScale', 'log')   % <-- THIS makes right axis logarithmic
set(gcf, 'Color', 'white')
set(gca, 'Color', 'white')
set(gca, 'XColor', 'black')
set(gca, 'GridColor', 'black')
title('SNDR and convergency time for each Vinc mismatch','Color','black');
legend('simulation averaged SNDR', 'simulation max SNDR','theoretical SNDR', 'adjusted SNDR', 'simulation converge times', 'theoratical converge time', 'adjusted converge time')
hold off;


