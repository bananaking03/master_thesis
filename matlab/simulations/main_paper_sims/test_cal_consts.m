clear; close all; clc;
% change to differantial system!!!!! Check distortie
Vhigh = 1;
Vlow = -1;
N_bits =10; % probeer meer bits
% cal_cycles = 5000;
cal_cycles = 5000;
N = (2048*2^-3) - 1; % fft size
fs = 48000;  % coherent sampling
f0 = (13/N)*fs;
% f0 = 990;
non_lin_parameters = [0 1 0];
cal_len = 40000;  % need increase for more bits
% N = 2048*2^4; % fft size
analyze_spesific = 1;
cal_cutoff = 0;

save_video = true;

L = 2^N_bits;
Vinc = (Vhigh - Vlow) / L; % 1 LSB

% cal_consts = [8*10^(-5) 9*10^(-5) 10^(-4) 2*10^(-4) 3*10^(-4) 4*10^(-4)];
% cal_consts = [10^(-3) 10^(-4) 10^(-5) 10^(-6) 10^(-7) 10^(-8)];
cal_consts = logspace(-8, -4, 10)';
% cal_consts = [10^-6];

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
H_delta_diff_history = zeros(L-1, (cal_cycles));
err = zeros(num_cases,1);
SNDRs = zeros(num_cases,1);

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
    [digi_out, thr, H_delta_diffs] = flash_adc_dither_sim(analog_in, cal_len, cal_cycles, ...
        cal_constant, cal_cutoff, init_thresholds, Vhigh, Vlow, Vinc, N_bits, non_lin_parameters);
    
    % err(i) = sum(ideal_digi_out ~= digi_out)/length(digi_out);
    SNDRs(i) = calculate_SNDR(digi_out(:,cal_cycles),analog_in,N);
    if (i==analyze_spesific)
        digi_out_spes = digi_out(:,cal_cycles);
        H_delta_diff_history = H_delta_diffs;
    end

    thresholds(:, (i-1)*(cal_cycles+1)+1:i*(cal_cycles+1)) = thr; % store final thresholds
    % H_delta_diff_history(:, (i-1)*(cal_cycles+1)+1:i*(cal_cycles+1)) = H_delta_diffs;
end

disp(err);
figure;
semilogx(cal_consts, SNDRs, '*');
hold on;
yline(calculate_SNDR(initial_digi_out,analog_in,N), '--r', 'initial SNDR');
yline(calculate_SNDR(ideal_digi_out,analog_in,N), '--r', 'ideal SNDR');
xlabel('cal constant');
ylabel('SNDR achieved in last iteration (dB)');
title('SNDR for each cal constant');
grid on;
hold off;

figure;
plot(analog_in);

figure;
plot(digi_out_spes);


%% --- Plot results ---

% Plot threshold positions vs comparator index
% figure;
% plot(1:L-1, thresholds, 'LineWidth', 1.5);
% xlabel('Comparator index');
% ylabel('Threshold Voltage (V)');
% title('Calibrated Thresholds for Different Calibration Lengths');
% legend(arrayfun(@(x) sprintf('Cal Len = %d', x), cal_consts, 'UniformOutput', false), ...
%     'Location', 'best');
% grid on;

%% === Plot 1: Differences between adjacent thresholds over time ===
figure('Color','w');
for i = 1:num_cases
    subplot(2, ceil(num_cases/2), i);
    block = thresholds(:, (i-1)*(cal_cycles+1)+1 : i*(cal_cycles+1));
    deltaT = diff(block, 1, 1); % differences between adjacent thresholds
    
    % Plot each ΔT line over cycles
    plot(deltaT', 'LineWidth', 1);
    xlabel('Calibration Iteration');
    ylabel('\Delta Threshold (V)');
    title(sprintf('ΔThreshold Evolution (Cal Len = %d)', cal_consts(i)));
    grid on;
end
sgtitle('Evolution of Adjacent Threshold Differences');

%% === Plot 2: Convergence metric (std of ΔThresholds) ===
tau = zeros(num_cases,1);
figure('Color','w');
hold on;
for i = 1:num_cases
    block = thresholds(:, (i-1)*(cal_cycles+1)+1 : i*(cal_cycles+1));
    threshold_error = block - ideal_thresholds; % differences between thresholds  !!!!!!!!!!!!!!!!!!!change to difference with ideal thresholds
    spread = std(threshold_error); % standard deviation of ΔT per iteration
    plot(0:cal_cycles, spread, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Cal const = %d', cal_consts(i)));

    % calculate convergence speed
    n = 1:cal_cycles;

    spread_inf = mean(spread(end-100:end));   % steady-state estimate
    y = abs(spread(2:end) - spread_inf);              % remove offset

    % remove zeros to avoid log issues
    valid = y > 0;
    n_fit = n(valid);
    y_fit = y(valid);
    
    p = polyfit(n_fit, log(y_fit), 1);
    
    tau(i) = -1/p(1);   % convergence time constant
end
xlabel('Calibration Iteration');
ylabel('Std of Adjacent ΔThresholds (V)');
title('Convergence of Threshold Spacing');
legend('Location', 'northeast');
grid on;

figure;
semilogx(cal_consts,tau*cal_len, '*');
hold on;
xlabel('Calibration constant');
ylabel('cycles to convergence');
title('Convergence speed (cal_len = 40000');
hold off;

% Optional: compare histograms of final digital outputs
%% --- Plot 1: Threshold evolution for each calibration length ---
figure;
for i = 1:num_cases
    subplot(2, ceil(num_cases/2), i);
    block = thresholds(:, (i-1)*(cal_cycles+1)+1 : i*(cal_cycles+1));
    plot(block', 'LineWidth', 1);
    title(sprintf('Cal Len = %d', cal_consts(i)));
    xlabel('Calibration Iteration');
    ylabel('Threshold Voltage (V)');
    grid on;
end
sgtitle('Threshold Evolution Across Calibration Cycles');

%% --- Plot frequency response
% FFT
Nfft = 1*N;
win    = hann(N)';            % Hann window
Wnorm  = sum(win) / N;        % Normalize amplitude for window loss
dbFloor = 0;                % Baseline level for all spectra
fft_db = @(x) max( to_dB( fft(x, Nfft)), dbFloor );
fft_in  = fft_db(analog_in(1:end)');
fft_out = fft_db((digi_out_spes(1:end)'-L/2)/(L/2));
fft_out_ideal = fft_db((ideal_digi_out(1:end)'-L/2)/(L/2));
f = (0:Nfft-1)*(fs/Nfft);

figure;
bar(f(1:Nfft/2), fft_in(1:Nfft/2),0.05,'FaceColor','b','EdgeColor','b', 'LineWidth', 0.8,'BaseValue',dbFloor); hold on;
bar(f(1:Nfft/2), fft_out(1:Nfft/2),0.05,'FaceColor','r','EdgeColor','r', 'LineWidth', 0.1,'FaceAlpha',0.7,'BaseValue',dbFloor);
title('Frequency Domain – Harmonics at 2f, 3f, …');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
% xlim([0 5000]); 
ylim([dbFloor max([fft_out fft_in 60])]);
ylabel('Magnitude (dB)');
legend('Input','Output');

figure;
bar(f(1:Nfft/2), fft_in(1:Nfft/2),0.05,'FaceColor','b','EdgeColor','b', 'LineWidth', 0.8,'BaseValue',dbFloor); hold on;
bar(f(1:Nfft/2), fft_out_ideal(1:Nfft/2),0.05,'FaceColor','g','EdgeColor','g', 'LineWidth', 0.3,'FaceAlpha',0.7,'BaseValue',dbFloor);
title('Frequency Domain – Harmonics at 2f, 3f, …');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
% xlim([0 5000]); 
ylim([dbFloor max([fft_out fft_in 60])]);
ylabel('Magnitude (dB)');
legend('Input','Output ideal');

%% --- Animation Setup ---  threshold evolution
case_to_show = analyze_spesific; % choose which calibration length to animate (index in cal_consts)
block = thresholds(:, (case_to_show-1)*(cal_cycles+1)+1 : case_to_show*(cal_cycles+1));
% block_H = H_delta_diff_history(:,(case_to_show-1)*(cal_cycles)+1:case_to_show*(cal_cycles));
block_H = H_delta_diff_history;

figure('Color','w');
tiledlayout(1,2);

% Left plot for thresholds
ax1 = nexttile;
x = 1:L-1;
h_thr = plot(x, block(:,1), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
hold(ax1,'on');
% f = @(x) -0.2*(x/(L-1))^2+x/(L-1);
% fplot(f);
ylim([Vlow Vhigh]); xlim([1 L-1]);
xlabel(ax1,'Comparator Index'); ylabel(ax1,'Threshold Voltage (V)');
title(ax1, sprintf('Threshold Evolution (Cal Len = %d)', cal_consts(case_to_show)));
grid(ax1,'on');

% Left plot for thresholds
ax2 = nexttile;
bar(x,block_H(:,1))

x = 1:L-1;
h_hist = plot(x, block_H(:,1), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5); 
% xlim([1 L-1]);
xlabel(ax2,'Comparator Index'); ylabel(ax2,'count');
title(ax2, sprintf('Histogram difference Evolution (Cal Len = %d)', cal_consts(case_to_show)));
grid(ax2,'on');

% Uncomment to record as a video
if save_video
    vid = VideoWriter('threshold_evolution.mp4', 'MPEG-4');
    vid.FrameRate = 30;
    open(vid);
end

% Animate through calibration cycles
for k = 1:size(block(:,2:end),2)
    % update thresholds
    set(h_thr, 'YData', block(:,k));
    title(sprintf('Threshold Evolution (Cal Len = %d) | Cycle %d/%d', ...
        cal_consts(case_to_show), k, cal_cycles));

    % update histogram
    set(h_hist, 'YData', block_H(:,k));
    title(sprintf('Histogram Evolution (Cal Len = %d) | Cycle %d/%d', ...
        cal_consts(case_to_show), k, cal_cycles));

    drawnow;
    if save_video
        writeVideo(vid, getframe(gcf));
    end
end

if save_video
    close(vid);
    disp('✅ Saved animation as threshold_evolution.mp4');
end