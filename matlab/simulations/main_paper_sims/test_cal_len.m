clear; close all; clc;

Vref = 1;
N_bits = 5;
cal_cycles = 1000;
fs = 10000;
non_lin_parameters = [0 1 0];

save_video = true;

L = 2^N_bits;
Vinc = Vref / L; % 1 LSB

% cal_lens = [10 100 1000 10000 100000 1000000];
cal_lens = [10000000];

% Initialize thresholds
init_thresholds = linspace(0, Vref, L+1)'; % L+1 edges
init_thresholds(4) = init_thresholds(4) + Vinc*0.9;
init_thresholds = init_thresholds(2:end-1); % remove 0 and Vref

% Preallocate storage
num_cases = length(cal_lens);
thresholds = zeros(L-1, (cal_cycles+1)*num_cases);
H_delta_diff_history = zeros(L-1, (cal_cycles)*num_cases);

% Create input vector (simulate same analog input for all)
cal_len_max = max(cal_lens);
% t = (0:1/fs:cal_len_max/fs*cal_cycles)';
% analog_in = 0.5 + 0.5*sin(2*pi*20*t);
t = (0:1/fs:(cal_len_max-1)/fs)';
analog_in = 0.5 + 0.5*sin(2*pi*20*t);
% analog_in = mod(t,1);

% Run calibration for each calibration length
for i = 1:num_cases
    cal_len = cal_lens(i);
    cal_constant = 0.1/(cal_len);
    
    % Call your calibration function
    [digi_out, thr, H_delta_diffs] = flash_adc_dither_sim(analog_in, cal_len, cal_cycles, ...
        cal_constant, init_thresholds, Vref, Vinc, N_bits, non_lin_parameters);
    
    if (i==1)
        a = digi_out;
    end

    thresholds(:, (i-1)*(cal_cycles+1)+1:i*(cal_cycles+1)) = thr; % store final thresholds
    H_delta_diff_history(:, (i-1)*(cal_cycles+1)+1:i*(cal_cycles+1)) = H_delta_diffs;
end

figure;
plot(analog_in);


%% --- Plot results ---

% Plot threshold positions vs comparator index
figure;
plot(1:L-1, thresholds, 'LineWidth', 1.5);
xlabel('Comparator index');
ylabel('Threshold Voltage (V)');
title('Calibrated Thresholds for Different Calibration Lengths');
legend(arrayfun(@(x) sprintf('Cal Len = %d', x), cal_lens, 'UniformOutput', false), ...
    'Location', 'best');
grid on;

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
    title(sprintf('ΔThreshold Evolution (Cal Len = %d)', cal_lens(i)));
    grid on;
end
sgtitle('Evolution of Adjacent Threshold Differences');

%% === Plot 2: Convergence metric (std of ΔThresholds) ===
figure('Color','w');
hold on;
for i = 1:num_cases
    block = thresholds(:, (i-1)*(cal_cycles+1)+1 : i*(cal_cycles+1));
    deltaT = diff(block, 1, 1); % differences between thresholds
    spread = std(deltaT); % standard deviation of ΔT per iteration
    plot(0:cal_cycles, spread, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Cal Len = %d', cal_lens(i)));
end
xlabel('Calibration Iteration');
ylabel('Std of Adjacent ΔThresholds (V)');
title('Convergence of Threshold Spacing');
legend('Location', 'northeast');
grid on;

% Optional: compare histograms of final digital outputs
%% --- Plot 1: Threshold evolution for each calibration length ---
figure;
for i = 1:num_cases
    subplot(2, ceil(num_cases/2), i);
    block = thresholds(:, (i-1)*(cal_cycles+1)+1 : i*(cal_cycles+1));
    plot(block', 'LineWidth', 1);
    title(sprintf('Cal Len = %d', cal_lens(i)));
    xlabel('Calibration Iteration');
    ylabel('Threshold Voltage (V)');
    grid on;
end
sgtitle('Threshold Evolution Across Calibration Cycles');

%% --- Plot 2: Final thresholds comparison ---
figure;
hold on;
for i = 1:num_cases
    % Get final column (last calibration state)
    final_thr = thresholds(:, i*(cal_cycles+1));
    plot(1:L-1, final_thr, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Cal Len = %d', cal_lens(i)));
end
xlabel('Comparator Index');
ylabel('Threshold Voltage (V)');
title('Final Thresholds for Different Calibration Lengths');
legend('Location', 'best');
grid on;

%% --- Plot 3 (optional): Convergence metric over calibration cycles ---
figure;
hold on;
for i = 1:num_cases
    block = thresholds(:, (i-1)*(cal_cycles+1)+1 : i*(cal_cycles+1));
    % Measure standard deviation of threshold spacing (as a rough convergence indicator)
    thr_spread = std(diff(block));
    plot(0:cal_cycles, thr_spread, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Cal Len = %d', cal_lens(i)));
end
xlabel('Calibration Iteration');
ylabel('Std of Threshold Spacing (V)');
title('Convergence Behavior of Threshold Calibration');
legend('Location', 'northeast');
grid on;

%% --- Animation Setup ---
case_to_show = 1; % choose which calibration length to animate (index in cal_lens)
block = thresholds(:, (case_to_show-1)*(cal_cycles+1)+1 : case_to_show*(cal_cycles+1));
block_H = H_delta_diff_history(:,(case_to_show-1)*(cal_cycles)+1:case_to_show*(cal_cycles));

figure('Color','w');
tiledlayout(1,2);

% Left plot for thresholds
ax1 = nexttile;
x = 1:L-1;
h_thr = plot(x, block(:,1), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
hold(ax1,'on');
% f = @(x) -0.2*(x/(L-1))^2+x/(L-1);
% fplot(f);
ylim([0 Vref]); xlim([1 L-1]);
xlabel(ax1,'Comparator Index'); ylabel(ax1,'Threshold Voltage (V)');
title(ax1, sprintf('Threshold Evolution (Cal Len = %d)', cal_lens(case_to_show)));
grid(ax1,'on');

% Left plot for thresholds
ax2 = nexttile;
bar(x,block_H(:,1))

x = 1:L-1;
h_hist = plot(x, block_H(:,1), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5); 
% xlim([1 L-1]);
xlabel(ax2,'Comparator Index'); ylabel(ax2,'count');
title(ax2, sprintf('Histogram difference Evolution (Cal Len = %d)', cal_lens(case_to_show)));
grid(ax2,'on');

% Uncomment to record as a video
if save_video
    vid = VideoWriter('threshold_evolution.mp4', 'MPEG-4');
    vid.FrameRate = 20;
    open(vid);
end

% Animate through calibration cycles
for k = 1:size(block(:,2:end),2)
    % update thresholds
    set(h_thr, 'YData', block(:,k));
    title(sprintf('Threshold Evolution (Cal Len = %d) | Cycle %d/%d', ...
        cal_lens(case_to_show), k, cal_cycles));

    % update histogram
    set(h_hist, 'YData', block_H(:,k));
    title(sprintf('Histogram Evolution (Cal Len = %d) | Cycle %d/%d', ...
        cal_lens(case_to_show), k, cal_cycles));

    drawnow;
    if save_video
        writeVideo(vid, getframe(gcf));
    end
end

if save_video
    close(vid);
    disp('✅ Saved animation as threshold_evolution.mp4');
end