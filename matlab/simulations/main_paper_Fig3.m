Vref = 1;
N_bits = 5;
non_lin_adc = 0; % level of nonlinearity of adc;
cal_cycles = 100;
cal_len = 100000;
fs =  10000;

cal_constant = cal_len/(10000000*cal_len);
L = 2^N_bits;
Vinc = Vref/((2^N_bits)); % 1 LSB 

thresholds = linspace(0, Vref, L+1)'; % L+1 edges (including endpoints)
thresholds = thresholds(2:end-1); % remove 0 and Vref (for comparators)
threshold_history = zeros(length(thresholds), cal_cycles);
disp(size(thresholds));

for i=1:cal_cycles
    t = (0:1/fs:cal_len/fs)';
    analog_in = 0.5 + 0.4*sin(2*pi*10*t);
    D = randi([0 1], length(analog_in),1);
    
    % add dither
    analog_in_inc = analog_in + D*Vinc;

    % apply non-linear amplification (no amp needed in matlab but simulates
    % non-linearity)
    analog_in_amp = -0.1*analog_in_inc.^2 + analog_in_inc;
    
    % adc
    adc_out = flash_adc(analog_in_amp,N_bits,Vref,non_lin_adc,thresholds);
    
    % create digital output
    digi_out = adc_out - D;
    
    % seperate incremented values and non-incremented values
    D_plus = digi_out(D == 1);
    D_min = digi_out(D == 0);
    
    % Compute histogram counts for each integer in range
    edges = 0:2^N_bits;
    H_plus = histcounts(D_plus, edges);
    H_min = histcounts(D_min, edges);

   % Compute the difference between the amount of occurances between
   % adc outputs of D_plus and D_min
   H_delta = H_plus - H_min; % 0 if linear (and enough cal_len)

   % Compute the difference between the difference in occurances in
   % subsequent adc outputs.
   H_delta_diff = zeros(2^N_bits-1,1);
   for j=1:2^N_bits-1
        H_delta_diff(j) = H_delta(j+1) - H_delta(j);
   end

   % Update thresholds to reduce nonlinearity
   thresholds = thresholds - cal_constant*H_delta_diff;

   % Save thresholds for visualization
   threshold_history(:, i) = thresholds;
end

figure;
plot(threshold_history', 'LineWidth', 1);
xlabel('Calibration iteration');
ylabel('Threshold voltage (V)');
title('Evolution of ADC thresholds during calibration');
grid on;

figure;
imagesc(1:cal_cycles, 1:length(thresholds), threshold_history);
colorbar;
xlabel('Calibration iteration');
ylabel('Threshold index');
title('Threshold evolution (heatmap view)');

% Compute total movement from initial thresholds
movement = vecnorm(diff(threshold_history, 1, 2)); % Euclidean norm per iteration

figure;
plot(movement, 'LineWidth', 1.5);
xlabel('Calibration iteration');
ylabel('Threshold movement (norm)');
title('Convergence of threshold calibration');
grid on;

% plot histograms
figure; 
subplot(2,1,1);
histogram(D_plus, 'BinMethod', 'integers');
xlabel('Value');
ylabel('Count');
title('Conversion with increment');
subplot(2,1,2);
histogram(D_min, 'BinMethod', 'integers');
xlabel('Value');
ylabel('Count');
title('Conversion without increment');
