function [digi_out, threshold_history, H_delta_diff_history] = flash_adc_dither_sim(input, cal_len, cal_cycles, cal_constant, init_thresholds, Vref, Vinc, N_bits, non_lin_f)

L = 2^N_bits;
H_delta_diff_history = zeros(L-1, cal_cycles);
thresholds = [init_thresholds; Vref];   % append final overflow threshold
threshold_history = zeros(L-1, cal_cycles+1);
threshold_history(:,1) = init_thresholds;
% digi_out = zeros(cal_len*cal_cycles,1);
digi_out = zeros(cal_len,1);
nonlin_fun = @(x) polyval(fliplr(non_lin_f),x);

for i=1:cal_cycles
%     analog_in = input((i-1)*cal_len+1:(i)*cal_len);
    analog_in = input;
    D = randi([0 1], cal_len,1);
    
    % add dither
    analog_in_inc = analog_in + D*Vinc;

    % apply non-linear amplification (no amp needed in matlab but simulates
    % non-linearity
    analog_in_amp = nonlin_fun(analog_in_inc)';
%     figure
%     plot(analog_in_amp)
%     xlim([0 10000])
    analog_in_amp = rescale(analog_in_amp,min(analog_in_inc),max(analog_in_inc));
    
    % adc
    adc_out = flash_adc(analog_in_amp,N_bits,Vref,0,thresholds(1:end));

%     disp(max(adc_out))
    
    % create digital output
%     digi_out((i-1)*cal_len+1:i*cal_len) = adc_out - D;
    digi_out = adc_out - D;
    
    % seperate incremented values and non-incremented values
    D_plus = digi_out(D == 1);
    D_min = digi_out(D == 0);
    
    % Compute histogram counts for each integer in range, one extra bin at
    % the end which gets discarded due to overflow
    edges = -0.5:1:(L + 0.5);   % L = 2^N_bits, so for 5-bit => -0.5:1:32.5
    H_plus = histcounts(D_plus, edges);
    H_min  = histcounts(D_min , edges);        % Hmin has bin 33????????????????????????????????

    binCenters = edges(1:end-1) + diff(edges)/2;

%     % Plot H_plus
%     figure;
%     b = bar(binCenters, H_plus);   % <-- REMOVED 'hist'
%     b.FaceColor = 'flat';          % enable per-bar colors
%     b.CData = repmat([0 0.447 0.741], numel(H_plus), 1);  % default color
%     b.CData(3,:) = [1 0 0];        % recolor bar #4 red
%     b.CData(2,:) = [1 0 0];        % recolor bar #3 red
%     xlabel('Value');
%     ylabel('Count');
%     title('Histogram from H\_plus');
%     
%     % Plot H_min
%     figure;
%     b = bar(binCenters, H_min);    % <-- REMOVED 'hist'
%     b.FaceColor = 'flat';
%     b.CData = repmat([0 0.447 0.741], numel(H_min), 1);   % default color
%     b.CData(3,:) = [1 0 0];        % recolor bar #4 red
%     b.CData(4,:) = [1 0 0];        % recolor bar #5 red
%     xlabel('Value');
%     ylabel('Count');
%     title('Histogram from H\_min');

%     
%     % normalize, then ignore overflow bin when computing diffs:
%     H_plus = H_plus / sum(H_plus);
%     H_min  = H_min;
    
    H_delta = H_plus(1:L) - H_min(1:L);   % use bins 0..31, discard overflow (bin 32)


   % Compute the difference between the difference in occurances in
   % subsequent adc outputs.
   H_delta_diff = diff(H_delta);

%    % Plot H_plus
%     figure;
%     b = bar(binCenters(1:end-1), H_delta);   % <-- REMOVED 'hist'
%     b.FaceColor = 'flat';          % enable per-bar colors
%     b.CData = repmat([0 0.447 0.741], numel(H_delta), 1);  % default color
%     b.CData(3,:) = [1 0 0];        % recolor bar #4 red
%     b.CData(2,:) = [1 0 0];        % recolor bar #3 red
%     xlabel('Value');
%     ylabel('Count');
%     title('Histogram from H\_delta');
%     
%     % Plot H_min
%     figure;
%     b = bar(binCenters(1:end-2), H_delta_diff);    % <-- REMOVED 'hist'
%     b.FaceColor = 'flat';
%     b.CData = repmat([0 0.447 0.741], numel(H_delta_diff), 1);   % default color
%     b.CData(3,:) = [1 0 0];        % recolor bar #4 red
%     b.CData(4,:) = [1 0 0];        % recolor bar #5 red
%     xlabel('Value');
%     ylabel('Count');
%     title('Histogram from H\_delta_diff');
%    disp(H_delta_diff);
%    H_delta_diff = zeros(2^N_bits-1,1);
%    for j=1:2^N_bits-1
%         H_delta_diff(j) = H_delta(j+1) - H_delta(j);
%    end

   % Update thresholds to reduce nonlinearity
   thresholds(1:end-1) = thresholds(1:end-1) - cal_constant*H_delta_diff(1:end)';

   % --- Enforce monotonic thresholds (prevent crossing) ---
    % Minimum spacing between thresholds
    min_step = 1e-6;  % you can adjust this value as needed
    
    for k = 2:length(thresholds)-1
        if thresholds(k) <= thresholds(k-1)
            thresholds(k) = thresholds(k-1) + min_step;
        end
        if thresholds(k) >= thresholds(k+1)
            thresholds(k) = thresholds(k+1) - min_step;
        end
    end
    % Ensure thresholds stay within valid range
    thresholds = max(thresholds, 0);
    thresholds(end) = Vref;  % maintain top reference

   % Save thresholds for visualization
   threshold_history(:, i+1) = thresholds(1:end-1);

   % Save histograms for visualization
   H_delta_diff_history(:,i+1) = H_delta(2:end)';

%    figure
%    plot(thresholds)
%    title('thresholds')
%    figure
%    plot(H_delta)
%    title('H_delta')
%    figure
%    plot(H_delta_diff)
%    title('H_delta_diff')
%    figure
%    plot(H_plus)
%    title('H_plus')
%    figure
%    plot(H_min)
%    title('H_min')
end

end

