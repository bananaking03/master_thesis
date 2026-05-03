function [digi_out, SNDRs, thresholds] = flash_adc_dither_sim_simple(input, cal_len, cal_cycles, cal_constant, cal_cutoff, init_thresholds, Vhigh, Vlow,  Vinc, N_bits, non_lin_f,N, input_test)

L = 2^N_bits;
thresholds = [init_thresholds; Vhigh];   % append final overflow threshold
cal_len = round(cal_len);
nonlin_fun = @(x) polyval(fliplr(non_lin_f),x);
SNDRs = zeros(cal_cycles,1);

error_matrix = 0.5 * eye(L) + ...
    -0.25 * diag(ones(L-1,1), 1) + ...
    -0.25 * diag(ones(L-1,1), -1);
% error_matrix(end,end) = 0;
% error_matrix(end,end-1) = 0;
% error_matrix(end-1,end) = 0;
% error_matrix(end-1,end-1) = 0;
% error_matrix(end-1,end-2) = 0;
error_matrix_inv = inv(error_matrix);

for i=1:cal_cycles
    % thresholds = LSB_value * round(thresholds/LSB_value);
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
    % analog_in_amp = rescale(analog_in_amp,min(analog_in_inc),max(analog_in_inc));
    
    % adc
    adc_out = flash_adc(analog_in_amp,N_bits,Vhigh,Vlow,thresholds(1:end));

%     disp(max(adc_out))
    
    % create digital output
%     digi_out((i-1)*cal_len+1:i*cal_len) = adc_out - D;
    LSB = (Vhigh - Vlow)/L;
    digi_out = adc_out - D * (Vinc / LSB);
    
    % seperate incremented values and non-incremented values
    D_plus = digi_out(D == 1);
    D_min = digi_out(D == 0);
    
    % Compute histogram counts for each integer in range, one extra bin at
    % the end which gets discarded due to overflow
    edges = -0.5:1:(L + 0.5);   % L = 2^N_bits, so for 5-bit => -0.5:1:32.5
    H_plus = histcounts(D_plus, edges);
    H_min  = histcounts(D_min , edges);        % Hmin has bin 33????????????????????????????????
    binCenters = edges(1:end-1) + diff(edges)/2;

    % Plot H_plus
    % figure;
    % b = bar(binCenters, H_plus);   % <-- REMOVED 'hist'
    % b.FaceColor = 'flat';          % enable per-bar colors
    % b.CData = repmat([0 0.447 0.741], numel(H_plus), 1);  % default color
    % b.CData(4,:) = [1 0 0];        % recolor bar #4 red
    % b.CData(3,:) = [1 0 0];        % recolor bar #3 red
    % xlabel('ADC output code');
    % ylabel('Count');
    % title('Histogram from H\_plus');

    % Plot H_min
    % figure;
    % b = bar(binCenters, H_min);    % <-- REMOVED 'hist'
    % b.FaceColor = 'flat';
    % b.CData = repmat([0 0.447 0.741], numel(H_min), 1);   % default color
    % b.CData(5,:) = [1 0 0];        % recolor bar #4 red
    % b.CData(6,:) = [1 0 0];        % recolor bar #5 red
    % xlabel('ADC output code');
    % ylabel('Count');
    % title('Histogram from H\_min');

%     
%     % normalize, then ignore overflow bin when computing diffs:
%     H_plus = H_plus / sum(H_plus);
%     H_min  = H_min;
    m = round(Vinc / ((Vhigh - Vlow)/L));
    
    H_delta = (H_plus(1:L) - H_min(1:L)); %./ (H_plus(1:L) + H_min(1:L));   % use bins 0..31, discard overflow (bin 32)
    H_delta = H_delta(:);

    H_delta(end-1:end) = [0 0];

   % Plot H_delta
%     figure;
%     b = bar(binCenters(1:end-1), H_delta);   % <-- REMOVED 'hist'
%     b.FaceColor = 'flat';          % enable per-bar colors
%     b.CData = repmat([0 0.447 0.741], numel(H_delta), 1);  % default color
%     b.CData(5,:) = [1 0 0];        % recolor bar #4 red
% %     b.CData(2,:) = [1 0 0];        % recolor bar #3 red
%     xlabel('ADC output code');
%     ylabel('Count');
%     title('Histogram from H\_delta');

    % apply matched filter  not matched but regularized deconvolution
    % g = [-1 2 -1];
    % alpha = 0.4; % alpha = 0.2 to 0.4
    % g = [alpha, 1+ 2*alpha, alpha] ;
    % g = [0.25 0.5 0.75 1 0.75 0.5 0.25];  % example shape
    % g = [0.5 1 0.5];
    % H_delta_matched = conv(H_delta, g, 'same');

    % H_delta(end-1:end) = [0 0];

    %-----------------------------------------------------------------------
    % H_delta_matched = H_delta;
    % % 
    % % mu = 0.0000005;
    % % 
    % for k = 1:10
    %     H_delta_matched = H_delta + 0.25*[0 H_delta_matched(1:end-1)] + 0.25*[H_delta_matched(2:end) 0];
    %     % H_delta_matched = (1-mu)*H_delta_matched + mu*(H_delta + [0 0 H_delta_matched(1:end-2)] + [H_delta_matched(2:end) 0] - [0 H_delta_matched(1:end-1)]); % see distribution when using 2 LSB dither
    % end
    % H_delta_matched = H_delta_matched * (sum(H_delta)/sum(H_delta_matched));
    % H_delta_matched(end-1:end) = [0 0];
%---------------------------------------------------------------------------

    H_delta_matched = error_matrix_inv*H_delta;
    % H_delta_matched = error_matrix\H_delta;

    % Plot H_delta_matched
%     figure;
%     b = bar(binCenters(1:end-1), H_delta_matched);   % <-- REMOVED 'hist'
%     b.FaceColor = 'flat';          % enable per-bar colors
%     b.CData = repmat([0 0.447 0.741], numel(H_delta_matched), 1);  % default color
%     b.CData(5,:) = [1 0 0];        % recolor bar #4 red
% %     b.CData(2,:) = [1 0 0];        % recolor bar #3 red
%     xlabel('ADC output code');
%     ylabel('Count');
%     title('Histogram from H\_delta_matched');

   % Update thresholds to reduce nonlinearity
   % thresholds(1:end-1) = thresholds(1:end-1) + ((abs(H_delta(1:end-1)) > cal_cutoff) .* cal_constant.*H_delta(1:end-1))./cal_len;
   thresholds(1:end-1) = thresholds(1:end-1) + ((abs(H_delta_matched(1:end-1)) > cal_cutoff) .* cal_constant.*H_delta_matched(1:end-1))./(cal_len);

   % H_grad = [H_delta(1); diff(H_delta(:))];
   % H_grad = [diff(H_delta(:)); H_delta(end)];

   % m = round(Vinc / ((Vhigh - Vlow)/L));
   % 
   %  H_deconv = H_delta;
   % 
   %  for k = 2:m
   %      H_deconv(1:end-k) = H_deconv(1:end-k) - H_delta(k+1:end);
   %  end

    % H_deconv = H_delta - [zeros(m,1); H_delta(1:end-m)];

    % H_grad = [H_deconv(1); diff(H_deconv)];
    % H_delta_avg = movmean(H_delta,m); 

    % Plot H_delta_avg
%     figure;
%     b = bar(binCenters(1:end-1), H_delta_avg);   % <-- REMOVED 'hist'
%     b.FaceColor = 'flat';          % enable per-bar colors
%     b.CData = repmat([0 0.447 0.741], numel(H_delta_avg), 1);  % default color
%     b.CData(5,:) = [1 0 0];        % recolor bar #4 red
% %     b.CData(2,:) = [1 0 0];        % recolor bar #3 red
%     xlabel('ADC output code');
%     ylabel('Count');
%     title('Histogram from H\_delta');
%-------------------------------------------------------------------------
   % valid_range = 2+m-1:L-m;  % ignore edges
   % 
   %  % H_grad = [H_delta(1); diff(H_delta(:))];
   %  % H_grad = [H_deconv(1); diff(H_deconv)];
   % 
   %  update = zeros(size(thresholds));
   %  update(valid_range) = H_delta_avg(valid_range);
   % 
   %  thresholds(1:end-1) = thresholds(1:end-1) + ...
   %      ((abs(H_delta_avg(1:end-1)) > cal_cutoff) .* cal_constant .* update(1:end-1)) ./ cal_len;
% ----------------------------------------------------------------------
   % --- Enforce monotonic thresholds (prevent crossing) ---
    % Minimum spacing between thresholds
    % min_step = 1e-10;  % you can adjust this value as needed
    % 
    % for k = 2:length(thresholds)-1
    %     if thresholds(k) <= thresholds(k-1)
    %         thresholds(k) = thresholds(k-1) + min_step;
    %     end
    %     if thresholds(k) >= thresholds(k+1)
    %         thresholds(k) = thresholds(k+1) - min_step;
    %     end
    % end
    % Ensure thresholds stay within valid range
    thresholds = max(thresholds, Vlow);
    thresholds(end) = Vhigh;  % maintain top reference

    digi_out_test = flash_adc(input_test,N_bits,Vhigh,Vlow,thresholds(1:end));
    
    if cal_len >= 10000
        SNDRs(i) = calculate_SNDR(digi_out_test(end - 10000+1:end), input_test(end - 10000+1:end),N);
    else
        SNDRs(i) = calculate_SNDR(digi_out_test(end - cal_len+1:end), input_test(end - cal_len+1:end),N);
    end
    
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

