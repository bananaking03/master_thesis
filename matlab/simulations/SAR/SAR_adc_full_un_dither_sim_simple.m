function [digi_out, SNDRs, thresholds] = SAR_adc_full_un_dither_sim_simple(input, cal_len, cal_cycles, ...
    cal_constant, cal_cutoff, init_thresholds, Vhigh, Vlow,  Vinc, N_bits, non_lin_f,N, ...
    input_test, pos_thresholds, L_algo, L_tot, N_bits_extra)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% cal_len = round(cal_len);
nonlin_fun = @(x) polyval(fliplr(non_lin_f),x);
SNDRs = zeros(cal_cycles,1);
thresholds = init_thresholds;
edges = -0.5:1:(L_algo + 0.5);   % L = 2^N_bits, so for 5-bit => -0.5:1:32.5
binCenters = edges(1:end-1) + diff(edges)/2;

L_algo = 2^N_bits;

DAC_select = 1:2^N_bits_extra:L_tot;

% DAC_select = round(DAC_select + randn(size(DAC_select))*1500);
% DAC_select = sort(DAC_select);
% DAC_select = min(max(DAC_select,1),L_tot);

% DAC_select = 1:(2^N_bits_extra)/2:L_tot/2;
% DAC_select(50) = round(DAC_select(50)+ 0.9*2^(N_bits_extra));


% start calibration
for i=1:cal_cycles

    analog_in = input;
    D = randi([0 1], cal_len,1);

    % add dither
    analog_in_inc = analog_in + D*Vinc;

    analog_in_amp = nonlin_fun(analog_in_inc)';

    % adc
    % adc_out = flash_adc(analog_in_amp,N_bits,Vhigh,Vlow,thresholds(1:4:end));
    adc_out = flash_adc(analog_in_amp,N_bits,Vhigh,Vlow,thresholds);

    LSB_algo = (Vhigh - Vlow)/L_algo;
    % digi_out = adc_out - D * (Vinc / LSB_algo);
    digi_out = adc_out - D;

    % seperate incremented values and non-incremented values
    D_plus = digi_out(D == 1);
    D_min = digi_out(D == 0);

    % Compute histogram counts for each integer in range, one extra bin at
    % the end which gets discarded due to overflow
    H_plus = histcounts(D_plus, edges);
    H_min  = histcounts(D_min , edges);        % Hmin has bin 33????????????????????????????????

    H_delta = (H_plus(1:L_algo) - H_min(1:L_algo)); %./ (H_plus(1:L) + H_min(1:L));   % use bins 0..31, discard overflow (bin 32)
    H_delta = H_delta(:);

    H_delta(end-2:end) = [0 0 0];

    %% update thresholds
    threshold_update = ((abs(H_delta(1:end-1)) > cal_cutoff) .* (L_tot/(Vhigh-Vlow)) .* cal_constant.*H_delta(1:end-1))./(cal_len);

    % DAC_select(1:4:end-4) = DAC_select(1:4:end-4) + threshold_update.';
    % DAC_select(end-3:end) = DAC_select(end-3:end) + threshold_update;

    DAC_select(1:end-1) = DAC_select(1:end-1) + threshold_update.';

    % Interpolate missing values
    % idx = (1:4:length(DAC_select));
    % DAC_select(1:end-3) = interp1(idx, DAC_select(1:4:end), 1:(length(DAC_select)-3), 'linear');
    % 
    % % clip it
    DAC_select = min(max(DAC_select,1),L_tot);

    % disp(DAC_select);
    thresholds = pos_thresholds(round(DAC_select));

    thresholds = max(thresholds, Vlow);
    % disp(thresholds(end));
    % thresholds(end) = Vhigh;  % maintain top reference

    digi_out_test = flash_adc(input_test,N_bits,Vhigh,Vlow,thresholds(1:end));
    
    if cal_len >= 10000
        SNDRs(i) = calculate_SNDR(digi_out_test(end - 10000+1:end), input_test(end - 10000+1:end),N);
    else
        SNDRs(i) = calculate_SNDR(digi_out_test(end - cal_len+1:end), input_test(end - cal_len+1:end),N);
    end
end