function sndr_dB = calculate_SNDR(x_adc,x_in,N)
% Ensure column vectors
    x_adc = x_adc(:);
    x_in  = x_in(:);
    % N = length(x_adc);
    % 
    % % Remove DC
    % x_adc = x_adc - mean(x_adc);
    % x_in  = x_in  - mean(x_in);
    % 
    % % Window (optional but recommended)
    % w = hann(N);
    % x_adc_w = x_adc .* w;
    % x_in_w  = x_in  .* w;

    % FFT
    X_adc = fft(x_adc,N) / N;
    X_in  = fft(x_in,N) / N;

    % Single-sided spectrum
    X_adc = X_adc(1:(N-1)/2);
    X_in  = X_in(1:(N-1)/2);

    % Frequency bin of fundamental (from input reference)
    [~, fund_bin] = max(abs(X_in));

    % Power spectrum
    P = abs(X_adc).^2;

    % Fundamental power (±1 bin for leakage safety)
    % fund_bins = fund_bin-1 : fund_bin+1;
    P_signal = P(fund_bin);

    % Total power excluding DC
    P_total = sum(P(2:end-1));

    % Noise + distortion power
    P_noise_dist = P_total - P_signal;

    % SNDR
    sndr_dB = 10*log10(P_signal / P_noise_dist);
end