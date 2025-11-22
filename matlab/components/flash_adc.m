function [digital_out] = flash_adc(analog_in, N_bits, Vref, error_level, thresholds)
% FLASH_ADC Simulates a Flash ADC with non-linearities.
%
%   [digital_out, thresholds] = flash_adc(analog_in, N_bits, Vref, error_level)
%
%   INPUTS:
%       analog_in   - Vector of analog input samples (in Volts)
%       N_bits      - ADC resolution in bits
%       Vref        - Reference voltage (0 to Vref)
%       error_level - Controls INL/DNL severity (0 = ideal, ~0.05 = small errors, >0.2 = severe)
%
%   OUTPUTS:
%       digital_out - Quantized digital output (integer codes)
%       thresholds  - Actual ADC threshold voltages (including errors)
%
%   EXAMPLE:
%       t = 0:1e-3:1;
%       x = 0.5 + 0.4*sin(2*pi*5*t);
%       [d, thr] = flash_adc(x, 8, 1, 0.05);
%       plot(t, x, t, d/255);
%
%   (c) 2025 GPT-5 ADC Modeling Toolkit

    % --- Input checks ---
    if nargin < 4
        error('Usage: flash_adc(analog_in, N_bits, Vref, error_level)');
    end
    analog_in = analog_in(:); % force column vector

    % --- Ideal ADC step size ---
    L = 2^N_bits;
%     ideal_thresholds = linspace(0, Vref, L+1); % L+1 edges (including endpoints)
%     ideal_thresholds = ideal_thresholds(2:end-1); % remove 0 and Vref (for comparators)

%     thresholds = ideal_thresholds;
    % --- Introduce DNL/INL errors ---
    % DNL errors: small random step variations
%     dnl = error_level * (randn(size(ideal_thresholds)) * (Vref / L));
%     thresholds = ideal_thresholds + cumsum(dnl);
% 
%     % Keep thresholds monotonic (simulate realistic ADC behavior)
%     thresholds = sort(thresholds);
%     thresholds = max(min(thresholds, Vref), 0);

    % --- Quantize input ---
    digital_out = zeros(size(analog_in));
%     disp(length(thresholds))

    for k = 1:length(analog_in)
        % Count how many thresholds the input exceeds
        % thresholds must have length = L  (i.e., 2^N_bits)
        % digital_out will then be 0..L
        
        digital_out(k) = sum(analog_in(k) > thresholds);   % now outputs 0..L
    end

    % Clip to valid range
%     digital_out(digital_out < 0) = 0;
%     digital_out(digital_out > L-1) = L-1;
end