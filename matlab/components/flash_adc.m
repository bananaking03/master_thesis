function [digital_out] = flash_adc(analog_in, N_bits, Vhigh, Vlow, thresholds)
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

    % --- Input checks ---
    if nargin < 4
        error('Usage: flash_adc(analog_in, N_bits, Vref, error_level)');
    end
    analog_in = analog_in(:); % force column vector

    % --- Ideal ADC step size ---
    % L = 2^N_bits;

    % --- Quantize input ---
    % digital_out = zeros(size(analog_in));

    % for k = 1:length(analog_in)
    %     % Count how many thresholds the input exceeds
    %     % thresholds must have length = L  (i.e., 2^N_bits)
    %     % digital_out will then be 0..L
    % 
    %     digital_out(k) = sum(analog_in(k) > thresholds);   % now outputs 0..L
    % 
    %     digital_out(k) = digital_out(k) + (analog_in(k) > Vhigh); % may need removal
    % end

    digital_out = sum(analog_in > thresholds.', 2);
    % digital_out = digital_out + (analog_in > Vhigh);

    % Clip to valid range
%     digital_out(digital_out < 0) = 0;
%     digital_out(digital_out > L-1) = L-1;
end