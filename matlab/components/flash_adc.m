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
    thresholds = thresholds(:).'; % row vector for comparison

    % --- Quantize input ---
    % Use chunking for large input arrays to avoid oversized temporary matrices.
    maxRowsPerChunk = 1e6;
    numSamples = numel(analog_in);
    digital_out = zeros(numSamples, 1, 'like', analog_in);

    if numSamples == 0
        return;
    end

    if numSamples * numel(thresholds) <= 1e8
        digital_out = sum(analog_in > thresholds, 2);
    else
        for startIdx = 1:maxRowsPerChunk:numSamples
            endIdx = min(startIdx + maxRowsPerChunk - 1, numSamples);
            digital_out(startIdx:endIdx) = sum(analog_in(startIdx:endIdx) > thresholds, 2);
        end
    end
    % digital_out = digital_out + (analog_in > Vhigh);

    % Clip to valid range
%     digital_out(digital_out < 0) = 0;
%     digital_out(digital_out > L-1) = L-1;
end