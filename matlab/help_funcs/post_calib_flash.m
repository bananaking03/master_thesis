function [digital_out, code, recon_levels] = post_calib_flash(analog_in, N_bits, Vhigh, Vlow, thresholds)

    analog_in = analog_in(:);

    L = 2^N_bits;

    % --- Compute code index (0 to L) ---
    code = sum(analog_in > thresholds.', 2);

    % --- Build reconstruction levels ---
    % Add end points
    full_thresholds = [Vlow; thresholds(:); Vhigh];

    % Midpoint of each quantization bin
    recon_levels = (full_thresholds(1:end-1) + full_thresholds(2:end)) / 2;

    % Clip code to valid range
    code(code < 0) = 0;
    code(code > L-1) = L-1;

    % --- Output reconstructed voltage instead of code ---
    digital_out = recon_levels(code + 1);
end