L = 2^4;
H_delta = zeros(2^4,1);
H_delta(3) = -1/2*100;
H_delta(4) = 1*100;
H_delta(5) = -1/2*100;

error_matrix = 0.5 * eye(L) + ...
    -0.25 * diag(ones(L-1,1), 1) + ...
    -0.25 * diag(ones(L-1,1), -1);
error_matrix_inv = inv(error_matrix);

error = error_matrix_inv*H_delta;

edges = -0.5:1:(2^4-1 + 0.5);   % L = 2^N_bits, so for 5-bit => -0.5:1:32.5
binCenters = edges(1:end-1) + diff(edges)/2;

figure;
b = bar(binCenters, H_delta);   % <-- REMOVED 'hist'
b.FaceColor = 'flat';          % enable per-bar colors
b.CData = repmat([0 0.447 0.741], numel(H_delta), 1);  % default color
b.CData(4,:) = [1 0 0];        % recolor bar #4 red
b.CData(3,:) = [1 0 0];        % recolor bar #3 red
xlabel('ADC output code');
ylabel('Count');
title('Histogram from H\_plus');

figure;
b = bar(binCenters, error);   % <-- REMOVED 'hist'
b.FaceColor = 'flat';          % enable per-bar colors
b.CData = repmat([0 0.447 0.741], numel(error), 1);  % default color
b.CData(4,:) = [1 0 0];        % recolor bar #4 red
b.CData(3,:) = [1 0 0];        % recolor bar #3 red
xlabel('ADC output code');
ylabel('Count');
title('Histogram from H\_plus');