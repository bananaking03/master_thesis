clear;
clc;
close all;

% Read exported SimVision data H_matched
data = readtable('histogram_matched_outputs_callen2000_schemcomp.csv', ...
    'VariableNamingRule', 'preserve');

% Display column names
% disp(data.Properties.VariableNames);

%% Convert into decimal

H_matched = zeros(height(data), 1, 'int64');

for i = 2:width(data)

    % Extract hexadecimal string
    hex_str = string(data{2,i});
    hex_str = regexprep(hex_str, '[^0-9A-Fa-f]', '');

    % Pad to 16 hexadecimal digits
    hex_str = pad(hex_str, 16, 'left', '0');

    % Split into two 32-bit parts
    high = uint64(hex2dec(extractBetween(hex_str, 1, 8)));
    low  = uint64(hex2dec(extractBetween(hex_str, 9, 16)));

    % Combine into 64-bit integer
    value = bitor(bitshift(high, 32), low);

    % Interpret as signed 64-bit two's complement
    H_matched(i-1) = typecast(value, 'int64');

end

%% Read exported SimVision data H_plus
data = readtable('histogram_plus_outputs_callen2000_schemcomp.csv', ...
    'VariableNamingRule', 'preserve');

% Keep only one row per timestamp
[~, idx] = unique(data.Var1, 'stable');
data = data(idx,:);

% Remove parentheses
data.Var2 = erase(string(data.Var2), '(');
data.Var257 = erase(string(data.Var257), ')');

timestamps = data.Var1;

H_plus = [ ...
    str2double(data.Var2), ...
    data{:,3:256}, ...
    str2double(data.Var257) ...
];

%% Read exported SimVision data H_min
data = readtable('histogram_min_outputs_callen2000_schemcomp.csv', ...
    'VariableNamingRule', 'preserve');

% Keep only one row per timestamp
[~, idx] = unique(data.Var1, 'stable');
data = data(idx,:);

% Remove parentheses
data.Var2 = erase(string(data.Var2), '(');
data.Var257 = erase(string(data.Var257), ')');

timestamps = data.Var1;

H_min = [ ...
    str2double(data.Var2), ...
    data{:,3:256}, ...
    str2double(data.Var257) ...
];

%%
H_diff = H_plus.' - H_min.';

L_algo = 256;
lambda_reg = 0.1;

error_matrix = 0.5 * eye(L_algo) + ...
    -0.25 * diag(ones(L_algo-1,1), 1) + ...
    -0.25 * diag(ones(L_algo-1,1), -1);

H_delta_matched_full = (error_matrix + lambda_reg*eye(L_algo))\H_diff;

%%

%% Compare H_matched (SimVision) vs H_delta_matched_full (MATLAB)

t_idx = size(H_delta_matched_full, 2);   % timestamp column to compare (last by default)
scale = 1;                               % set to 2^frac_bits if H_matched is fixed-point

H_matched_vec = double(H_matched(:)) / scale;
H_delta_vec   = H_delta_matched_full(:, t_idx);

% Sanity check on lengths
assert(numel(H_matched_vec) == numel(H_delta_vec), ...
    'Length mismatch: H_matched has %d elements, H_delta has %d.', ...
    numel(H_matched_vec), numel(H_delta_vec));

bins = 1:L_algo;

figure('Name', 'Matched filter comparison', 'Position', [100 100 1200 400]);

ax1 = subplot(1,2,1);
bar(bins, H_matched_vec);
title('H\_matched (SimVision)');
xlabel('Bin'); ylabel('Value');
xlim([0 L_algo+1]); grid on;

ax2 = subplot(1,2,2);
bar(bins, H_delta_vec);
title(sprintf('H\\_delta\\_matched\\_full (MATLAB, t = %g)', timestamps(t_idx)));
xlabel('Bin'); ylabel('Value');
xlim([0 L_algo+1]); grid on;

linkaxes([ax1 ax2], 'x');   % zooming one plot zooms the other

% Quick numeric check
fprintf('Max abs difference: %g\n', max(abs(H_matched_vec - H_delta_vec)));