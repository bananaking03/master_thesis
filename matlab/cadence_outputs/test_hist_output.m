% clear;
% clc;
% close all;
% 
% H_plus_file = 'histogram_plus_callen10000_va.csv';
% H_min_file = 'histogram_min_callen10000_va.csv';
% H_matched_file = 'histogram_matched_callen10000_va.csv';
% 
% % Read exported SimVision data H_matched
% data = readtable(H_matched_file, ...
%     'VariableNamingRule', 'preserve');
% 
% % Display column names
% % disp(data.Properties.VariableNames);
% 
% %% Convert into decimal
% 
% H_matched = zeros(height(data), 1, 'int64');
% 
% for i = 2:width(data)
% 
%     % Extract hexadecimal string
%     hex_str = string(data{2,i});
%     hex_str = regexprep(hex_str, '[^0-9A-Fa-f]', '');
% 
%     % Pad to 16 hexadecimal digits
%     hex_str = pad(hex_str, 16, 'left', '0');
% 
%     % Split into two 32-bit parts
%     high = uint64(hex2dec(extractBetween(hex_str, 1, 8)));
%     low  = uint64(hex2dec(extractBetween(hex_str, 9, 16)));
% 
%     % Combine into 64-bit integer
%     value = bitor(bitshift(high, 32), low);
% 
%     % Interpret as signed 64-bit two's complement
%     H_matched(i-1) = typecast(value, 'int64');
% 
% end
% 
% %% Read exported SimVision data H_plus
% data = readtable(H_plus_file, ...
%     'VariableNamingRule', 'preserve');
% 
% % Keep only one row per timestamp
% % [~, idx] = unique(data.Var1, 'stable');
% data = data(2,:);
% 
% % Remove parentheses
% data.Var2 = erase(string(data.Var2), '(');
% data.Var257 = erase(string(data.Var257), ')');
% 
% timestamps = data.Var1;
% 
% H_plus = [ ...
%     str2double(data.Var2), ...
%     data{:,3:256}, ...
%     str2double(data.Var257) ...
% ];
% 
% %% Read exported SimVision data H_min
% data = readtable(H_min_file, ...
%     'VariableNamingRule', 'preserve');
% 
% % Keep only one row per timestamp
% [~, idx] = unique(data.Var1, 'stable');
% data = data(idx,:);
% 
% % Remove parentheses
% data.Var2 = erase(string(data.Var2), '(');
% data.Var257 = erase(string(data.Var257), ')');
% 
% timestamps = data.Var1;
% 
% H_min = [ ...
%     str2double(data.Var2), ...
%     data{:,3:256}, ...
%     str2double(data.Var257) ...
% ];
% 
% %%
% H_diff = H_plus.' - H_min.';
% 
% L_algo = 256;
% lambda_reg = 0.1;
% 
% error_matrix = 0.5 * eye(L_algo) + ...
%     -0.25 * diag(ones(L_algo-1,1), 1) + ...
%     -0.25 * diag(ones(L_algo-1,1), -1);
% 
% H_delta_matched_full = (error_matrix + lambda_reg*eye(L_algo))\H_diff;
% 
% %%
% 
% %% Compare H_matched (SimVision) vs H_delta_matched_full (MATLAB)
% 
% t_idx = size(H_delta_matched_full, 2);   % timestamp column to compare (last by default)
% scale = 1;                               % set to 2^frac_bits if H_matched is fixed-point
% 
% H_matched_vec = double(H_matched(:)) / scale;
% H_delta_vec   = H_delta_matched_full(:, t_idx);
% 
% % Sanity check on lengths
% assert(numel(H_matched_vec) == numel(H_delta_vec), ...
%     'Length mismatch: H_matched has %d elements, H_delta has %d.', ...
%     numel(H_matched_vec), numel(H_delta_vec));
% 
% bins = 1:L_algo;
% 
% figure('Name', 'Matched filter comparison', 'Position', [100 100 1200 400]);
% 
% ax1 = subplot(1,2,1);
% bar(bins, H_matched_vec);
% title('H\_matched (SimVision)');
% xlabel('Bin'); ylabel('Value');
% xlim([0 L_algo+1]); grid on;
% 
% ax2 = subplot(1,2,2);
% bar(bins, H_delta_vec);
% title(sprintf('H\\_delta\\_matched\\_full (MATLAB, t = %g)', timestamps(t_idx)));
% xlabel('Bin'); ylabel('Value');
% xlim([0 L_algo+1]); grid on;
% 
% linkaxes([ax1 ax2], 'x');   % zooming one plot zooms the other
% 
% % Quick numeric check
% fprintf('Max abs difference: %g\n', max(abs(H_matched_vec - H_delta_vec)));

clear;
clc;
close all;

%% File names

H_plus_file = 'histogram_plus_callen10000_va.csv';
H_min_file = 'histogram_min_callen10000_va.csv';
H_matched_file = 'histogram_matched_callen10000_va.csv';

%% Read exported SimVision data

[H_plus, timestamps_plus] = read_histogram(H_plus_file);
[H_min, timestamps_min] = read_histogram(H_min_file);
[H_matched, timestamps_matched] = read_histogram(H_matched_file);

H_matched(40) = 0;
H_matched(162) = 0;

% Select last timestamp
H_plus = H_plus(end,:);
H_min = H_min(end,:);
H_matched = H_matched(end,:);

timestamps = timestamps_plus(end);

%% Calculate histogram difference

H_diff = double(H_plus(:)) - double(H_min(:));

L_algo = 256;
lambda_reg = 0.1;

error_matrix = 0.5 * eye(L_algo) ...
    - 0.25 * diag(ones(L_algo-1,1), 1) ...
    - 0.25 * diag(ones(L_algo-1,1), -1);

H_delta_matched_full = ...
    (error_matrix + lambda_reg*eye(L_algo)) \ H_diff;

%% Compare H_matched (SimVision) vs H_delta_matched_full (MATLAB)

scale = 1;  % Set to 2^frac_bits if H_matched is fixed-point

H_matched_vec = double(H_matched(:)) / scale;
H_delta_vec = H_delta_matched_full(:);

% Sanity check on lengths
assert(numel(H_matched_vec) == numel(H_delta_vec), ...
    'Length mismatch: H_matched has %d elements, H_delta has %d.', ...
    numel(H_matched_vec), numel(H_delta_vec));

bins = 0:L_algo-1;

%% Plot comparison

figure('Name', 'Matched filter comparison', ...
    'Position', [100 100 1200 400], ...
    'Color', 'w');

% SimVision result
ax1 = subplot(1,2,1);

bar(bins, H_matched_vec);

title('H_{matched} (SimVision)', 'Color', 'k');
xlabel('Bin', 'Color', 'k');
ylabel('Value', 'Color', 'k');

xlim([-1 L_algo]);
grid on;
box on;

% MATLAB result
ax2 = subplot(1,2,2);

bar(bins, H_delta_vec);

title(sprintf('H_{delta,matched} (MATLAB, t = %g)', timestamps), ...
    'Color', 'k');

xlabel('Bin', 'Color', 'k');
ylabel('Value', 'Color', 'k');

xlim([-1 L_algo]);
grid on;
box on;

% White background and black text
set([ax1 ax2], ...
    'Color', 'w', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'FontSize', 12, ...
    'LineWidth', 1);

% Use same axis limits for easier comparison
linkaxes([ax1 ax2], 'xy');

%% Numerical comparison

difference = H_matched_vec - H_delta_vec;

fprintf('\nComparison results:\n');
fprintf('Maximum absolute difference: %g\n', max(abs(difference)));
fprintf('Mean absolute difference:    %g\n', mean(abs(difference)));
fprintf('RMSE:                        %g\n', sqrt(mean(difference.^2)));

%% Plot difference

figure('Name', 'Difference between SimVision and MATLAB', ...
    'Position', [100 100 900 400], ...
    'Color', 'w');

plot(bins, difference, '-o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 3);

xlabel('Bin', 'Color', 'k');
ylabel('Difference', 'Color', 'k');
title('Difference between SimVision and MATLAB', 'Color', 'k');

grid on;
box on;

set(gca, ...
    'Color', 'w', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'FontSize', 12, ...
    'LineWidth', 1);

xlim([-1 L_algo]);

%% ============================================================
%  Function: Read SimVision histogram CSV
%  Supports:
%  - Old format with parentheses and Var1, Var2, ...
%  - New format with signal names as column headers
%  - Hexadecimal values
%  - Signed 64-bit two's complement
%  - Multiple timestamps
% ============================================================

function [H, timestamps] = read_histogram(filename)

    %% Read CSV

    data = readtable(filename, ...
        'VariableNamingRule', 'preserve');

    % First column contains timestamps
    timestamps = data{:,1};

    %% Remove duplicate timestamps

    [~, idx] = unique(timestamps, 'stable');

    data = data(idx,:);
    timestamps = timestamps(idx);

    %% Determine histogram size

    N_bins = width(data) - 1;
    N_times = height(data);

    % Allocate signed 64-bit output
    H = zeros(N_times, N_bins, 'int64');

    %% Detect export format

    first_value = strtrim(string(data{1,2}));

    % Old format has parentheses
    old_format = startsWith(first_value, '(');

    %% Convert histogram values

    for i = 1:N_times

        for j = 2:width(data)

            % Extract value
            value_str = strtrim(string(data{i,j}));

            % Remove parentheses
            value_str = erase(value_str, '(');
            value_str = erase(value_str, ')');

            % Remove whitespace
            value_str = strtrim(value_str);

            if old_format

                %% Old format: decimal values

                H(i,j-1) = int64(str2double(value_str));

            else

                %% New format: hexadecimal values

                % Remove non-hexadecimal characters
                value_str = regexprep(value_str, '[^0-9A-Fa-f]', '');

                % Check hexadecimal length
                if strlength(value_str) > 16
                    error('Hexadecimal value exceeds 64 bits.');
                end

                % Pad to 16 hexadecimal digits
                value_str = pad(value_str, 16, 'left', '0');

                % Split into two 32-bit parts
                high = uint64(hex2dec(extractBetween(value_str, 1, 8)));
                low  = uint64(hex2dec(extractBetween(value_str, 9, 16)));

                % Combine into 64-bit unsigned integer
                value = bitor(bitshift(high, 32), low);

                % Interpret as signed two's complement
                H(i,j-1) = typecast(value, 'int64');

            end
        end
    end

    %% Display information

    fprintf('\nLoaded: %s\n', filename);
    fprintf('Number of timestamps: %d\n', N_times);
    fprintf('Number of histogram bins: %d\n', N_bins);

    if old_format
        fprintf('Detected format: Old (decimal with parentheses)\n');
    else
        fprintf('Detected format: New (hexadecimal)\n');
    end

end