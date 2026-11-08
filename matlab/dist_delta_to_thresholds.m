N = 8;
L = 2^N;

delta_reg_file = 'deltareg_15_dis_va.csv';

[delta_reg, timestamps] = read_histogram(delta_reg_file);

% Select last timestamp
delta_reg = delta_reg(end,:);

thresholds = zeros(1,L);

for i=1:length(delta_reg)
    for j=1:i
        thresholds(i) = thresholds(i) + delta_reg(j);
    end
end

bins = 1:L;

figure;
bar(bins,thresholds);













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
