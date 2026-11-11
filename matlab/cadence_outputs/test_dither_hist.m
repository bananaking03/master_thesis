%% Test rather the outputs of the ADC are going into the correct histogram bins in verilog

result_file = 'result_dith_callen10000_dis_va';
dither_file = 'dither_callen10000_va';
H_plus_file = 'histogram_plus_callen10000_va';
H_min_file = 'histogram_min_callen10000_va';

%% Create histograms from results
data = read_histogram(result_file);

% Extract columns
data_col1 = data(1:10000,1);
data_col2 = floor(double(data(1:10000,2))/4);

% Initialize histograms (256 bins)
H_min  = zeros(256,1);
H_plus = zeros(256,1);

% Loop over all rows except the first one
for i = 2:length(data_col2)

    % Delayed first column
    if data_col1(i-1) == 0

        % H_min
        H_min(data_col2(i)+1) = H_min(data_col2(i)+1) + 1;

    elseif data_col1(i-1) == 1

        % H_plus (subtract 1, wrapping around 0-255)
        value = mod(data_col2(i)-1,256);
        H_plus(value+1) = H_plus(value+1) + 1;

    end
end

%% Read verilog histograms

data = read_histogram(H_plus_file);

H_plus_v = data(1,:)';

data = read_histogram(H_min_file);

H_min_v = data(1,:)';

H_plus_delta = double(H_plus_v) - H_plus;
H_min_delta = double(H_min_v) - H_min;
%% Plot histograms

% Results plots
figure('Color','w');

x = 0:255;

subplot(2,1,1);
bar(x, H_min);
title('H_{min}');
xlabel('ADC code');
ylabel('Count');
grid on;
xlim([0 255]);
set(gca,'Color','w','XColor','k','YColor','k');

subplot(2,1,2);
bar(x, H_plus);
title('H_{plus}');
xlabel('ADC code');
ylabel('Count');
grid on;
xlim([0 255]);
set(gca,'Color','w','XColor','k','YColor','k');

sgtitle('Dither Histograms','Color','k');

% verilog plots
figure('Color','w');

subplot(2,1,1);
bar(x, H_min_v);
title('H_{min}');
xlabel('ADC code');
ylabel('Count');
grid on;
xlim([0 255]);
set(gca,'Color','w','XColor','k','YColor','k');

subplot(2,1,2);
bar(x, H_plus_v);
title('H_{plus}');
xlabel('ADC code');
ylabel('Count');
grid on;
xlim([0 255]);
set(gca,'Color','w','XColor','k','YColor','k');

sgtitle('Dither Histograms verilog','Color','k');

% difference plots
figure('Color','w');

subplot(2,1,1);
bar(x, H_min_delta);
title('H_{min}');
xlabel('ADC code');
ylabel('Count');
grid on;
xlim([0 255]);
set(gca,'Color','w','XColor','k','YColor','k');

subplot(2,1,2);
bar(x, H_plus_delta);
title('H_{plus}');
xlabel('ADC code');
ylabel('Count');
grid on;
xlim([0 255]);
set(gca,'Color','w','XColor','k','YColor','k');

sgtitle('Dither Histograms verilog','Color','k');