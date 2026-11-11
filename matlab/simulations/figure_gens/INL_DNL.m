clear;
clc;
close all;

%% ADC parameters
N_bits = 4;
N_codes = 2^N_bits;
Vmin = -1;
Vmax = 1;

LSB = (Vmax - Vmin)/N_codes;

%% Generate ideal thresholds
% Includes the two endpoints of the ADC range
ideal_thresholds = linspace(Vmin, Vmax, N_codes+1);

%% Generate non-ideal thresholds
% Errors are expressed in LSB
threshold_errors = [0.05 0.10 0.18 0.12 -0.05 ...
                   -0.20 -0.35 -0.15 0.10 0.30 ...
                    0.20 0.05 -0.10 -0.05 0.08];

actual_thresholds = ideal_thresholds;
actual_thresholds(2:end-1) = ...
    ideal_thresholds(2:end-1) + threshold_errors*LSB;

%% Calculate DNL
% DNL = (Actual code width / Ideal LSB) - 1
actual_widths = diff(actual_thresholds);

DNL = actual_widths/LSB - 1;

%% Calculate INL
% INL = (Actual threshold - Ideal threshold) / LSB
% Endpoint thresholds are excluded
INL = (actual_thresholds(2:end-1) - ...
       ideal_thresholds(2:end-1))/LSB;

%% Select examples to highlight
selected_code = 8;       % Code for DNL illustration
selected_threshold = 8; % Transition for INL illustration

% MATLAB indices
dnl_idx = selected_code + 1;
inl_idx = selected_threshold + 1;

%% Plot
figure('Color','w','Position',[100 100 1100 850]);


%% 1. ADC transfer characteristic (zoomed in)
subplot(3,1,1);
hold on;

stairs(ideal_thresholds, [0:N_codes-1 N_codes-1], ...
    '--k', 'LineWidth', 2);

stairs(actual_thresholds, [0:N_codes-1 N_codes-1], ...
    'b', 'LineWidth', 2.5);

% Highlight INL deviation
plot([ideal_thresholds(inl_idx) actual_thresholds(inl_idx)], ...
     [selected_threshold selected_threshold], ...
     'r-o', 'LineWidth', 2.5, 'MarkerFaceColor','r');

% Highlight DNL code width
y_dnl = selected_code - 0.5;

plot([actual_thresholds(dnl_idx) actual_thresholds(dnl_idx+1)], ...
     [y_dnl y_dnl], 'm', 'LineWidth', 4);

% Show ideal LSB width for comparison
plot([ideal_thresholds(dnl_idx) ideal_thresholds(dnl_idx+1)], ...
     [y_dnl-0.35 y_dnl-0.35], '--k', 'LineWidth', 2,'Color','yellow');

grid on;
box on;

title('ADC Transfer Characteristic (Zoomed In)','Color','k');
xlabel('Input voltage (V)');
ylabel('ADC output code');

lgd = legend('Ideal ADC','Non-ideal ADC', ...
    'INL deviation','Actual code width','Ideal code width', ...
    'Location','northwest');

set(lgd, 'Color','w', ...
         'TextColor','k', ...
         'EdgeColor','k');

% Zoom around the selected thresholds
x_center = actual_thresholds(dnl_idx);
xlim([x_center-2.5*LSB x_center+2.5*LSB]);

ylim([selected_code-2 selected_code+2]);


%% 2. Differential Non-Linearity
subplot(3,1,2);
hold on;

bar(0:N_codes-1, DNL, 'FaceColor',[0.2 0.6 0.9]);

% Highlight selected code
bar(selected_code-1, DNL(dnl_idx), ...
    'FaceColor','m');

yline(0,'k','LineWidth',1);
yline(-1,'r--','LineWidth',1.5);

text(selected_code-1, DNL(dnl_idx)+0.08, ...
    sprintf('%.2f LSB',DNL(dnl_idx)), ...
    'Color','m','HorizontalAlignment','center');

title('Differential Non-Linearity (DNL)','Color','k');
xlabel('ADC code');
ylabel('DNL (LSB)');

grid on;
box on;
xlim([-0.5 N_codes-0.5]);

%% 3. Integral Non-Linearity
subplot(3,1,3);
hold on;

bar(1:N_codes-1, INL, 'FaceColor',[0.2 0.7 0.4]);

% Highlight selected threshold
bar(selected_threshold, INL(selected_threshold), ...
    'FaceColor','r');

yline(0,'k','LineWidth',1);

text(selected_threshold, ...
    INL(selected_threshold)+0.08*sign(INL(selected_threshold)+0.001), ...
    sprintf('%.2f LSB',INL(selected_threshold)), ...
    'Color','r','HorizontalAlignment','center');

title('Integral Non-Linearity (INL)','Color','k');
xlabel('Threshold number');
ylabel('INL (LSB)');

grid on;
box on;
xlim([0 N_codes]);

%% Formatting
set(findall(gcf,'Type','axes'), ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k', ...
    'FontSize',11);

sgtitle('ADC Threshold Errors: INL and DNL','Color','k');
