clear;
clc;
close all;

filename = "resolution_time.matlab";

%% Load exported Cadence data

data = readtable(filename, ...
    'FileType', 'text', ...
    'Delimiter', ',', ...
    'VariableNamingRule', 'preserve');

V_diff = data{:,1};
res_time = data{:,2};

%% Plot resolution time vs differential input voltage

figure('Color', 'w');

plot(V_diff * 1e3, res_time * 1e12, '-o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6);

xlabel('Differential Input Voltage (mV)', 'Color', 'k');
ylabel('Resolution Time (ps)', 'Color', 'k');
title('Comparator Resolution Time vs Differential Input Voltage', 'Color', 'k');

set(gca, ...
    'Color', 'w', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'FontSize', 12, ...
    'LineWidth', 1);

xlim([min(V_diff), max(V_diff)] * 1e3);

%% Zoomed-in plot around V_diff = 0

figure('Color', 'w');

plot(V_diff * 1e3, res_time * 1e12, '-o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6);

xlabel('Differential Input Voltage (mV)', 'Color', 'k');
ylabel('Resolution Time (ps)', 'Color', 'k');
title('Comparator Resolution Time (Zoomed)', 'Color', 'k');

grid on;
box on;

set(gca, ...
    'Color', 'w', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'FontSize', 12, ...
    'LineWidth', 1);

% Zoom into the region around zero
xlim([-10 10]);

%% Zoomed-in plot around V_diff = 0

figure('Color', 'w');

plot(V_diff * 1e3, res_time * 1e12, '-o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6);

xlabel('Differential Input Voltage (mV)', 'Color', 'k');
ylabel('Resolution Time (ps)', 'Color', 'k');
title('Comparator Resolution Time (Zoomed)', 'Color', 'k');

grid on;
box on;

set(gca, ...
    'Color', 'w', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'FontSize', 12, ...
    'LineWidth', 1);

% Zoom into the region around zero
xlim([-5 5]);