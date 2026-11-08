
clear;
clc;
close all;

%% Parameters
Vmax = 1;
N = 1000;

Vin = linspace(-1.5*Vmax, 1.5*Vmax, N);

%% 1. Ideal linear system
Vout_linear = Vin;

%% 2. Nonlinear system (third-order distortion)
a1 = 0.6;
a3 = 0.4;

Vout_nonlinear = a1*Vin + a3*Vin.^3;

%% 3. Clipping system
Vclip = Vmax;

Vout_clipped = min(max(Vin, -Vclip), Vclip);

%% Plot transfer characteristics
figure('Color', 'w');
hold on;

plot(Vin, Vout_linear, 'b', 'LineWidth', 2);
plot(Vin, Vout_nonlinear, 'r', 'LineWidth', 2);
plot(Vin, Vout_clipped, 'm', 'LineWidth', 2);

grid on;
box on;

xlabel('Input voltage');
ylabel('Output voltage');
title('Transfer Characteristics', 'Color', 'k');

lgd = legend('Ideal linear', ...
       'Third-order distortion', ...
       'Clipping', ...
       'Location', 'northwest');

set(lgd, 'TextColor', 'k', ...
         'Color', 'w', ...
         'EdgeColor', 'k');

xlim([-1.5 1.5]);
ylim([-1.25 1.25]);

% White background and black axes/text
set(gca, 'Color', 'w', ...
         'XColor', 'k', ...
         'YColor', 'k', ...
         'FontSize', 12);

set(gcf, 'Color', 'w');


%% adc version

%% Parameters
Vmax = 1;                 % ADC input range: -Vmax to +Vmax
N_bits = 4;               % ADC resolution
N_levels = 2^N_bits;      % Number of quantization levels
N = 10000;

Vin = linspace(-1.5*Vmax, 1.5*Vmax, N);

LSB = 2*Vmax/N_levels;

%% 1. Linear system + ADC
Vout_linear = Vin;

% ADC quantization
adc_input_linear = min(max(Vout_linear, -Vmax), Vmax);
code_linear = floor((adc_input_linear + Vmax)/LSB);
code_linear = min(max(code_linear, 0), N_levels-1);

Vout_linear_adc = -Vmax + (code_linear + 0.5)*LSB;

%% 2. Nonlinear system + ADC
a1 = 0.6;
a3 = 0.4;

Vout_nonlinear = a1*Vin + a3*Vin.^3;

% ADC quantization
adc_input_nonlinear = min(max(Vout_nonlinear, -Vmax), Vmax);
code_nonlinear = floor((adc_input_nonlinear + Vmax)/LSB);
code_nonlinear = min(max(code_nonlinear, 0), N_levels-1);

Vout_nonlinear_adc = -Vmax + (code_nonlinear + 0.5)*LSB;

%% 3. Clipping system + ADC
Vclip = 0.75*Vmax;

Vout_clipped = min(max(Vin, -Vclip), Vclip);

% ADC quantization
adc_input_clipped = min(max(Vout_clipped, -Vmax), Vmax);
code_clipped = floor((adc_input_clipped + Vmax)/LSB);
code_clipped = min(max(code_clipped, 0), N_levels-1);

Vout_clipped_adc = -Vmax + (code_clipped + 0.5)*LSB;

%% Plot transfer characteristics
figure('Color', 'w');
hold on;

stairs(Vin, Vout_linear_adc, 'b', 'LineWidth', 1.8);
stairs(Vin, Vout_nonlinear_adc, 'r', 'LineWidth', 1.8);
stairs(Vin, Vout_clipped_adc, 'm', 'LineWidth', 1.8);

grid on;
box on;

xlabel('Input voltage', 'Color', 'k');
ylabel('ADC output voltage', 'Color', 'k');
title('ADC Transfer Characteristics', 'Color', 'k');

lgd = legend('Linear ADC', ...
             'Third-order nonlinear ADC', ...
             'Clipped ADC', ...
             'Location', 'northwest');

set(lgd, 'TextColor', 'k', ...
         'Color', 'w', ...
         'EdgeColor', 'k');

set(gca, 'Color', 'w', ...
         'XColor', 'k', ...
         'YColor', 'k', ...
         'FontSize', 12);

set(gcf, 'Color', 'w');

xlim([-1.5 1.5]);
ylim([-1.1 1.1]);

