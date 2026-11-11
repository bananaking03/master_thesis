
clear;
clc;
close all;

%% Parameters
N_bits = 6;
Vref = 1;
Vcm = Vref/2;

% Differential input voltages
Vin_plus  = 0.45;
Vin_minus = 0.72;

%% Determine which input is lower
if Vin_plus < Vin_minus
    active_DAC = "DAC+";
    Vlow = Vin_plus;
    Vhigh = Vin_minus;
    active_plus = true;
else
    active_DAC = "DAC-";
    Vlow = Vin_minus;
    Vhigh = Vin_plus;
    active_plus = false;
end

%% SAR conversion
% Only the DAC connected to the lower input is active
DAC_correction = zeros(1,N_bits+1);

current_code = 0;

for k = 1:N_bits

    % Binary-weighted trial correction
    trial_code = current_code + 2^(N_bits-k);

    trial_correction = trial_code/(2^N_bits)*Vref;

    % Compare the raised lower input with the higher input
    if Vlow + trial_correction <= Vhigh
        current_code = trial_code;
    end

    DAC_correction(k+1) = ...
        current_code/(2^N_bits)*Vref;
end

%% Generate comparator input waveforms
time = 0:N_bits;

% Higher input remains constant
Vhigh_plot = Vhigh*ones(size(time));

% Lower input is raised by the active DAC
Vlow_plot = Vlow + DAC_correction;

% Assign to differential inputs
if active_plus
    Vcomp_plus = Vlow_plot;
    Vcomp_minus = Vhigh_plot;
else
    Vcomp_plus = Vhigh_plot;
    Vcomp_minus = Vlow_plot;
end

%% Plot
figure('Color','w','Position',[100 100 1000 550]);
hold on;

plot(time,Vcomp_plus,'r','LineWidth',2.5);
plot(time,Vcomp_minus,'m','LineWidth',2.5);

plot(time,Vcomp_plus,'ro','MarkerFaceColor','r');
plot(time,Vcomp_minus,'mo','MarkerFaceColor','m');

grid on;
box on;

xlabel('SAR clock cycle');
ylabel('Comparator input voltage (V)');
title('Differential SAR ADC Operation','Color','k');

lgd = legend('V_{comp+}','V_{comp-}', ...
    'Location','best');

set(lgd,'Color','w', ...
        'TextColor','k', ...
        'EdgeColor','k');

set(gca,'Color','w', ...
        'XColor','k', ...
        'YColor','k', ...
        'FontSize',12);

xlim([0 N_bits]);
ylim([0 Vref*1.1]);
xticks(0:N_bits);

%% Results
fprintf('Active DAC: %s\n',active_DAC);
fprintf('Inactive DAC: %s\n', ...
    string(~active_plus) .* "DAC+" + ...
    string(active_plus) .* "DAC-");
fprintf('Initial differential voltage: %.4f V\n', ...
    abs(Vin_plus-Vin_minus));
fprintf('Final differential voltage: %.4f V\n', ...
    abs(Vcomp_plus(end)-Vcomp_minus(end)));
fprintf('Final DAC correction: %.4f V\n', ...
    DAC_correction(end));
