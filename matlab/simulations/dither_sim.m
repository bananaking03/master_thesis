L = 2^3;
Vinc = Vref / L; % 1 LSB

init_thresholds = linspace(0, Vref, L+1)'; % L+1 edges
init_thresholds(4) = init_thresholds(4) + Vinc*0.9;
init_thresholds = init_thresholds(2:end-1); % remove 0 and Vref
    
t = 0:0.001:1;
% add dither
t_dith = t + Vinc;

output = flash_adc(t,3,1,0,init_thresholds);
output_dith = flash_adc(t_dith,3,1,0,init_thresholds);
output_dith_comp = output_dith - 1;

figure;
plot(output, 'LineWidth', 2)
hold on
plot(output_dith_comp, 'g--', 'LineWidth', 2)

% Add legend
legend('Original output', 'Dithered output', 'Dithered output compensated')
xlabel('Time index')
ylabel('ADC output code')
title('Flash ADC Output with Dithering')
grid on