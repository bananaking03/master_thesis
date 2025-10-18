t = 0:1e-4:1;
vref = 1;
non_lin_adc = 1; % level of nonlinearity of adc;

lin_output_flash_adc = flash_adc(t, 5,vref, non_lin_adc);
figure;
plot(t,lin_output_flash_adc);