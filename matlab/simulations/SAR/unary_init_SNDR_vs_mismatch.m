cal_len = 5000;
% N = 5001;
% fs = 2000;
% f0 = (123/N) * fs;
N = (2048*2^-3)-1; % fft size
fs = 48000;  % coherent sampling
f0 = (13/N)*fs;
LSC = logspace(-6,-20,16);
Vlow= -1;
Vhigh = 1;
N_bits = 10;

A_vt = 2*10^-5; % technology constant
cap_dens = 1; % F per m^2
AreaLSC = LSC./cap_dens;
sigmaLSC = A_vt*LSC0./sqrt(AreaLSC);
mismatch_LSC = A_vt./sqrt(AreaLSC);
                                                                % seperate
                                                                % DAC SNDR
                                                                % from ADC
                                                                % SNDR
% analog_in = 0.5 + 0.5*sin(2*pi*20*t);
t = (0:1/fs:(cal_len-1)/fs)';
analog_in =sin(2*pi*f0*t);

SNDRs = zeros(1,length(LSC));

for i = 1:length(LSC)
    % Create the caps
    caps = derrive_caps_unary(2^N_bits,mismatch_LSC(i),LSC(i));

    thresholds = cap_to_thr_unary(caps,Vlow,Vhigh);

    analog_out = flash_adc(analog_in,N_bits,Vhigh,Vlow,thresholds);

    SNDRs(i) = calculate_SNDR(analog_out,analog_in,N);
    
    % % thresholds from caps
    % init_thresholds = cap_to_thr(caps(1:N_bits+1),Vlow,Vhigh);
    % 
    % [initial_digi_out] = flash_adc(analog_in,N_bits,Vhigh,Vlow,init_thresholds);
    % 
    % % SNDRs(N_bits,i) = calculate_SNDR(initial_digi_out,analog_in,N);
    % % convert codes to analog using ideal DAC
    % ideal_dac_out = Vlow + (initial_digi_out/(2^N_bits - 1))*(Vhigh - Vlow);
    % 
    % % compute SNDR using analog signals
    % SNDRs(N_bits,i) = calculate_SNDR(ideal_dac_out,analog_in,N);
end

figure;
semilogx(LSC,SNDRs);
hold on;
grid on;
xlabel('LSC');
ylabel('SNDR (dB)');
title('SNDR vs LSC capacitance');
ylim([0 80])
% add a legend
% legend(arrayfun(@(n) sprintf('N_{bits} = %d', n), 1:12, 'UniformOutput', false));