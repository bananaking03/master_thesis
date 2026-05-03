cal_len = 5000;
% N = 5001;
% fs = 2000;
% f0 = (123/N) * fs;
N = (2048*2^-3)-1; % fft size
fs = 48000;  % coherent sampling
f0 = (13/N)*fs;
LSC = 10^-15;
                                                                % seperate
                                                                % DAC SNDR
                                                                % from ADC
                                                                % SNDR
% analog_in = 0.5 + 0.5*sin(2*pi*20*t);
t = (0:1/fs:(cal_len-1)/fs)';
analog_in =sin(2*pi*f0*t);

SNDRs = zeros(12,1001);

figure;
hold on;
grid on;
xlabel('Mismatch LSC');
ylabel('SNDR (dB)');
title('SNDR vs Mismatch LSC for Varying N_bits');
ylim([0 80])
for N_bits = 1:12
    i = 1;
    for mismatch_LSC = 0:0.001:1
        % Create the caps
        caps = derrive_caps(N_bits,mismatch_LSC,LSC);

        % dac_levels = cap_to_thr(caps(1:N_bits+1),Vlow,Vhigh);
%%
        L = 2^N_bits;
        Ctot = sum(caps);
        Vrange = Vhigh - Vlow;
        
        dac_levels = zeros(L,1);
        
        for code = 0:L-1
            
            bits = bitget(code, (length(caps)-1):-1:1);   % MSB first
            
            % Differential bottom-plate switching model
            charge = sum((bits - 0.5) .* caps(2:end));
            
            dac_levels(code+1) = Vrange * charge / Ctot;
            
        end
        
        % Sort in case mismatch breaks monotonicity
        dac_levels = sort(dac_levels);
%%
        analog_out = zeros(size(analog_in));

        for k = 1:length(analog_in)
        
            idx = find(analog_in(k) < init_thresholds,1);
        
            if isempty(idx)
                idx = length(dac_levels);
            end
        
            idx = min(idx,length(dac_levels));   % <-- prevents overflow
        
            analog_out(k) = dac_levels(idx);
        
        end

        SNDRs(N_bits,i) = calculate_SNDR(analog_out,analog_in,N);
        
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
        i = i+1;
    end
    plot(0:0.001:1,SNDRs(N_bits,:));
end
% add a legend
legend(arrayfun(@(n) sprintf('N_{bits} = %d', n), 1:12, 'UniformOutput', false));