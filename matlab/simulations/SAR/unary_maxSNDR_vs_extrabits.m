cal_len = 5000;
N = (2048*2^-3)-1; % fft size
fs = 48000;  % coherent sampling
f0 = (13/N)*fs;
LSC0 = 1.5*10^-20;
mismlinkatch_LSC0 = 0.25;
Vlow = -1;
Vhigh = 1;
Bits_tested = 10;
N_bits = Bits_tested;
% extra_caps_tested = 0:14;
extra_caps_tested = 0:8;

A_vt = 2*10^-5; % technology constant
cap_dens = 1; % F per m^2
AreaLSC = LSC0/cap_dens;
sigmaLSC = A_vt*LSC0/sqrt(AreaLSC);
mismatch_LSC0 = A_vt/sqrt(AreaLSC);

L = 2^N_bits;

% analog_in = 0.5 + 0.5*sin(2*pi*20*t);
t = (0:1/fs:(cal_len-1)/fs)';
analog_in =sin(2*pi*f0*t);

% for N_bits = Bits_tested
%     % ideal thresholds
ideal_thresholds = linspace(Vlow, Vhigh, L+1)';
ideal_thresholds = ideal_thresholds(2:end-1);

[ideal_digi_out] = flash_adc(analog_in,N_bits,Vhigh,Vlow,ideal_thresholds);
ideal_SNDRs = calculate_SNDR(ideal_digi_out,analog_in,N);
%     for N_bits_caps = extra_caps_tested
%         %% 
% 
%         L = 2^N_bits;
%         L_caps = 2^N_bits_caps;
% 
%         % Create the caps (decrease value of smallest capacitor)
%         caps = derrive_caps(N_bits + N_bits_caps,mismatch_LSC0 * sqrt(L_caps),LSC0/(L_caps));
% 
%         pos_thresholds = cap_to_thr(caps,Vlow,Vhigh); % possible thresholds
%         [~, idx] = min(abs(pos_thresholds(:) - ideal_thresholds(:).'), [], 1);
%         best_thresholds = pos_thresholds(idx);
% 
%         [best_digi_out] = flash_adc(analog_in,N_bits,Vhigh,Vlow,best_thresholds);
% 
%         SNDRs(N_bits,N_bits_caps+1) = calculate_SNDR(best_digi_out,analog_in,N);
%     end
% end

MC_runs = 500;
SNDRs = zeros(length(extra_caps_tested), MC_runs);
SNDRs_interpol = zeros(length(extra_caps_tested), MC_runs);
SNDRs_interpol2 = zeros(length(extra_caps_tested), MC_runs);
SNDRs_interpol3 = zeros(length(extra_caps_tested), MC_runs);

for mc = 1:MC_runs

    for N_bits_caps = extra_caps_tested

        L_caps = 2^(N_bits_caps+N_bits);
        L_caps_extra = 2^N_bits_caps;
        % 
        % % Create the caps with mismatch
        % caps = derrive_caps_unary(L_caps, ...
        %     mismatch_LSC0, LSC0);
        % 
        % pos_thresholds = cap_to_thr_unary(caps,Vlow,Vhigh);

        %---------------------------------------------------

        caps_DACun = derrive_caps_unary(L, ...
            mismatch_LSC0, LSC0);

        caps_DACbin = derrive_caps(N_bits_caps,mismatch_LSC0*sqrt(LSC0/(LSC0/(2^(N_bits_caps+1)))),LSC0/(2^(N_bits_caps+1)));

        pos_thresholds_un = cap_to_thr_unary(caps_DACun,Vlow,Vhigh);
        % pos_thresholds_un = pos_thresholds_un(:);

        pos_thresholds_bin = cap_to_thr(caps_DACbin,0,(Vhigh-Vlow)/(L));
        % pos_thresholds_bin = pos_thresholds_bin(:).';

        % C = pos_thresholds_un + pos_thresholds_bin;
        % pos_thresholds = C(:).';

        pos_thresholds = reshape(pos_thresholds_bin(:) + pos_thresholds_un(:).', 1, []);
        pos_thresholds = pos_thresholds(1:end-L_caps_extra+1);

        %---------------------------------------------------------------

        % clear pos_thresholds_bin
        % clear pos_thresholds_un

        % for i = 1:
        %     pos_thresholds((i-1)*L_caps+1:) = pos_thresholds_un + pos_thresholds_bin;
        % end

        best_thresholds = zeros(length(ideal_thresholds),1);

        [~, idx] = min(abs(pos_thresholds(:) - ideal_thresholds(:).'), [], 1);
        best_thresholds = pos_thresholds(idx);

        % interpolate the best thresholds, halve of thrs interpolated 
        interpolthr = zeros(length(best_thresholds),1);
        interpolthr(1:2:end) = best_thresholds(1:2:end);
        interpol_idx = (idx(3:2:end) + idx(1:2:end-2))./2;
        interpolthr(2:2:end) = pos_thresholds(round(interpol_idx)); 
        
        % % interpolate the best thresholds, 2 thirds thrs interpolated
        % interpolthr2 = zeros(length(best_thresholds),1);
        % interpolthr2(1:3:end) = best_thresholds(1:3:end);
        % interpol_1idx2 = idx(1:3:end-3) + (idx(4:3:end) - idx(1:3:end-3))./3;
        % interpol_2idx2 = idx(1:3:end-3) + 2.*(idx(4:3:end) - idx(1:3:end-3))./3;
        % interpolthr2(2:3:end) = pos_thresholds(round(interpol_1idx2));
        % interpolthr2(3:3:end) = pos_thresholds(round(interpol_2idx2));
        % interpolthr3(end-mod((L-1),3):end) = best_thresholds(end-mod((L-1),3):end);

        % interpolate the best thresholds, 3 fourths thrs interpolated
        interpolthr3 = zeros(length(best_thresholds),1);
        interpolthr3(1:4:end) = best_thresholds(1:4:end);
        interpol_1idx3 = idx(1:4:end-4) + (idx(5:4:end) - idx(1:4:end-4))./4;
        interpol_2idx3 = idx(1:4:end-4) + 2.*(idx(5:4:end) - idx(1:4:end-4))./4;
        interpol_3idx3 = idx(1:4:end-4) + 3.*(idx(5:4:end) - idx(1:4:end-4))./4;
        interpolthr3(2:4:end-2) = pos_thresholds(round(interpol_1idx3));
        interpolthr3(3:4:end-2) = pos_thresholds(round(interpol_2idx3));
        interpolthr3(4:4:end-2) = pos_thresholds(round(interpol_3idx3));
        % interpolthr3(end-mod((L-1),4):end) = best_thresholds(end-mod((L-1),4):end);

        [best_digi_out] = flash_adc(analog_in,N_bits,Vhigh,Vlow,best_thresholds.');
        [interpol_digi_out] = flash_adc(analog_in,N_bits,Vhigh,Vlow,interpolthr);
        % [interpol_digi_out2] = flash_adc(analog_in,N_bits,Vhigh,Vlow,interpolthr2);
        [interpol_digi_out3] = flash_adc(analog_in,N_bits,Vhigh,Vlow,interpolthr3);

        SNDRs(N_bits_caps+1, mc) = calculate_SNDR(best_digi_out,analog_in,N);
        SNDRs_interpol(N_bits_caps+1, mc) = calculate_SNDR(interpol_digi_out,analog_in,N);
        % SNDRs_interpol2(N_bits_caps+1, mc) = calculate_SNDR(interpol_digi_out2,analog_in,N);
        SNDRs_interpol3(N_bits_caps+1, mc) = calculate_SNDR(interpol_digi_out3,analog_in,N);
    end
end
%% 

%  big steps because of MSB's with mismatch that get added to other
%  thresholds
% figure;
% hold on
% for i = Bits_tested
%     plot(extra_caps_tested, SNDRs(i,:));
%     yline(ideal_SNDRs(i));
% end
% xlabel('amount of extra capacitors');
% ylabel('SNDR (dB)')
% title('SNDR vs amount of capacitors')
% hold off;
% grid on;
% legend(arrayfun(@(x) sprintf('N_{bits} = %d', x), Bits_tested, 'UniformOutput', false));
%% 

figure;
% boxchart(repelem(extra_caps_tested,MC_runs), SNDRs(:))
boxchart(SNDRs(:,1:119)')
hold on
boxchart(SNDRs_interpol(:,1:119)')
% boxchart(SNDRs_interpol2')
boxchart(SNDRs_interpol3(:,1:119)')
xlabel('2log(capacitors)+1')
ylabel('SNDR (dB)')
title('SNDR vs amount of capacitors')
yline(ideal_SNDRs);

% legend('best SNDRs','interpolated SNDRs', '2/3 interpolated', '3/4 interpolated', 'max SNDR')
legend('best SNDRs','interpolated SNDRs', '3/4 interpolated', 'max SNDR')
grid on
hold off