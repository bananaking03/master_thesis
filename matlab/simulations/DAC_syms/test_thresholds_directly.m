cap_tot = 10^-12;

N_bits = 10;
N_bits_caps = 8; 
L_caps = 2^10;
cap_dens = 2;
A_vt = 10^-12;
Vlow = -1;

LSC0 = cap_tot/L_caps;
AreaLSC = LSC0/cap_dens;
sigmaLSC = A_vt*LSC0/sqrt(AreaLSC);
mismatch_LSC0 = A_vt/sqrt(AreaLSC);
mismatch_LSC0 = 0;
ideal_thr = linspace(-1,1,2^10);
% % Create the caps with mismatch
L_caps = 2^(N_bits_caps+N_bits);
caps_full_un = derrive_caps_unary(L_caps, ...
    mismatch_LSC0, LSC0);
pos_thresholds_full_un = cap_to_thr_unary(caps_full_un,Vlow,Vhigh).';

caps_DACun = derrive_caps_unary(L, ...
            mismatch_LSC0, LSC0);

L_caps_extra = 2^N_bits_caps;
% caps_DACbin = derrive_caps(N_bits_caps,mismatch_LSC0*sqrt(LSC0/(8*LSC0/(2^(N_bits_caps)))),10240*LSC0/(2^(N_bits_caps)));
caps_DACbin = derrive_caps_unary(L_caps_extra,mismatch_LSC0*sqrt(LSC0/(LSC0/(2^(N_bits_caps+1)))),LSC0/(2^(N_bits_caps+1)));

% pos_thresholds_un = cap_to_thr_unary(caps_DACun,Vlow,Vhigh);
caps = caps_DACun;

Vrange = Vhigh - Vlow;
Vmid   = (Vhigh + Vlow) / 2;

C_tot = sum(caps);

pos_thresholds_un = zeros(length(caps)+1,1);

if (C_tot == 0)
    pos_thresholds_un(:) = Vmid;
    return;
end

for i = 0:length(caps)

    % Normalized cumulative charge
    frac = sum(caps(1:i)) / C_tot;

    % Center thresholds around midpoint
    pos_thresholds_un(i+1) = Vmid + Vrange * (frac - 0.5);
end

% pos_thresholds_bin = cap_to_thr(caps_DACbin,-(Vhigh-Vlow)/(2*L),(Vhigh-Vlow)/(2*L));
pos_thresholds_bin = cap_to_thr_unary(caps_DACbin,-1*(Vhigh-Vlow)/(2*L),1*(Vhigh-Vlow)/(2*L));
% pos_thresholds_bin = flip(pos_thresholds_bin);
pos_thresholds_bin = pos_thresholds_bin(:).';

pos_thresholds = reshape(pos_thresholds_bin(1:end-1) + pos_thresholds_un(:), 1, []);
pos_thresholds = pos_thresholds(1:end-L_caps_extra+1);

if (N_bits_caps == 0)
    pos_thresholds = pos_thresholds_un.';
end

pos_thresholds_half_un = sort(pos_thresholds);

pos_thresholds_half_un = pos_thresholds_half_un(pos_thresholds_half_un >= -1 & pos_thresholds_half_un <= 1);
diff = pos_thresholds_half_un - pos_thresholds_full_un;

disp(sum(diff));

%% 
figure;
[~, idx] = min(abs(pos_thresholds_half_un(:) - ideal_thr(:).'), [], 1);
        best_thresholds_half_un = pos_thresholds_half_un(idx);
diff = ideal_thr - best_thresholds_half_un;
plot(diff);

figure;
[~, idx] = min(abs(pos_thresholds_full_un(:) - ideal_thr(:).'), [], 1);
        best_thresholds_full_un = pos_thresholds_full_un(idx);
diff = ideal_thr - best_thresholds_full_un;
plot(diff);
