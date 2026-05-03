function [thresholds, dac_levels] = cap_to_thr(caps, Vlow, Vhigh)

% caps  : capacitor array [C_MSB ... C_LSB]
% Vlow  :
% Vhigh :
%
% thresholds : decision thresholds
% dac_levels : differential DAC levels

N = length(caps) -1;
L = 2^N;

Ctot = sum(caps);
Vrange = Vhigh - Vlow;

dac_levels = zeros(L,1);

for code = 0:L-1
    
    bits = bitget(code, N:-1:1);   % MSB first
    
    % Differential bottom-plate switching model
    charge = sum(bits .* caps(2:end));
    
    dac_levels(code+1) = Vrange * charge / Ctot;

    if (Ctot == 0)
        dac_levels(code+1) = 0;
    end
end

% Sort in case mismatch breaks monotonicity
dac_levels = sort(dac_levels);

% Thresholds = midpoints between adjacent codes
% thresholds = (dac_levels(1:end-1) + dac_levels(2:end)) / 2;
% 
% % add Vlow and Vhigh
% thresholds = [Vlow; thresholds; Vhigh];

thresholds = dac_levels(:);

end