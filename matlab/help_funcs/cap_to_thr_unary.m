% function [thresholds] = cap_to_thr_unary(caps, Vlow, Vhigh)
% %CAP_TO_THR_UNARY Summary of this function goes here
% %   Detailed explanation goes here
% Vrange = Vhigh - Vlow;
% C_tot = sum(caps);
% 
% thresholds = zeros(length(caps)+1,1);
% thresholds(1) = Vlow;
% for i = 1:length(caps)
%     thresholds(i+1) = Vlow + Vrange * sum(caps(1:i))/C_tot;
% end
% end

function [thresholds] = cap_to_thr_unary(caps, Vlow, Vhigh)
%CAP_TO_THR_UNARY Unary capacitor DAC thresholds centered around midpoint

Vrange = Vhigh - Vlow;
Vmid   = (Vhigh + Vlow) / 2;

C_tot = sum(caps);

thresholds = zeros(length(caps)+1,1);

if (C_tot == 0)
    thresholds(:) = Vmid;
    return;
end

for i = 0:length(caps)

    % Normalized cumulative charge
    frac = sum(caps(1:i)) / C_tot;

    % Center thresholds around midpoint
    thresholds(i+1) = Vmid + Vrange * (frac - 0.5);
end

end

% function [thresholds] = cap_to_thr_unary(caps, Vlow, Vhigh)
% 
% Vrange = Vhigh - Vlow;
% C_tot  = sum(caps);
% 
% Ncodes = length(caps) + 1;
% 
% thresholds = zeros(Ncodes,1);
% 
% if (C_tot == 0)
%     return;
% end
% 
% cumcaps = [0; cumsum(caps(:))];
% 
% % Raw DAC levels
% levels = Vrange * cumcaps / C_tot;
% 
% % Shift so middle code is exactly zero
% mid_idx = ceil(Ncodes/2);
% levels = levels - levels(mid_idx);
% 
% thresholds = levels;
% 
% end