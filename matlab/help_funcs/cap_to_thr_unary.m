function [thresholds] = cap_to_thr_unary(caps, Vlow, Vhigh)
%CAP_TO_THR_UNARY Summary of this function goes here
%   Detailed explanation goes here
Vrange = Vhigh - Vlow;
C_tot = sum(caps);

thresholds = zeros(length(caps)+1,1);
thresholds(1) = Vlow;
for i = 1:length(caps)
    thresholds(i+1) = Vlow + Vrange * sum(caps(1:i))/C_tot;
end
end