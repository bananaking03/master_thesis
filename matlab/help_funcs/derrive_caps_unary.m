function [caps] = derrive_caps_unary(amount_of_caps,mismatch,LSC)
%DERRIVE_CAPS_UNIRARY Summary of this function goes here
%   Detailed explanation goes here
% sigma = LSC*mismatch;
% 
% caps = LSC.*ones(1,amount_of_caps) + sigma .* randn([1, amount_of_caps]);
% 
% caps(caps < 0) = LSC/1000;

sigma = LSC * mismatch;

caps = zeros(1, amount_of_caps);

for k = 1:amount_of_caps
    new_val = -1;   % force first pass into loop
    while new_val < 0
        new_val = LSC + sigma * randn();
    end
    caps(k) = new_val;
end

end