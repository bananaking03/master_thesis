function [caps] = derrive_caps(N_bits,mismatch_LSC,LSC)
%DERRIVE_CAPS Summary of this function goes here
%   LSC = least significant capacitor
%   Ac = technoligy constant

caps = zeros(1,N_bits+1);

caps(1) = LSC;
caps(2) = LSC;
for i = 2:N_bits
    caps(i+1) = LSC *2^(i-1);
end

C1 = caps(1);

% Derive mismatch for all capacitors
mismatch = mismatch_LSC * sqrt(C1 ./ caps);

% Convert to sigma
sigma_C = caps .* mismatch;

% Add Gaussian noise
% caps = caps + sigma_C .* randn(size(caps));
% 
% % clip caps so thay can't be negative
% caps(caps < 0) = 0;

% Add Gaussian noise, retry until non-negative
for k = 1:length(caps)
    new_val = -1;   % force first entry into loop
    while new_val < 0
        new_val = caps(k) + sigma_C(k) * randn();
    end
    caps(k) = new_val;
end

end