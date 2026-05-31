function y = shift_ones(x, shifts)
%SHIFT_ONES Shift the positions of 1's in a vector
%
% x       : input binary vector
% shifts  : shift amount for each 1 in x
%            positive = right
%            negative = left
%
% y       : output vector

N = length(x);

% Find positions of the 1's
idx = find(x == 1);

% Check matching lengths
if length(idx) ~= length(shifts)
    error('Number of shifts must match number of 1s in x');
end

% Apply shifts
new_idx = idx + shifts(:).';

% Keep only valid indices
valid = (new_idx >= 1) & (new_idx <= N);
new_idx = new_idx(valid);

% Create output vector
y = zeros(size(x));
y(new_idx) = 1;

end