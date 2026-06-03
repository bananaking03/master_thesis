function dither = gen_dither_sequence(n, seed)
%GEN_DITHER_SEQUENCE Generate the same dither sequence as dither_gen.v
%   dither = gen_dither_sequence(n) returns an n-element vector of 0/1
%   values matching the Verilog LFSR used in dither_gen.v.
%
%   dither = gen_dither_sequence(n, seed) lets you override the 16-bit
%   LFSR seed. The default seed is hex2dec('ACE1'), which matches the
%   Verilog reset value.

    if nargin < 2
        seed = hex2dec('ACE1');
    end
    assert(isscalar(n) && n >= 0 && n == floor(n), 'n must be a non-negative integer');
    assert(isscalar(seed) && seed == floor(seed) && seed >= 0 && seed < 2^16, 'seed must be a 16-bit integer');

    lfsr = uint16(seed);
    dither = false(1, n);

    for k = 1:n
        % output the current LSB first, then update the LFSR
        dither(k) = bitand(lfsr, uint16(1));
        newbit = bitxor(bitxor(bitxor(bitget(lfsr, 16), bitget(lfsr, 14)), bitget(lfsr, 13)), bitget(lfsr, 11));
        lfsr = bitshift(lfsr, -1);
        lfsr = bitor(lfsr, bitshift(uint16(newbit), 15));
    end
end
