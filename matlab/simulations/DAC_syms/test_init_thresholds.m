N_bits = 12;

caps = derrive_caps(N_bits,1*sqrt(2^,10^(-15));

thresholds = cap_to_thr(caps,-1,1);

figure;
% samples = -1:2/(2^N_bits):1;
plot(samples,thresholds,'*')
hold on
fplot(@(x) x);
xlim([-1 1])