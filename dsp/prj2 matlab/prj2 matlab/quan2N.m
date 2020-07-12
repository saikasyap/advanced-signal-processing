function quan_b = quan2N(N,b)
% Function quan2N quantifies filter coefficients b to N bits
% Scale b such that b/2^mb has values between -1 and +1
mb = ceil(log(max(abs(b)))/log(2));
norm_b = b/2^mb; 
% Multiply the normalized b by 2^N, round to the nearest integer,
%     then divide by 2^N to get the N-bit quantization of b
% Re-multiply by 2^mb to restore the gain to the quantized b 
quan_b = 2^mb*round(norm_b*2^N)/2^N;
