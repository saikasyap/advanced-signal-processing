%question 3 a(ii) 
[b_b a_b] = butter(3,[0.25 0.30],'bandpass'); 

v = [b_b a_b];
v_q = quan2N(12,v);
b_q = v_q(1:7);
a_q = v_q(8:14);


[h,w] = freqz(b_q,a_q);% quantized
[h1,w1] = freqz(b_b,a_b);
% unquantized
figure()
plot(w1/(2*pi),20*log10(abs(h1)),'-.b')
%magnitude
hold on
%quantized filter freq response
plot(w/(2*pi),20*log10(abs(h)))
%unquantized filter freq response
title('magnitude of Frequency response');
xlabel('Frequency [Hz]'), ylabel('Amplitude [dB]')
 axis([0 0.5 -55 5]), grid
legend('unquan =20*log10(abs(h1))','quan = 20*log10(abs(h))');
 