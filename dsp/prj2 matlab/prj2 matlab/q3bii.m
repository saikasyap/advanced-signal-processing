%question 3 b (ii) 

[b_e a_e] = ellip(3,1,50,[0.25 0.30],'bandpass');
v = [b_e a_e];
v_q = quan2N(11,v);
b_q = v_q(1:7);
a_q = v_q(8:14);

[h,w] = freqz(b_q,a_q);
[h1,w1] = freqz(b_e,a_e);
figure()

%unquantized filter freq response
plot(w1/(2*pi),20*log10(abs(h1)),'-.r')
%magnitude
hold on
%quantized filter freq response
plot(w/(2*pi),20*log10(abs(h)))
title('magnitude of Frequency response');
xlabel('Frequency [Hz]'), ylabel('Amplitude [dB]')
 axis([0 0.5 -55 5]), grid
legend('unquan =20*log10(abs(h1))','quan = 20*log10(abs(h))');

 