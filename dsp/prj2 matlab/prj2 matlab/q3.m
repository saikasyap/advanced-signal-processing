%question 3
%question a(i)

[b_b a_b] = butter(3,[0.25 0.30],'bandpass'); 

[h,w] = freqz(b_b,a_b);
%magnitude
figure()
%unquantized filter freq response
subplot(2,1,1)
plot(w/(2*pi),20*log10(abs(h)))
title('magnitude of Frequency response');
xlabel('Frequency [Hz]'), ylabel('Amplitude [dB]')
axis([0 0.5 -55 5]), grid
%phase
subplot(2,1,2)
plot(w/(2*pi), 360/(2*pi)*angle(h))
xlabel('Frequency[Hz]'), ylabel('phase[degrees]')
axis([0 0.5 -100 100]), grid