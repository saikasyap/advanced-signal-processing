% questionn 4 c
% Generation of Sine Wave 

Fs = 44100;
t =0:1/Fs:400/Fs;
f= 440;

x = sin(2*pi*f*t);
% using plot function to plot the signal
figure;
plot(x);
title('Sine wave versus time');
xlabel('time');
ylabel('Magnitude');
