% Generation of Sine Wave 
r =0.9;
Wo = 0.019*pi;
Wi = 0.0159*pi;
Wii = 00239*pi;
Fs = 44100;
t =0:1/Fs:400/Fs;
f= 440;
x = sin(2*pi*f*t);
% using plot function to plot the signal
%figure;
%plot(x);

%using filter function 

a = [1 -2*cos(Wo) 1];
b= [1 -2*r*cos(Wo) r^2];

y1 = filter(a,b,x);
c = [1 -2*cos(Wi) 1];
d = [1 -2*r*cos(Wi) r^2];
y2 = filter(c,d,y1);
e = [1 -2*cos(Wii) 1];
f = [1 -2*r*cos(Wii) r^2];
y = filter(e,f,y2);

figure;
plot(y);
title('output of the filter')
