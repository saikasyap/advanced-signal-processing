% questionn 4 d
% Generation of Sine Wave 
r =0.8;

Fs = 44100;
t =0:1/Fs:400/Fs;
f= 4040;
% Wo value 
Wo = 2*pi*f/Fs;

x = sin(2*pi*f*t);

%using filter function 

d = [1 -2*cos(Wo) 1];
c = [1 -2*r*cos(Wo) r^2];

y = filter(d,c,x);
figure;
plot(y);
hold on;
plot(x);
plot(y);
xlabel('time');
ylabel('Magnitude');
title('output of the filter on the same axes as the input signal using the hold command')    
 % Plot the X values vs. the Y values 


