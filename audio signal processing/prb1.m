% problem Set 10 8.37

n = [0:1:255];
WN =exp(-j*2*pi/256);

WN1 = WN.^31*n;
WN2 = WN.^-32*n;

WN3 = WN.^-1*n;

hn = (1/256)*(WN1-WN2)/(1-WN3);


stem(256,hn)
title('Impulse response')