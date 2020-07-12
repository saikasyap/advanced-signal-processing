clear all
a1 = [1 -(0.5-1i*0.5)]; a2 = [1 -(0.5+1i*0.5)]; a = conv(a1,a2); 
b1 = [1 (1-1i*0.75)/1.25]; b2 = [1 (1+1i*0.75)/1.25]; b = conv(b1,b2); 
zplane(b,a); title(' Hmin(z) Poles and zeros') 
figure 
freqz(b,a); title('Hmin(z) Frequency Response') 
