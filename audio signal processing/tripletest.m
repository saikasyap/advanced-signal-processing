% Question 5 d  :
% Triple Notch Filter Example 
j=1*i;

r = 0.9; 
Wo = 0.019*pi;
r1 = 0.9; 
Wii = 0.3*pi;
r2 = 0.9; 
Wi = 0.9*pi;

%numerator ----> a 
a1 = [1 -exp(j*Wo)]; a2 = [1 -exp(-j*Wo)];
a5 = conv(a1,a2)
a3 = [1 -exp(j*Wi)]; a4 = [1 -exp(-j*Wi)];
a6 = conv(a3,a4);
a7 = [1 -exp(j*Wii)]; a8 = [1 -exp(-j*Wii)];
a9 = conv(a7,a8);
a = conv(a9,conv(a6,a5));
% Denominator ---> b
b1 = [1 -r*exp(j*Wo)]; b2 = [1 -r*exp(-j*Wo)];
b3 = conv(b1,b2);
b4 = [1 -r1*exp(j*Wi)]; b5 = [1 -r1*exp(-j*Wi)];
b6 = conv(b4,b5);
b7 = [1 -r2*exp(j*Wii)]; b8 = [1 -r2*exp(-j*Wii)];
b9 = conv(b7,b8);
b = conv(b3,conv(b6,b9));

figure
zplane(a,b); title('Hnotch(z) Poles and zeros')
figure
freqz(a,b); title('Hnotch(z) Frequency Response')


