% question 5 d
r =0.9;
%intiating angular frequency values for triple notch
Wo = 0.019*pi;
Wi = 0.0159*pi;
Wii = 00239*pi;
load('River_Short.mat')
soundsc('River_Short')
a = [1 -2*cos(Wo) 1];
b = [1 -2*r*cos(Wo) r^2];
%filtering out the sound using notch filter
z1 = filter(a,b,zp);
soundsc('zp')

c = [1 -2*cos(Wi) 1];
d = [1 -2*r*cos(Wi) r^2];

z2 = filter(c,d,z1);
e = [1 -2*cos(Wii) 1];
f = [1 -2*r*cos(Wii) r^2];
z = filter(e,f,z2);

figure;
plot(zp);
title('zp signal');

figure;
plot(z1);
title('after single notch');

figure;
plot(z2);
title('after double notch');

figure;
plot(z);
title('output of the filter after triple notch')