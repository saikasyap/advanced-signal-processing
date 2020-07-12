%question 5 d and e
% defining values of r :
r=0.9;
%r=0.95
%r=0.99
%r=0.999
load('River_Short.mat')
%intiating angular frerquency values for triple notch
Wo = 2*pi*440/fs;
Wi = 2*pi*(440+88)/fs;
Wii = 2*pi*(440-88)/fs;
load('River_Short.mat')
%soundsc(zp,fs)
a = [1 -2*cos(Wo) 1];
b = [1 -2*r*cos(Wo) r^2];
% filtering out the sound using notch filter
z1 = filter(a,b,zp);
c = [1 -2*cos(Wi) 1];
%soundsc(z1,fs)
d = [1 -2*r*cos(Wi) r^2];
z2 = filter(c,d,z1);
%soundsc(z2,fs)

%soundsc(z1,fs)
e = [1 -2*cos(Wii) 1];
f = [1 -2*r*cos(Wii) r^2];
z = filter(e,f,z2);
% hearing the filitered signal
soundsc(z,fs)
N = length(xp);
Zp = fft(zp,N);
%spectrum of Zp signal
% figure;
% plot(20*log10(abs(Zp(1:N/8))))
% title('zp signal :spectrum');
% figure;
%spectrum of the signal after triple notch 
bin_vals = [0 : N-1];
fax_Hz = bin_vals*Fs/N;
N_2 = ceil(N/20);
Z = abs(fft(z));
plot(fax_Hz(1:N_2),20*log10(Z(1:N_2)))
title('z signal : filtered spectrum')
grid