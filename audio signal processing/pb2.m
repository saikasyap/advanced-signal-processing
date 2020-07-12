%problem 8.39
n = [0:1:31];
x1 = cos(0.2*pi*n);
x2 = sin(0.3*pi*n);
%upsampling
x3 = upsample(x1,2);
x4 = upsample(x2,2);


%ploting y 

X1 = fft(x1,32); X2 = fft(x2,32);
Y = X1.*X2;
y = ifft(Y,32);


subplot(2,1,1);
stem(y)
title('Y response')

N = 64;

X3 = fft(x3,N); X4 = fft(x4,N);
X5 = X3.*X4;
x5 = ifft(X5,N);

subplot (2,1,2);
stem(x5)
title('X5 response')

