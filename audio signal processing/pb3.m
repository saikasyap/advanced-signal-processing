% problem 8.52
n = [0:1:7];
xn = cos(0.33*pi*n);


XN = fft(xn,8);


wn = 2*xn;

YK = 2*XN;
yn = xn+ cos(0.33*pi*(n-4));

yk = ifft(YK,n);

subplot(3,1,1);
stem(xn)
title ('xn');

subplot(3,1,2);
stem(wn)
title ( 'wn')


subplot(3,1,3);
stem(yn)
title ( 'yn')


