% question 2 a

b = [1 -1.13]; % numerator
a =[1 -1.5 0.9];% denominator
zplane(b,a); title('H(z) Poles and zeros unquantized')
fvtool(b,a);title('H(z)using fvtool function')
freqz(b,a); title('H(z) Frequency Response')
 [z,p,k] = tf2zp(b,a)