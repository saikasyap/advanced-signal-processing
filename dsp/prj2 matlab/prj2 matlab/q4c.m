% numerator coeff
b_11 =  1.0000;
b_12 = 2.0000 - 0.0000i;
b_13 = 1.0000 - 0.0000i;
b_21 = 1.0000;
b_22 = 0.0000*i;
b_23 = -1.0000 - 0.0000i;
b_31 = 1.0000;
b_32 = -2.0000;
b_33 = 1.0000;

c_1 = [b_11 b_12 b_13];
c_2 = [b_21 b_22 b_23];
c_3 = [b_31 b_32 b_33];

num = conv(c_1,conv(c_2,c_3));

% denominator coeff
a_11 =  1.0000;
a_12 = -1.1459;
a_13 = 0.9204;
a_21 = 1.0000;
a_22 = -1.3506;
a_23 =  0.9289;
a_31 = 1.0000;
a_32 = -1.2079;
a_33 = 0.8541;

d_1 = [a_11 a_12 a_13];
d_2 = [a_21 a_22 a_23];
d_3 = [a_31 a_32 a_33];

den = conv(d_1,conv(d_2,d_3));

% Group Quantization
v_gb = [num den];
v_gb = quan2N(16,v_gb);
num_qb = v_gb(1:7);
den_qb = v_gb(8:14); 

% individual quantization 

% Giving value of N here
N = 30 ;
%num coeff
qb_11 = quan2N(N,b_11);
qb_12 = quan2N(N,b_12);
qb_13 = quan2N(N,b_13);
qb_21 = quan2N(N,b_21);
qb_22 = quan2N(N,b_22);
qb_23 = quan2N(N,b_23);
qb_31 = quan2N(N,b_31);
qb_32 = quan2N(N,b_32);
qb_33 = quan2N(N,b_33);
%den coeff
qa_11 = quan2N(N,a_11);
qa_12 = quan2N(N,a_12);
qa_13 = quan2N(N,a_13);
qa_21 = quan2N(N,a_21);
qa_22 = quan2N(N,a_22);
qa_23 = quan2N(N,a_23);
qa_31 = quan2N(N,a_31);
qa_32 = quan2N(N,a_32);
qa_33 = quan2N(N,a_33);


C_1 = [qb_11 qb_12 qb_13];
C_2 = [qb_21 qb_22 qb_23];
C_3 = [qb_31 qb_32 qb_33];

numer = conv(C_1,conv(C_2,C_3));

D_1 = [qa_11 qa_12 qa_13];
D_2 = [qa_21 qa_22 qa_23];
D_3 = [qa_31 qa_32 qa_33];

denom = conv(D_1,conv(D_2,D_3));


hold all
freqz(numer,denom);

freqz(num,den);

[H,W] = freqz(numer,denom);
[H1,W1] = freqz(num,den);

figure()
%unquantized filter freq response
plot(W1/(2*pi),(abs(H1)))
%magnitude
 hold on
% %quantized filter freq response
 plot(W/(2*pi),(abs(H)))
% 
 title('magnitude of Frequency response');
 xlabel('Frequency [Hz]'), ylabel('Amplitude [dB]')
axis([0 0.5 0 1]), grid
% 
 


 
