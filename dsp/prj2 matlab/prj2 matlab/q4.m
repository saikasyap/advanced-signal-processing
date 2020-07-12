%Question 4
[b_b a_b] = butter(3,[0.25 0.30],'bandpass');
% question a i
b_b
a_b
% question a ii
b_r = roots(b_b);

%roots
b1 = b_r(1);
b2 = b_r(2);
b3 = b_r(3);
b4 = b_r(4);
b5 = b_r(5);
b6 = b_r(6);

% iii
q1 =poly([b1 b2]);
q2= poly([b3 b4]);
q3 = poly([b5 b6]); 



b_11 = q1(1);
b_12 = q1(2);
b_13 = q1(3);
b_21 = q2(1);
b_22 = q2(2);
b_23 = q2(3);
b_31 = q3(1);
b_32 = q3(2);
b_33 = q3(3);
%  Question iv 
a_r = roots(a_b);
%roots
a1 = a_r(1);
a2 = a_r(2);
a3 = a_r(3);
a4 = a_r(4);
a5 = a_r(5);
a6 = a_r(6);

p1 = poly([a1 a2]);
p2 = poly([a3 a4]);
p3 =poly([a5 a6]);


a_11 =  p1(1);
a_12 = p1(2);
a_13 = p1(3);
a_21 = p2(1);
a_22 = p2(2);
a_23 =  p2(3);
a_31 = p3(1);
a_32 = p3(2);
a_33 = p3(3);

c_1 = [b_11 b_12 b_13];
c_2 = [b_21 b_22 b_23];
c_3 = [b_31 b_32 b_33];

num = conv(c_1,conv(c_2,c_3));
% denominator coeff


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
N = 3;
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

[H,W] = freqz(numer,denom);
[H1,W1] = freqz(num,den);

figure()
%unquantized filter freq response
plot(W1/(2*pi),(abs(H1)),'-.r')
%magnitude
 hold on
% %quantized filter freq response
 plot(W/(2*pi),(abs(H)))
% 
 title('magnitude of Frequency response');
 xlabel('Frequency [Hz]'), ylabel('Amplitude ')
axis([0 0.5 0 5000]), grid
legend('unquantized','quantized');




