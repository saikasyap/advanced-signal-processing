[b_b a_b] = ellip(3,1,50,[0.25 0.30],'bandpass');
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
[H3,W3] = freqz(num_qb,den_qb);

% %quantized filter freq response
 plot(W3/(2*pi),(abs(H3)))
 xlabel('Frequency [Hz]'), ylabel('Amplitude ')
axis([0 0.5 0 500]), grid