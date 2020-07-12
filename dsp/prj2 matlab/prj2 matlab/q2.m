%Question 2b 
%For N = 3 bits
%individual quantization vs group Quantization 
%Individual
N=3;
bq = quan2N(N,-1.13);
bq1 = [1 bq];
 aq1 = quan2N(N,-1.5);
 aq2 = quan2N(N,0.9); 
 aq = [1 aq1 aq2];
 % Ploting pole zeros
 figure
 subplot(2,1,1)
 zplane(bq1,aq); title('H(z) Poles and zeros quantized individually')
 [z,p,k] = tf2zp(bq1,aq)

 %For Group quantization 
 
b = [1 -1.13];
a = [1 -1.5 .9];
v = [b a];
v_q = quan2N(N,v);
b_q = v_q(1:2);
a_q = v_q(3:5); 
% Ploting pole zeros
subplot(2,1,2)
 zplane(b_q,a_q); title('H(z) Poles and zeros quantized in group')
  [z1,p1,k1] = tf2zp(b_q,a_q)
