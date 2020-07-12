clc;
clear all;
d = 0.5; 
INRrange = [-40 -25 -20 -17.5 -15 -12.5 -10 0 10 20 30 40];
N = 50 ; % Number of sensors
D = [0:1:N-1].';
vv = exp(j*2*pi*d*D*0); % Broad Side
% replica vector 
u=[-1:0.001:1]; 
V = exp(j*2*pi*d*D*u);
% Interferer location
ui = 3/N;
Vi = exp(-j*2*pi*d*D*ui);
% sensor spacing half wavelength wrt wc
% INR = 40 db => Si = 100, Sw =1;
%ECM for recieved signal is

 k1= 1;
for INR = 10.^(INRrange/10)
   
 Sx = INR*(Vi*Vi') + eye(N);
dl =1 ; % Number of Strong planewavesignals (interferers);
% Finding Eigen vectors and eigen values
[evec1,egval] = eig(Sx);
% Sorting them in descending order
[egval,ind] = sort(diag((egval)),'descend'); % sort eigenvalues in descending order 
evec = evec1(:,ind); % arrange eigenvectors in same order
% Finding Sdmr
S1 = egval(1,:)*(evec(:,1)*evec(:,1)');
sn = (1/(N-dl))*(sum(egval(2:N)));
S2= zeros(50,50);
for i = dl+1:N
S2 = S2+(sn*(evec(:,i)*evec(:,i)'));
end
Sdmr = S1+S2;
% Weight vector of dmr
Wd1 = inv(vv'*(inv(Sdmr))*vv)*((inv(Sdmr))*vv);  
ND(k1) =  abs(Wd1'*Vi)^2;
Nde(k1) = 10*log10(ND(k1));
WNGen(k1) = 1 /(Wd1'*Wd1);

k1 = k1+1;

end

% Weight vector of conventional
Wc = ones(N,1)/N;
Ndcbf = 10*log10(abs(Wc'*Vi)^2);



[Ndl1] = NDINR(INRrange,5,Vi);
Ndl2 =NDINR(INRrange,50,Vi);
 Ndl3   =NDINR(INRrange,500,Vi);
  Ndl4   =NDINR(INRrange,5000,Vi);


figure 
 plot((INRrange),Nde,'x');
 grid on
  hold on
 plot((INRrange),Ndcbf*ones(size(INRrange)),'d');
  plot((INRrange),Ndl1,'o');
% 
  plot((INRrange),Ndl2,'*');
 plot((INRrange),Ndl3,'s');
  plot((INRrange),Ndl4,'+');

hold off
title(' DMR as a function of INR for the cannonical single interferer example')
xlabel('10log10(INR)')
  ylabel('10log10(NDens)')
 xlim([-40 40])
 ylim([-140 0])
 legend('{\it c} = 10','{\it c} = 1','{\it c} = 0.1','{\it c} = 0.01','{\it Ensemble}')