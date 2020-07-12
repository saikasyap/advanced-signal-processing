% Random Matrix Theory Model for DMR Notchdepth
%%%By Sai Kasyap%%%%%
% p--> Nx1 column vector
% D --> Number of narrowband planewave Signals(Interferers)
%N ---> Number of Elements in ULA
%Vm--> Replica Vector in the look direction
% Cosine is found as mentioned in the cox paper
% FIgure 5 in RMT model
clc;
clear all;
%   Detailed explanation goes here
d = 0.5;     
N = 50 ; % Number of sensors
D = [0:1:N-1].';
vm = exp(j*2*pi*d*D*0); % Broad Side
% replica vector 
u=[-1:0.001:1]; 
V = exp(j*2*pi*d*D*u);
% Interferer location
ui = 0.06;
% intialising SCM as zero intially
Sc = zeros(N,N);
S2= zeros(50,50);
%nt = 3000;
% Number of snapshots
L=[2 3 10 20 30 100 200 300 1000 2000 3000 10000];

INRrange = [0 10 20 30 40];
Nd = zeros(length(INRrange),length(L));
WNG = zeros(length(INRrange),length(L));
INR = 10.^(INRrange/10);

ui = 0.06;
vi = exp(-j*2*pi*d*D*ui);

% number of monte carlo trails
nt = 3000;




% number of monte carlo trails
for m=1:length(INR)
    disp(['loop ' int2str(m) ' of 5 ...'])
    for q=1:length(L)
        for k1 = 1:nt

b = sqrt(INR(m)/2)*(randn(1,L(q))+j*randn(1,L(q))); % complex circular gaussian RV
n = sqrt(1/2)*(randn(N,L(q))+j*randn(N,L(q)));
p = vi*b+n;
S  = p*p';  
SCM = S/L(q); % Structured Covariance matrix

% Finding eigen values and eigen vectors

[Sevec1,Seval]=eig(SCM);

% Sorting them in descending order
[Seval,ind] = sort(diag((Seval)),'descend');   % sort eigenvalues in descending order 
  Sevec = Sevec1(:,ind);         % arrange eigenvectors in same order

% finding the estimated noise power..
% Dl - Number of planewaves assumed to be 1.
dl =1;
sn =  (L(q)/L(q)-1)*(1/(N-dl))*(sum(Seval(2:N)));
% Finding s- dmr
S1 = Seval(1,:)*(Sevec(:,1)*Sevec(:,1)');

for i = dl+1:N
S2 = S2+(sn*(Sevec(:,i)*Sevec(:,i)'));
end
Sdmr = S1+S2;

mcosq = gencos(Sevec(:,1),vm);
% % Need to find cosq b/w SCM of ei,vm

gw = (Seval(1,:)-sn)/Seval(1,:);

Wdnum = vm-(gw*Sevec(:,1)*Sevec(:,1)'*vm);
Wdden = vm'*vm*(1-gw*mcosq);
Wdmr = Wdnum/Wdden;
 Ne(k1) = (abs(Wdmr'*vi).^2);


end
 Nd(m,q) = mean(Ne);
end
end


Notchdepth = 10*log10(Nd);


cosq =gencos(vm,vi);
sinq = (1-cosq);

tanq = sinq/cosq;
cotq = cosq/sinq;


for q = 1: length(L)
    
c(q) =N/L(q);
end

%CBF
NDcbf = 10*log10(cosq);
%Ensmeble Notchdepth
NDens = NDcbf-20*log10(N*INR*sinq+1);
for m=1:length(INR)
    disp(['loop ' int2str(m) ' of 5 ...'])
    for q=1:length(L)
       
%Defining Trignometric functions
% Defining Break points
	 c1(m) =(cotq/INR(m)); 
     c2(m) = (INR(m)*tanq);
	c3(m) =(1+N*INR(m)*sinq);
    
    if ((N*INR)<0.1)
	warning('for weak interferer ');
	ND = NDcbf*ones(size(c));
end
    if (INR<1)
	warning('RMT model requires INR>>1');
end
if (N*INR*sinq<1)
	warning('RMT model requires N*INR*sin2>>1');
end
if ((N*INR)^2< c(q))
    warning('RMT model requires (N*INR)^2>c');
end
if (c(q) < c1(m))
    ND(m,q) = NDens(m);
end
   if (c(q) == c1(m))
     ND(m,q) = NDens(m,q) +10*log10(sqrt(2)); 
   end
 if ((c1(m))<c(q))&&(c(q)< c2(m))
     O(q) = 10*log10(c(q));
   r2(m,q) =  NDens(m)+10*log10(c(q))-10*log10(c1(m));
   ND(m,q) = r2(m,q);
 end
 if (c(q) ==c2(m))
     ND(m,q) = r2(m,q)+10*log10(sqrt(2));
 end
 if (c(q)>c2(m)&&c(q)<c3(m))
   r3(m,q) = NDens(m)+20*log10(c(q))-20*log10(c2(m))-10*log(c1(m));
 ND(m,q) = r3(m,q);
%r3(m,q) = r2(m,q)+20*log10(c(q))-20*log10(c2(m));
%ND(m,q) = r2(m,q)+20*log10(c(q))-20*log10(c2(m));
 end
 
   if (c(q)==c3(m))
      ND(m,q) = r3(m,q)-10*log10(sqrt(2));
   end
  if (c(q) > c3(m))
    ND(m,q)= NDcbf;
  end
    end
end

figure
semilogx(L,ND(1,:),'--',L,ND(2,:),'--',L,ND(3,:),'--',L,ND(4,:),'--',L,ND(5,:),'--');
hold on
 semilogx(L,Notchdepth(1,:),'o',L,Notchdepth(2,:),'d',L,Notchdepth(3,:),'s',L,Notchdepth(4,:),'*',L,Notchdepth(5,:),'^');
title('RMT Predictions ND vs Snapshots')
grid
xlabel('L')
ylabel('10log10(ND)')  
ylim([-140 0])
legend('{\it INR} = 0 dB','3{\it INR} = 10 dB','{\it INR} = 20 dB','{\it INR} = 0 dB','{\it INR} = 40 dB')


  
