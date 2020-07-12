% Figure 5 in the paper
% Comparison of model predictions(dashed lines) and Monte Carlo Trails
% for mean Notch depth vs INR as the sensor ratio varies from c= 10 to 0.01

% Random Matrix Theory Model for DMR Notchdepth
%%%By Sai Kasyap%%%%%
% p--> Nx1 column vector
% D --> Number of narrowband planewave Signals(Interferers)
%N ---> Number of Elements in ULA
%Vm--> Replica Vector in the look direction
% Cosine is found as mentioned in the cox paper

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
L=[5 50 500 5000 50000];

INRrange = [-40:10:40];
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
 Ne(k1) = (abs(Wdmr'*vi)^2);


end
 Nd(m,q) = mean(Ne);
end
end


Notchdepth = 10*log10(Nd);

% RMT PRedications .. Start From here
% Notchdepths vs Snapshots

cosq =gencos(vm,vi);
sinq = sqrt(1-cosq);

% Eigen vector 1
E1 = vi/sqrt(N);

% Orthogonal residual

Et = (vm -(E1'*vm)*E1)/(norm(vm -(E1'*vm)*E1));

%Condition for RMT Results\

% Calculate eigen vectors from the SCM

%vm = sqrt(N)*cos(vm,vi)*E1+sqrt(N)*sin(vm,vi)*Et;

%defining the projections as as p1 and p2

p1 =abs(Sevec(:,1)'*E1).^2;

p2 = abs(Sevec(:,1)'*Et).^2;

% Notch Depth RMT Predictions


NDcbf = cosq;

% Need to define the phase of E1*Sevec(:,1)*Sevec(:,1)'*Et say ph;


ang =  E1'*(Sevec(:,1)*Sevec(:,1)')*Et;
ph = angle(ang);

tanq = sinq/sqrt(cosq);
cotq = sqrt(cosq)/sinq;

% Constant Aspect Ratio
for m=1:length(INR)
    disp(['loop ' int2str(m) ' of 5 ...'])
for q = 1:length(L)
c(q) = N/L(q);

A1(m,q) = (1-ph*tanq*sqrt(INR(m)*sqrt(c(q))));
A2(m,q) = (1-ph*cotq*((INR(m))^-0.5)*sqrt(c(q)));
A3(m,q) = 1 + N*INR(m)*(sinq)^2+c(q);

% Need to plot Notch depth vs INR

% Plot for Notch depth Vs c;

NDr(m,q) = (NDcbf)*abs(A1(m,q)*A2(m,q)).^2/abs(A3(m,q)).^2;
end
end

NDrmt = 10*log10(NDr);

k1 =1;
for INR = 10.^(INRrange/10)
   

    S5 = INR*(vi*vi') + eye(N);
dl =1 ; % Number of Strong planewavesignals (interferers);
% Finding Eigen vectors and eigen values
[evec1,egval] = eig(S5);
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
Wd1 = inv(vm'*(inv(Sdmr))*vm)*((inv(Sdmr))*vm);  
ND(k1) =  abs(Wd1'*vi)^2;
NdENS(k1) = 10*log10(ND(k1));
k1 = k1+1;
end




figure 
hold on

%semilogx(L,Notchdepth(1,:),'--',L,Notchdepth(2,:),'--');
% Exact values
 plot(INRrange,Notchdepth(:,1),'--',INRrange,Notchdepth(:,2),'--',INRrange,Notchdepth(:,3),'--',INRrange,Notchdepth(:,4),'--',INRrange,Notchdepth(:,5),'--',INRrange,NdENS,'-x');
% Predictions
 plot(INRrange,NDrmt(:,1),'o',INRrange,NDrmt(:,2),'d',INRrange,NDrmt(:,3),'s',INRrange,NDrmt(:,4),'*',INRrange,NDrmt(:,5),'^');

%semilogx(L,NDrmt);


title('RMT Predictions ND vs Snapshots')

grid

xlabel('INR range')
  ylabel('10log10(ND)')
 
 ylim([-140 0])

legend('{\it c} = 10','3{\it c} = 1','{\it c} = 0.01','{\it c} = 0.001','{\it Ensemble}')

