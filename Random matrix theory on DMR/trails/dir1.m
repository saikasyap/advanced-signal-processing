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
L=[2 3 10 20 30 100 200 300 1000 2000 3000 10000 20000 30000 50000];

INRrange = [0 10 20 30 40];
Nd = zeros(length(INRrange),length(L));
WNG = zeros(length(INRrange),length(L));
INR = 10.^(INRrange/10);

ui = 0.06;
vi = exp(-j*2*pi*d*D*ui);

% number of monte carlo trails
nt = 3000;




% % number of monte carlo trails
% for m=1:length(INR)
%     disp(['loop ' int2str(m) ' of 5 ...'])
%     for q=1:length(L)
%         for k1 = 1:nt
% 
% b = sqrt(INR(m)/2)*(randn(1,L(q))+j*randn(1,L(q))); % complex circular gaussian RV
% n = sqrt(1/2)*(randn(N,L(q))+j*randn(N,L(q)));
% p = vi*b+n;
% S  = p*p';  
% SCM = S/L(q); % Structured Covariance matrix
% 
% % Finding eigen values and eigen vectors
% 
% [Sevec1,Seval]=eig(SCM);
% 
% % Sorting them in descending order
% [Seval,ind] = sort(diag((Seval)),'descend');   % sort eigenvalues in descending order 
%   Sevec = Sevec1(:,ind);         % arrange eigenvectors in same order
% 
% % finding the estimated noise power..
% % Dl - Number of planewaves assumed to be 1.
% dl =1;
% sn =  (L(q)/L(q)-1)*(1/(N-dl))*(sum(Seval(2:N)));
% % Finding s- dmr
% S1 = Seval(1,:)*(Sevec(:,1)*Sevec(:,1)');
% 
% for i = dl+1:N
% S2 = S2+(sn*(Sevec(:,i)*Sevec(:,i)'));
% end
% Sdmr = S1+S2;
% 
% mcosq = gencos(Sevec(:,1),vm);
% % % Need to find cosq b/w SCM of ei,vm
% 
% gw = (Seval(1,:)-sn)/Seval(1,:);
% 
% Wdnum = vm-(gw*Sevec(:,1)*Sevec(:,1)'*vm);
% Wdden = vm'*vm*(1-gw*mcosq);
% Wdmr = Wdnum/Wdden;
%  Ne(k1) = (abs(Wdmr'*vi)^2);
% 
% 
%         end
%       
%       Nd(m,q) = mean(Ne);
% end
% end
% 
% 
% Notchdepth = 10*log10(Nd);

% RMT PRedications .. Start From here
% Notchdepths vs Snapshots

cosq =gencos(vm,vi);
maxL=1e6;

for q = 1 : length(L)
    c = N/L(q)
end
for ind=1:length(INR)
  bpasymptc(:,ind)=DMR_asympt(c,INR(ind),N,vm,vi);
end


figure 
hold on

% semilogx(L,Notchdepth(1,:),'r',L,Notchdepth(2,:),'g',L,Notchdepth(3,:),'b',L,Notchdepth(4,:),'k',L,Notchdepth(5,:),'m');
 
%semilogx(L,Notchdepth(1,:),'-o',L,Notchdepth(2,:),'-d');
semilogx(L,bpasymptc(:,1),'o',L,bpasymptc(:,2),'-d',L,bpasymptc(:,3),'s',L,bpasymptc(:,4),'-*',L,bpasymptc(:,5),'^');

title('Comparison of model predications ND vs  L ')

grid

xlabel('L')
  ylabel('10log10(ND)')
 
 ylim([-140 0])

