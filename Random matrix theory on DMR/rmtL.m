% Random Matrix theory model

clc;
clear all;
%   Detailed explanation goes here
d = 0.5;     
N = 50 ; % Number of sensors
D = [0:1:N-1].';
vm = exp(j*2*pi*d*D*0); % Broad Side
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
INR = 10.^(INRrange/10);
ui = 0.06;
vi = exp(-j*2*pi*d*D*ui);


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
ph1 = angle(ang);
ph = exp(j*ph1);
tanq = sinq/sqrt(cosq);
cotq = sqrt(cosq)/sinq;

% Constant Aspect Ratio
for m=1:length(INR)
    disp(['loop ' int2str(m) ' of 5 ...'])
for q = 1:length(L)
c(q) = N/L(q);

A1(m,q) = (1-ph*tanq*sqrt(INR(m)*sqrt(c(q))));
A2(m,q) = (1-ph*cotq*(inv(sqrt(INR(m))))*sqrt(c(q)));
A3(m,q) = 1 + N*INR(m)*(sinq)^2+c(q);

% Need to plot Notch depth vs INR

% Plot for Notch depth Vs c;

NDr(m,q) = (NDcbf)*abs(A1(m,q)*A2(m,q)).^2/abs(A3(m,q)).^2;
end
end

NDrmt = 10*log10(NDr);
