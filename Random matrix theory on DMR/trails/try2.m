clc;
clear all;
d = 0.5;     
N = 50 ; % Number of sensors
D = [0:1:N-1].';
vm = exp(j*2*pi*d*D*0); % Broad Side
% replica vector 
% Interferer location
ui = 0.06;
% intialising SCM as zero intially
vi = exp(-j*2*pi*d*D*ui);
nt = 3000;
% Number of snapshots
INRrange = [0 10 20 30 40];
INR = 10.^(INRrange/10);
L = [2 3 5 20 30 50 200 300 500 2000 3000 5000 20000 30000 50000 1e6];
ui = 0.06;
Vi = exp(-j*2*pi*d*D*ui);


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
ylim([-140 0]);
xlim([0 1e6]);
legend('{\it INR} = 0','{\it INR} = 10','{\it INR} = 20','{\it INR} = 30','{\it INR} = 40')
 title(' RMT predictions Notchdepth vs L')
xlabel('Snapshots')
  ylabel('10log10(ND)')
 
