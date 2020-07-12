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
INRrange = [-40 -25 -20 -17.5 -15 -12.5 -10 0 10 20 30 40];
INR = 10.^(INRrange/10);
L = [6 50 500 5000 50000];
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

for t = 1:length(INR)
if INRrange(t) < 10*log10(1/(N*sinq))

    Ndl(t) = 10*log10(abs(cosq));
  %constant offset
end
% At Break point 
if INRrange(t) == 10*log10(1/(N*sinq))

Ndl(t) = 10*log10(cosq)+10*log10(sqrt(2));
end

    if INRrange(t) >= 10*log10(1/(N*sinq))
   
    Ndl(t) = 10*log10(cosq) + -20*log10((sinq*N*10.^(INRrange(t)/10)));
    end
end

NDens = Ndl;


for m=1:length(INR)
    disp(['loop ' int2str(m) ' of 5 ...'])
    for q=1:length(L)
       
         if (c(q)>=1/sinq^2)
           inrphase = sqrt(c(q))/N;
          
	ND1(INR<=inrphase) = NDcbf;
	ND1(INR>inrphase) = NDcbf-10*log10(INR(INR>inrphase))+10*log10(inrphase);
       
       
     else
           inr1 = 1/(N*sinq);
			inr2(q) = ((1+c(q))).^2/(c(q)*tanq);
            

 
       if (INR(m)<inr1)
           ND(m,q) = NDcbf;
       end
       if (INR(m)==inr1)
           ND(m,q) = NDcbf-10*log(sqrt(2));
       end
       if(INR(m)>inr1&&INR(m)<inr2(q))
           r20(m,q)= NDcbf-20*log10(INR(m))+20*log10(inr1);
           ND(m,q) = NDcbf-20*log10(INR(m))+20*log10(inr1);
       end
       if(INR(m)==inr2(q))
           ND(m,q) =r20(m,q)+10*log(sqrt(2));
       end
       if(INR(m)>inr2(q))
           ND(m,q) = NDcbf-10*log10(INR(m))-10*log10(inr2(q))+20*log10(inr1);
       end
       end
    end   
     
   
end
figure
plot(INRrange,ND1,'--',INRrange,ND(:,2),'--',INRrange,ND(:,3),'--',INRrange,ND(:,4),'--',INRrange,ND(:,5),'--',INRrange,NDens,'--');
ylim([-140 0]);
xlim([-40 40]);
legend('{\it c} = 10','{\it c} = 1','{\it c} = 0.1','{\it c} = 0.01','{\it c} = 0.001','{\it Ensemble}')
 title(' RMT predictions Notchdepth vs INR')
xlabel('10log10(INR)')
  ylabel('10log10(ND)')
       