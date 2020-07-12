% Random Matrix Theory Model for DMR Notchdepth
%%%By Sai Kasyap%%%%%
% p--> Nx1 column vector
% D --> Number of narrowband planewave Signals(Interferers)
%N ---> Number of Elements in ULA
%Vm--> Replica Vector in the look direction
% Cosine is found as mentioned in the cox paper

clc;
clear all;

%Defining Trignometric functions
cosq =gencos(vm,vi);
sinq = (1-cosq);

tanq = sinq/cosq;
cotq = cosq/sinq;





	c1 = cot2/INR;
	c3 = INR*tan2;
	c4 = N*INR*sin2;
	tenlogregion = (c1<=c)&(c<c3);
	twentylogregion = (c3<=c)&(c<c4);

	BPasympt = NaN*ones(size(c));
	BPasympt(c<c1) = BPens;
	BPasympt(tenlogregion) = BPens+10*log10(c(tenlogregion))-10*log10(c1);
	BPasympt(twentylogregion) = BPens+10*log10(c3/c1)+...
		20*log10(c(twentylogregion))-20*log10(c3);
	BPasympt(c4<=c) = BPconv;