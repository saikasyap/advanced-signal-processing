function [mV,V,u] = lowrankrep(K,I)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Input
% Select K best eigen vectors depending on variance
% 
% I is a image matrix which is in 2-d vector after using store data
%ouput
%mV---> mean vector
%V---> Eigen vectors
%u1 --> set of eigen vectors defining face class

[m,n,numIMG] = size(I);

A = reshape(I,[m*n numIMG]);

% Removing Mean prior to computing eigen decomposition of the sample
% covariance  matri

[nrow,ncol] = size(A);
%mean vector mV
mV =mean(A,2)*ones(1,ncol); 
D = A-mV;

% Sample covariance matrix
C = D'*D;
% calculating svd
[U,S,V] = svd(C,0);
% U - eigen vector of A'*A;
%S - eigen values 
[S,ind] = sort(diag(abs(S)),'descend');   % sort eigenvalues in descending order 
   V      = V (:,ind);          % arrange eigenvectors in same order
   
 for i=1:size(V,2) %access each column
kk=V(:,i);
temp=sqrt(sum(kk.^2));
V(:,i)=V(:,i)./temp;
 end

   u = D*V;
%    normalising
%    u = u/norm(u);

%ncol  - number of images
%D - difference of Image and mean

for i = 1 : ncol
     kk=u(:,i);
temp=sqrt(sum(kk.^2));
u(:,i)=u(:,i)./temp;  
end  
% U is considered for projection not V

end

