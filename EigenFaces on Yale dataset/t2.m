% function to compute low rank representation of image data

% Converting the 3 d vector in to 2 - d vector

clc;
% load('projldata.m');
desiredposes = fieldnames(facedata);

[I,L] = storedata((1:15),desiredposes,facedata);
% I is matrix is of the form 195 x 161 x 165

[m,n,numIMG] = size(I);

A = reshape(I,[m*n numIMG]);

% Removing Mean prior to computing eigen decomposition of the sample
% covariance  matri

[nrow,ncol] = size(A);

D = A-(mean(A,2)*ones(1,ncol));
Me = mean(A,2);
% Finding mean image
Mimg =  reshape(Me,m,n);

% Sample covariance matrix
C = D'*D;
% calculating svd
[U,S,V] = svd(C,0);
% U - eigen vector of A'*A;
%S - eigen values 


[S,ind] = sort(diag(abs(S)),'descend');   % sort eigenvalues in descending order 
   Vs     = V (:,ind);          % arrange eigenvectors in same order

   %normalising eigen vectors
for i=1:size(V,2) %access each column
kk=Vs(:,i);
temp=sqrt(sum(kk.^2));
Vs(:,i)=Vs(:,i)./temp;

end

   
   u = D*V;
   
   % normalising
  
% calculation of wieght vectors


for i = 1 :size(u,2)
    kk=u(:,i);
temp=sqrt(sum(kk.^2));
u(:,i)=u(:,i)./temp; 
    
end 




% selecting best K Eigen vectors
% K = 40;
% Dimensionality reduction.  
% fprintf('Creating lower dimensional subspace\n') 
% u1 = u(:, 1:K);
% D1 = D(:, 1:K);

% Caculating Eigen Vectors

%finding min error through eculdian distance


% display the eigenvalues for C
for j = length(S);
    S1(j) = S(j)^2;
end

normalised_evalues = S / sum(S);
figure(1)
 plot(cumsum(normalised_evalues));
title('choosing best K eigen vectors')
xlabel('No. of eigen values'), ylabel('Variance accounted for');
xlim([1 170]), ylim([0 1.05]), grid on;


% 
%plot eigen vectors
num_eigenfaces = 15;
%display the eigenvectors
figure;
figure(2);
K = 38;
u = u(:,1:K);
D1 = D(:, 1:K);
Vs1 = Vs(:,1:K);
%Caculating Wieght vector

w = u'*D1;

for i=1:size(u,2)
img=reshape(u(:,i),m,n);
img=histeq(img,255);
 subplot(2, ceil(K/2), i);
imshow(img)
drawnow;

if i==3
title('Eigenfaces','fontsize',18)
end
end

%for a single image taken or tested

%Reconstruction for a test image
[I1,L1] = storedata(15,{'normal'},facedata);

[irow,icol,b] = size(L1);

B = reshape(L1,irow*icol,b);
figure (3)
subplot(1,2,1)
imagesc(L1(:,:,1));   colormap('gray');
[mV2,V2,u2 ] = lowrankrep(10,L1);
title('orignal Image')
 % m is the mean image, u is the eigenvector
  
  projection = u(:,1:K)'*(B-Me);
ReshapedImage = u(:,1:K)*projection+Me;


ReshapedImage = reshape(ReshapedImage,irow,icol);
subplot(1,2,2)
%show the reconstructed image.
imagesc(ReshapedImage(:,:,1)); colormap('gray');

title('Image after Reconstruction');






%Storing all the distances of Training data set


for i = 1:numIMG
[m1(i) ,m2(i) r(i) ] = eigenfacedetector(I(:,:,i), Me, D, u);
end





















