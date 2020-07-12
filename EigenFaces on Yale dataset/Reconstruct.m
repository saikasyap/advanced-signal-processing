function [ output_args ] = Reconstruct(u,L1,K,Me )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% u -->Eigen face vector
% L1 --> Input Image
% K ---> Choosing K best Eigen Values
% 
% Me --> Mean vector
% 



[irow,icol,b] = size(L1);

B = reshape(L1,irow*icol,b);
figure ;
subplot(1,2,1)
imagesc(L1(:,:,1));  colormap('gray');
%[mV2,V2,u2 ] = lowrankrep(10,L1);
  %m is the mean image, u is the eigenvector
  title('orignal Image')
  projection = u(:,1:K)'*(B-Me);
ReshapedImage = u(:,1:K)*projection+Me;

%ReshapedImage = u*V(:,1:K);

ReshapedImage = reshape(ReshapedImage,irow,icol);

%show the reconstructed image.
subplot(1,2,2)
imagesc(ReshapedImage(:,:,1)); colormap('gray');

title('Image after Reconstruction');

end

