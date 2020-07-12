function [m1,m2,Recogn_index] = eigenfacedetector(TestImage, mV, D, u)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
% TestImage -> image which has to be tested
% Me --> Mean vector
%u ---> Eigen FAce vectors
% DataminusMean --> D

% Images are projected to face space by multiplying V and egn face vector

L = size(u,2);
PrjImg = [];
for i = 1 : L
    PrjImg(:,i) = u'*D(:,i);
  
end

temp = TestImage(:,:,1);

[irow icol] = size(temp);
InImage = reshape(temp,irow*icol,1);
Diff = InImage-mV; % Centered test image
ProjectedTestImage = u'*Diff; % Test image feature vector

% Calculating Euclidean distances


for i = 1 : L
    q = PrjImg(:,i);
      Euc_dist(:,i) =norm( ProjectedTestImage - q ) ;
end

[m1, Recogn_index] = min(Euc_dist);

m2 = max(Euc_dist);


% 
% OutputName = strcat(int2str(Recogn_index),'.jpg');


end

