function [ Data,data] = storedata(idface ,desiredposes,facedata)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Data for all images in training set
%data for set of images/image for a particular individual.

nt = length(facedata);

if nargin<2
    desiredposes=fieldnames(facedata);
end
nposes=length(desiredposes);
i = 1;
for faceid = 1:nt
    for poseid = 1:nposes
        X = facedata(faceid).(desiredposes{poseid});
        [nrow,ncol]=size(X);
    Data(:,:, i) = reshape(X,nrow,ncol);
    i = i+1;
    end
end
%to plot a particular face of an individual
data = Data(:,:,idface);

end

    
%  figure;                 % plot 
%      imagesc(Data(:,:,1));
%         truesize;
%         colormap(gray);
