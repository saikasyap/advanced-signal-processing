
clc;

for i = 1 :13
    imagesc(supplementaldata(:,:,i));
truesize;
        colormap(gray);
   
        
        pause(0.2)                      % short pause between plots
end
    


for i = 1:numIMG
[m1(i) ,m2(i) r(i) ] = eigenfacedetector(I(:,:,i), Me, D, u);
end

for j = 1: 13
    [m3(j) ,m4(j) r2(j) ] = eigenfacedetector(supplementaldata(:,:,j), Me, D, u);
end


threshold1 = 1.5199e+04;
threshold2 = 15488.9510553266;

%Reconginition Procedure
for i = 1:13
    
    
    if m3(i)<threshold1
    fprintf('the image is present in database');
    ig = i % to know which images are wrong detections
    Reconstruct(u,supplementaldata(:,:,i),K,Me )
    end
    if threshold2<m3(i)<threshold1
       imageno= i
         fprintf('the image is not  present in database ');
    end
    if m3(i)> threshold2
        imgno = i
         fprintf('the image is not  present in database and not a face image');
         end
end
