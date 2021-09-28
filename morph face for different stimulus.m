%% initialize
clc;
close all
clear;
I1=imread('Capture.jpg');
I2=imread('Capture99.jpg');

num_rep=9;


counter=0;
initiall_stim=[-0.5 -.3 -.14 0 .14 .3 .5];
initial_stim=repmat(initiall_stim,num_rep,1,1);
%% Stimulus input
for index_stim=1:numel(initial_stim)

 clearvars -except index_stim initial_stim counter w I2 I1 out_image
stim=initial_stim(index_stim)
 
w1=(1+stim)/2;
w2=(1-stim)/2;
I_morphed=w1.*I1+w2.*I2;
mu =stim;
sigma = mu.*.2;

stimulus=mvnrnd(mu,sigma^2,3);
w_one=(1+stimulus)/2;
w_two=(1-stimulus)/2;

w_n=[w_one(1) w_two(1)];
w_m=[w_one(2) w_two(2)];
w_re=[w_one(3) w_two(3)];
w_le=[w_one(3) w_two(3)];
w{index_stim}=[w_m;w_n;w_re;w_le]; %M_N_RE_LE




image1=I1; 
image2=I2; 
I1=im2double(I1);
I2=im2double(I2);
I_morphed=im2double(I_morphed);
% imshow(I_morphed)
 %% detect mouth
MouthDetect = vision.CascadeObjectDetector('Mouth','MergeThreshold',50);
pixel_mouth1=step(MouthDetect,image1);
pixel_mouth2=step(MouthDetect,image2);

for i=1:numel(pixel_mouth1)

    S_mouth(i)=ceil(mean([pixel_mouth1(i) pixel_mouth2(i)]));
end

%% detect nose

NoseDetect= vision.CascadeObjectDetector('Nose','MergeThreshold',100);
pixel_nose1=step(NoseDetect,image1);
pixel_nose2=step(NoseDetect,image2);

for i=1:numel(pixel_mouth1)

    S_nose(i)=ceil(mean([pixel_nose1(i) pixel_nose2(i)]));
end
%% detect righteye and detect leftteye

RightEyeDetect= vision.CascadeObjectDetector('RightEye','MergeThreshold',65);
pixel_RE1=step(RightEyeDetect,image1);
pixel_RE2=step(RightEyeDetect,image2);

for i=1:numel(pixel_mouth1)

    S_RE(i)=ceil(mean([pixel_RE1(i) pixel_RE2(i)]));
end

LeftEyeDetect= vision.CascadeObjectDetector('LeftEye','MergeThreshold',100);
pixel_LE1=step(LeftEyeDetect,image1);
pixel_LE2=step(LeftEyeDetect,image2);

for i=1:numel(pixel_mouth1)

    S_LE(i)=ceil(mean([pixel_LE1(i) pixel_LE2(i)]));
end
S=[S_mouth;S_nose;S_RE;S_LE]; %% M_N_RE_LE


for i=1:4

   X(i,1)=S(i,1);
   Y(i,1)=S(i,2);
   X(i,2)=S(i,1)+S(i,3);
   Y(i,2)=S(i,2)+S(i,4);


       
      
end






for i=1:size(I_morphed,1)
  
    for j=1:size(I_morphed,2)
    
        for index=1:4
            
      if ismember(i,[X(index,1):X(index,2)])==1 && ismember(j,[Y(index,1):Y(index,2)])==1  
    I_morphed(i,j,:)=w{index_stim}(index,1).*I1(i,j,:)+w{index_stim}(index,2).*I2(i,j,:);
    
      end
      
        end
    end
end

% imshow(I_morphed)
% figure
image=I_morphed;

k=10;
z=1;
for i=1:2*k+1
    for j=1:2*k+1
        h(i,j)=(1/(2*pi*(z^2)))*(exp(-1*(((i-k-1)^2)+((j-k-1)^2))/(2*(z^2))));
    end
end
out=imfilter(image,h);




    out_image(:,:, index_stim)=out;

    imshow(out_image(:,:, index_stim))
    pause(.1);



end

save('out_image.mat','out_image','w')
for i=1:size(out_image,3)
image=out_image(:,:,i)
 image=imresize(image,[256 256]);
file=sprintf('image_%d.png', i);
imwrite(image,file);
end