clear;
clear all;


imagename = 'cow3.jpg';  
object = imread(imagename);
bg3 = im2double(imread('bg3.jpg'));
%bg3 = imresize(bg3, [206 309]);
%imwrite(bg3, 'bg3.jpg');
[oi, oj] = size(bg3);
sourceIM = object;
for i=1:oi
    for j=1:oj
        if bg3(i,j) ~= 1
            sourceIM(i,j,:) = 0;
        end
    end
end
ffsize = 1000000;% number of fill-front points
[h,w,z] = size(sourceIM);
figure(1);
imshow(sourceIM);

thresholdIM = object;
level = graythresh(thresholdIM);
boundary = im2bw(thresholdIM,level);

[m n] = size(boundary);

%expand the object boundary twice
pad = padarray(boundary,[1,1],'symmetric','both');
for i = 2:m+1
    for j = 2:n+1
        if((pad(i,j) == 1) && ((pad(i+1,j) == 0) || (pad(i-1,j) == 0) ...
            || (pad(i,j+1) == 0) || (pad(i,j-1) == 0)))
            boundary(i-1,j-1) = 0;
        else
            boundary(i-1,j-1) = boundary(i-1,j-1);
        end
    end
end

figure;
imshow(boundary);%threshold image

pad = padarray(boundary,[1,1],'symmetric','both');
for i = 2:m+1
    for j = 2:n+1
        if((pad(i,j) == 1) && ((pad(i+1,j) == 0) || (pad(i-1,j) == 0) ...
            || (pad(i,j+1) == 0) || (pad(i,j-1) == 0)))
            boundary(i-1,j-1) = 0;
        else
            boundary(i-1,j-1) = boundary(i-1,j-1);
        end
    end
end

figure;
imshow(boundary);%threshold image

%get the fill-front index


pad = padarray(boundary,[1,1],'symmetric','both');

index = zeros(m,n);
temp = 0;%count the number of fill-front points
for i = 2:m+1
    for j = 2:n+1
        if(pad(i,j) == 1 && (pad(i+1,j) == 0 || pad(i-1,j) == 0 || pad(i,j+1) == 0 || pad(i,j-1) == 0 ))
            index(i-1,j-1) = 1;
            temp = temp +1;
        else
            index(i-1,j-1) = 0;
        end
    end
end

firstffsize = temp;
[firstffx,firstffy] = find(index);


%figure(3);
%imshow(index);%fill-front point

%for data
gsourceIM = rgb2gray(sourceIM);%get gray scale image
gradt = im2single(gsourceIM);
[fx,fy] = gradient(gradt);
gradtnorm = sqrt((fx .^2) + (fy .^2));
data = gradtnorm + log10(1+gradtnorm);
%figure(4);
%imshow(data);

%calaulate  confidence for all pixels
patch = 11;%the patch size
pad = padarray(boundary,[floor(patch/2),floor(patch/2)],'symmetric','both');
confidence = ones(m,n);
for i = 1+floor(patch/2):m+floor(patch/2)
    for j = 1+floor(patch/2):n+floor(patch/2)
        confidence(i-floor(patch/2),j-floor(patch/2)) = (sum(sum(pad(i-floor(patch/2):i+floor(patch/2),j-floor(patch/2):j+floor(patch/2)))))/(patch^2);
    end
end
%figure(5);
%imshow(confidence);

%transfer image from rgb to lab 
colorTransform = makecform('srgb2lab');
sourceIM = applycform(sourceIM , colorTransform);

% start inpainting
temp = 1;

while(temp>0)
%STEP1:get the fill-front index
[ffx,ffy] = find(index); %save fill-front index
ffsize = temp;

%STEP2:get all points priority
priority = confidence .* data; 
%figure(6);
%imshow(priority);

% STEP3:find the highest priority of fill front point (curX,curY)
if(ffsize~=0)
    max = 0; 
    for i=1:ffsize;
        if(priority(ffx(i),ffy(i))>=max)
            curX = ffx(i);
            curY = ffy(i);
            max = priority(ffx(i),ffy(i));
        end
    end
    ffsize = ffsize - 1;%at least one pixel will be inpaint
end

%STEP4:find the best match patch
target = zeros(patch,patch,3);
match = zeros(patch,patch,3);
best = zeros(patch,patch,3);
for u=1:patch;
    for v =1:patch;
        target(u,v,:) = sourceIM(curX-floor(patch/2)+u+1,curY-floor(patch/2)+v-1,:);
    end
end
%find min SSD
min = 100000000000000;
sX = 1;
sY = 1;
i = 1;
while(i<(m-1-patch))
    j = 1; 
    while( j < (n-1-patch))
        out = 0;
        for u=1:patch;
            for v = 1:patch;
                if(boundary(i+u-1,j+v-1)==0 || object(i+u-1,j+v-1)~=255)
                    out = 1;
                end
                match(u,v,:) = sourceIM(i+u-1,j+v-1,:);
            end
        end
        % norm
        if(out~=1)
            tempScore = (target - match) .^2;
            ssd = 0;
            for u = 1:patch;
                for v = 1:patch;
                    if(boundary(curX-floor(patch/2)+u-1,curY-floor(patch/2)+v-1)~=0)
                        ssd = ssd + sqrt(tempScore(u,v,1)^2 + tempScore(u,v,2)^2 + tempScore(u,v,3)^2);
                    end
                end
            end
            corrScore=ssd/(patch*patch);
            %corrScore = sqrt(corrScore(1)^2 + corrScore(2)^2 + corrScore(3)^2 );
            if(corrScore <= min)
                min = corrScore;
                best = match;
                sX = i;
                sY = j;
            end
        end
        j = j + patch;
    end
    i = i + patch;
end
disp(temp);
%STEP5:copy this patch to target
for u=0:patch-1;
    for v =0:patch-1;
        if(boundary(curX-floor(patch/2)+u,curY-floor(patch/2)+v)==0 || index(curX-floor(patch/2)+u,curY-floor(patch/2)+v)==1 || object(curX-floor(patch/2)+u,curY-floor(patch/2)+v)~=255 )
            sourceIM(curX-floor(patch/2)+u,curY-floor(patch/2)+v,:) = best(u+1,v+1,:);
        end
    end
end
% update boundary and pad and object
lastBD = boundary;
for i = curX-floor(patch/2):curX+floor(patch/2);
    for j = curY-floor(patch/2):curY+floor(patch/2);
        object(i,j) = 255;
        boundary(i,j) = 1;
    end
end
pad = padarray(boundary,[1,1],'symmetric','both');
temp = 0;%count the number of fill-front points
for i = 2:m+1
    for j = 2:n+1
        if(pad(i,j) == 1 && (pad(i+1,j) == 0 || pad(i-1,j) == 0 || ...
            pad(i,j+1) == 0 || pad(i,j-1) == 0 ))
            index(i-1,j-1) = 1;
            temp = temp +1;
        else
            index(i-1,j-1) = 0;
        end
    end
end

%STEP6: update confidence
for i = curX-floor(patch/2):curX+floor(patch/2);
    for j = curY-floor(patch/2):curY+floor(patch/2);
        if(lastBD(i,j)==0)
            confidence(i,j) = (sum(sum(index(i-floor(patch/2):i+floor(patch/2),j-floor(patch/2):j+floor(patch/2)))))/(patch^2);
        end
    end
end


end
% inpainting end

%smooth the original boundary

%transfer image from lab to rgb 

colorTransform = makecform('lab2srgb');
sourceIM = applycform(sourceIM , colorTransform);
figure(7);
imshow(sourceIM);
%imwrite(sourceIM,'C:\Users\steven\Desktop\final3.jpg','jpg')



