clear;
clear all;
%imageName = 'bungee.png';
imageName = 'air';
originImage = imread([imageName, '.jpg']);

%如果要重畫mask，把以下的註解取消
%if you want to redraw mask, cancel the comment below
mask = im2uint8(roipoly(originImage));
imwrite(mask, [imageName, '-mask.jpg']);
%mask = imread('bungee-mask.png');


%figure;imshow(img);
%figure;imshow(mask);
[m,n] = size(mask);
boundary = zeros(m,n);

%patch的大小
%size of patch
patch = 9;

source = originImage;
for i=1:m
    for j=1:n
        if mask(i,j) == 255
            source(i,j,:) = 0;
        end
    end
end
figure;imshow(source);title('source');


% 計算boundary
% compute the bundary
for i=2:m-1
    for j=2:n-1
        
        for k=-1:0
            for l=-1:0
                if mask(i+k,j+l) ~= mask(i,j)
                   boundary(i,j) = 1;
                end
            end
        end
    end
end

%index 為mask內所有點座標
%index is all coordinate(x,u) in mask
index = [];
for i=1:m
    for j=1:n
        if mask(i,j) == 255
            index = [index; i j];
        end
    end
end
%figure;imshow(boundary);title('boundary');

%計算data(直接套用gradient)
%compute "data"(use matlab function "gradient")
grayOriginImage = rgb2gray(originImage);
grayOriginImage = im2single(grayOriginImage);
[fx,fy] = gradient(grayOriginImage);
gradtnorm = sqrt((fx .^2) + (fy .^2));
data = gradtnorm + log10(1+gradtnorm);
%figure;imshow(fx);
%figure;imshow(fy);
%figure;imshow(data);title('data');



%計算confidence (patch裡面原圖與mask的比例)
% compute "confidence 
confidence = double(1-mask);
%figure;imshow(confidence);

%indexM是indexy總數目
%indexM is the length of index

[indexM, indexN] = size(index);

patchIndex = floor(patch/2);

for i=1:indexM %每個index跑一次 %run every index
   count = 0;
   for px=-patchIndex:patchIndex
       for py=-patchIndex:patchIndex
           if mask(index(i,1)+px, index(i,2)+py) == 0
               count = count + 1;
           end
       end
   end
   confidence(index(i,1), index(i,2)) = double(count/(patch*patch));
end

priority = data.* confidence;
[priX, priY] = size(priority);
figure;imshow(priority);title('priority');

%開始inpainting了
%START TO INPAINTING

while ~isequal(mask, zeros(m,n) ) %as research 等到 omega = 0為止
                                  %omega 就是這裡的mask
                                  %as sesearch, until omdga = 0
                                  %omega is mask here
    %找最大的priority
    %find maximum priority
    max = -1;
    
    %high priority X, Y 座標
    %high priority coordinate (hpX, hpY)
    hpX = 0;
    hpY = 0;
    for i=1:indexM
        if priority(index(i,1), index(i,2) ) > max
            max = priority(index(i,1),index(i,2));
            hpX = index(i,1);
            hpY = index(i,2);
        end
    end

    %決定highPriority 是哪塊
    %find highPriority
    highPriority = uint8(zeros(patch, patch, 3));
    for i=1:patch
       for j=1:patch
           highPriority(i,j,:) = source(hpX+i-patchIndex-1, hpY+j-patchIndex-1, :);
       end
    end
    figure;imshow(highPriority);title('highPriority');
    
    
    %找最相似的那塊
    %find the most similar target(patch X patch)
    
    
    %避免target 碰到原圖邊界
    %to avoid "target" hit the border of originImage
    for i=patchIndex+1:m-patchIndex-1
        for j=patchIndex+1:n-patchIndex-1
            
            target = uint8(zeros(patch, patch, 3));
            %fill the target
            for px=-patchIndex:patchIndex
                for py=-patchIndex:patchIndex
                    target(px+patchIndex+1, py+patchIndex+1, :) = originImage(i+px, j+py,:);
                    %target(1~patch)
                end
            end
            figure;imshow(target);title('target');
            
            
        end
    end

end

