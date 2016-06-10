clear;
clear all;
imageName = 'bungee.png';
originImage = imread(imageName);

%如果要重畫mask，把以下的註解取消
%mask = im2uint8(roipoly(img));
mask = imread('bungee-mask.png');
%figure;imshow(img);
%figure;imshow(mask);
[m,n] = size(mask);
boundary = zeros(m,n);

source = originImage;
for i=1:m
    for j=1:n
        if mask(i,j) == 255
            source(i,j,:) = 0;
        end
    end
end
figure;imshow(source);title('source');


% 計算boundary 和boundary的index
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
grayOriginImage = rgb2gray(originImage);
grayOriginImage = im2single(grayOriginImage);
[fx,fy] = gradient(grayOriginImage);
gradtnorm = sqrt((fx .^2) + (fy .^2));
data = gradtnorm + log10(1+gradtnorm);
%figure;imshow(fx);
%figure;imshow(fy);
%figure;imshow(data);title('data');

%patch的大小
patch = 13;

%計算confidence (patch裡面原圖與mask的比例)
confidence = double(1-mask);
%figure;imshow(confidence);

%indexM是boundary總共的index數目
[indexM, indexN] = size(index);

patchIndex = floor(patch/2);

for i=1:indexM %每個boundary 的 index跑一次
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

while ~isequal(mask, zeros(m,n) ) %as research 等到 ? = 0為止
                                  %? 就是這裡的mask
    %找最大的priority
    max = -1;
    %high priority X, Y 座標
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
    highPriority = uint8(zeros(patch, patch, 3));
    for i=1:patch
       for j=1:patch
           highPriority(i,j,:) = source(hpX+i-patchIndex-1, hpY+j-patchIndex-1, :);
       end
    end
    figure;imshow(highPriority);
    
    
    %找最相似的那塊
    target = zeros(patch, patch);
    for i=patchIndex:m-patchIndex
        for j=patchIndex:n-patchIndex
            
        end
    end

end
