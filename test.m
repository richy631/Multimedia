clear;
clear all;
%imageName = 'bungee.png';
imageName = 'air';
originImage = imread([imageName, '.jpg']);

%如果要重畫mask，把以下的註解取消
%if you want to redraw mask, cancel the comment below
%mask = im2uint8(roipoly(originImage));
%imwrite(mask, [imageName, '-mask.png']);
mask = imread([imageName, '-mask.png']);


%figure;imshow(img);
%figure;imshow(mask);title('mask');
[m,n] = size(mask);
boundary = zeros(m,n);

%patch的大小
%size of patch
patch = 9;

maskedImage = originImage;
for i=1:m
    for j=1:n
        if mask(i,j) == 255
            maskedImage(i,j,:) = 0;
        end
    end
end
%figure;imshow(maskedImage);title('maskedImage');


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
%figure;imshow(boundary);title('boundary');




%計算data(直接套用gradient)
%compute "data"(use matlab function "gradient")
grayMaskedImage = rgb2gray(maskedImage);
grayMaskedImage = im2single(grayMaskedImage);
[fx,fy] = gradient(grayMaskedImage);
gradtnorm = sqrt((fx .^2) + (fy .^2));
data = gradtnorm + log10(1+gradtnorm);
%figure;imshow(fx);
%figure;imshow(fy);
%figure;imshow(data);title('data');


%change RGB to LAB
colorTransform = makecform('srgb2lab');
maskedImage = applycform(maskedImage , colorTransform);



%計算confidence (patch裡面原圖與mask的比例)
% compute "confidence 
confidence = double(1-mask);
%figure;imshow(confidence);title('confidence');



patchIndex = floor(patch/2);

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

%indexM是indexy總數目
%indexM is the length of index

[indexM, indexN] = size(index);

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


%figure;imshow(priority);title('priority');

%開始inpainting了
%START TO INPAINTING

while ~isequal(mask, zeros(m,n) ) %as research 等到 omega = 0為止
                                  %omega 就是這裡的mask
                                  %as sesearch, until omdga = 0
                                  %omega is mask here
	
    %STEP1 find all index in mask
    
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

    %indexM是indexy總數目
    %indexM is the length of index

    [indexM, indexN] = size(index);
    
    %STEP2 compute the priority
    
    priority = data.* confidence;
    [priX, priY] = size(priority);
    
    %STEP3 find the patch with maximum priority
    
    %找最大的priority的patch
    %everytime in the while loop will find a highPriority
    %find the patch with maximum priority
    max = 0;
    
    %high priority X, Y 座標
    %high priority coordinate (hpX, hpY)
    hpX = 0;
    hpY = 0;
    for i=1:indexM
        if priority(index(i,1), index(i,2) ) >= max
            max = priority(index(i,1),index(i,2));
            hpX = index(i,1);
            hpY = index(i,2);
        end
    end

    %決定highPriority 是哪塊
    %find highPriority
    highPriority = uint8(zeros(patch, patch, 3));
    for i= -patchIndex:patchIndex
       for j= -patchIndex:patchIndex
           highPriority(i+patchIndex+1,j+patchIndex+1,:) = maskedImage(hpX+i, hpY+j, :);
       end
    end
    %figure(8);imshow(highPriority);title('highPriority');
    
    %STEP4 find the minimum exemplar
    
    %找最相似的那塊
    %find the most similar target(patch X patch)
    
    
    %避免 碰到原圖邊界
    min = 99999999;
    %to avoid "target" hit the border of originImage
    targetX = 1;
    targetY = 1;
    target = uint8(zeros(patch, patch, 3));
    
    %(tx,ty) is the center of patch*patch(81)
    %so (tx,ty) should in the patchIndex(4) to m-patchIndex(4) 
    for tx=hpX-floor((patch*patch)/2):hpX+floor((patch*patch)/2)
        %out of boundary
        if tx <= patchIndex || tx >= m-patchIndex
            continue;
        end
        for ty=hpY-floor((patch*patch)/2):hpY+floor((patch*patch)/2)
            %out of boundary
            if ty <= patchIndex || ty >= n-patchIndex
                continue;
            end

            %in mask
            %inside (tx, ty) patch
            %there can't be any mask == white inside
            
            
            ismask = 0;
            tmptarget = uint8(zeros(patch, patch, 3));
            %fill the target
            for px=-patchIndex:patchIndex
                for py=-patchIndex:patchIndex
                    if mask(tx+px, ty+py) == 255
                        ismask = 1;
                        break;
                    else
                        tmptarget(px+patchIndex+1, py+patchIndex+1, :) = maskedImage(tx+px, ty+py,:);
                    end
                end
                if ismask == 1
                    break;
                end
            end
            
            if ismask == 1
                continue;
            end
            %figure(7);imshow(tmptarget);title('tmptarget');
            %figure(8);imshow(highPriority);title('highPriority');
            
            
            %compute SSD
            SSD = 0;
            for x=1:patch
                for y=1:patch
                    if mask(hpX-patchIndex-1 + x, hpY-patchIndex-1 + y) == 0
                        L = highPriority(x,y,1) - tmptarget(x,y,1); 
                        a = highPriority(x,y,2) - tmptarget(x,y,2);
                        b = highPriority(x,y,3) - tmptarget(x,y,3);
                        SSD = SSD + L*L + a*a + b*b;
                    end
                end
            end
            
            if SSD < min
               min = SSD;
               targetX = tx;
               targetY = ty;
               target = tmptarget;
            end
        end
    end
    %show targetX and targetY in the command line
    targetX
    targetY
    colorTransform = makecform('lab2srgb');
    targetShow = applycform(target , colorTransform);
    highPriorityShow = applycform(highPriority , colorTransform);
    figure(7);imshow(targetShow);title('target');
    figure(8);imshow(highPriorityShow);title('highPriority');
    
    %STEP5 copy image data from origin image
    
    %ready to copy the target to the highPriority
    for i=-patchIndex:patchIndex
        for j=-patchIndex:patchIndex
            maskedImage(hpX+i, hpY+j,:) = target(1+patchIndex+i, 1+patchIndex + j,:);
            %update the mask;
            mask(hpX+i, hpY+j) = 0;
        end 
    end
    
    %for showing result
    tmpImage = maskedImage;
    colorTransform = makecform('lab2srgb');
    tmpImage = applycform(tmpImage , colorTransform);
    figure(9);imshow(tmpImage);title('tmpImage');
    figure(10);imshow(mask);title('after fill in mask');
    
    %figure;imshow(mask);title('mask after:');
    

    
    %STEP6 update the confidence
    
    % compute "confidence 
    confidence = double(1-mask);
    %figure;imshow(confidence);title('confidence');
    
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
    fprintf('HAHAHAHAHAHAHAHAHA\n');

end
colorTransform = makecform('lab2srgb');
maskedImage = applycform(maskedImage , colorTransform);
figure(3);imshow(maskedImage);

