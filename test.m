clear;
clear all;
%imageName = 'bungee.png';
imageName = 'air';
originImage = imread([imageName, '.jpg']);

%�p�G�n���emask�A��H�U�����Ѩ���
%if you want to redraw mask, cancel the comment below
%mask = im2uint8(roipoly(originImage));
%imwrite(mask, [imageName, '-mask.png']);
mask = imread([imageName, '-mask.png']);


%figure;imshow(img);
figure;imshow(mask);title('mask');
[m,n] = size(mask);
boundary = zeros(m,n);

%patch���j�p
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


% �p��boundary
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

%index ��mask���Ҧ��I�y��
%index is all coordinate(x,u) in mask
index = [];
for i=1:m
    for j=1:n
        if mask(i,j) == 255
            index = [index; i j];
        end
    end
end
figure;imshow(boundary);title('boundary');

%�p��data(�����M��gradient)
%compute "data"(use matlab function "gradient")
grayOriginImage = rgb2gray(originImage);
grayOriginImage = im2single(grayOriginImage);
[fx,fy] = gradient(grayOriginImage);
gradtnorm = sqrt((fx .^2) + (fy .^2));
data = gradtnorm + log10(1+gradtnorm);
%figure;imshow(fx);
%figure;imshow(fy);
%figure;imshow(data);title('data');



%�p��confidence (patch�̭���ϻPmask�����)
% compute "confidence 
confidence = double(1-mask);
%figure;imshow(confidence);title('confidence');

%indexM�Oindexy�`�ƥ�
%indexM is the length of index

[indexM, indexN] = size(index);

patchIndex = floor(patch/2);

for i=1:indexM %�C��index�]�@�� %run every index
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
%figure;imshow(priority);title('priority');

%�}�linpainting�F
%START TO INPAINTING

while ~isequal(mask, zeros(m,n) ) %as research ���� omega = 0����
                                  %omega �N�O�o�̪�mask
                                  %as sesearch, until omdga = 0
                                  %omega is mask here
    %��̤j��priority
    %find maximum priority
    max = -1;
    
    %high priority X, Y �y��
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

    %�M�whighPriority �O����
    %find highPriority
    highPriority = uint8(zeros(patch, patch, 3));
    for i=1:patch
       for j=1:patch
           highPriority(i,j,:) = source(hpX+i-patchIndex-1, hpY+j-patchIndex-1, :);
       end
    end
    %figure;imshow(highPriority);title('highPriority');
    
    
    %��̬ۦ�������
    %find the most similar target(patch X patch)
    
    
    %�קKtarget �I�������
    min = 99999999;
    %to avoid "target" hit the border of originImage
    targetX = 0;
    targetY = 0;
    target = uint8(zeros(patch, patch, 3));
    for tx=hpX-floor((patch*patch)/2):hpX+floor((patch*patch)/2)
        %out of boundary
        if tx < 1 || tx > m
            continue;
        end
        for ty=hpY-floor((patch*patch)/2):hpY+floor((patch*patch)/2)
            %out of boundary
            if ty < 1 || tx > n
                continue;
            end
            
            %in mask
            %inside (tx, ty) patch
            %there can't be any mask == white inside
            ismask = 0;
            for mi=tx-patchIndex:tx+patchIndex
                for mj=ty-patchIndex:ty+patchIndex
                    if mask(mi,mj) == 255
                        ismask = 1;
                    end
                end
            end
            if ismask == 1
                continue;
            end
            
            
            SSD = 0;
            tmptarget = uint8(zeros(patch, patch, 3));
            %fill the target
            for px=-patchIndex:patchIndex
                for py=-patchIndex:patchIndex
                    tmptarget(px+patchIndex+1, py+patchIndex+1, :) = originImage(tx+px, ty+py,:);
                    %target(1~patch)
                end
            end
            %figure;imshow(tmptarget);title('tmptarget');
            
            %compute SSD
            for x=1:patch
                for y=1:patch
                    R = highPriority(x,y,1) - tmptarget(x,y,1); 
                    G = highPriority(x,y,2) - tmptarget(x,y,2);
                    B = highPriority(x,y,3) - tmptarget(x,y,3);
                    SSD = SSD + R*R + G*G + B*B;
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
    figure;imshow(target);title('target');
    fprintf('HAHAHAHAHAHAHAHAHA\n');
    

end

