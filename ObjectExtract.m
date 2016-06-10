I = im2double(imread('test3.jpeg'));
temp = im2double(imread('test3.jpeg'));
I = rgb2gray(I);
%I = I(:,:,3);
figure, imshow(I)
[y,x] = ginput(1);
x = round(x);
y = round(y);
J = regiongrowing(I,x,y,0.2); 
new = I+J;
[a,b,c] = size(I);
kk = zeros(a,b,3);
for i=1:a;
    for j=1:b;
        if(new(i,j)<1)
            kk(i,j,:)=temp(i,j,:);
        end
    end
end
figure, imshow(kk)

[y1,x1] = ginput(1);
x1 = round(x1);
y1 = round(y1);
[y2,x2] = ginput(1);
x2 = round(x2);
y2 = round(y2);
I = kk;
kk = zeros(a,b,3);
cow = zeros(a,b,3);
gcow = zeros(a,b);
%I = I(:,:,1) - 0.8*I(:,:,2);
for i=1:a;
    for j=1:b;
        if(I(i,j)<=0)
            cow(i,j,:)=temp(i,j,:);
            kk(i,j,:)=temp(246,192,:);
            %cow(i,j,:)= 0;
            gcow(i,j) = 0;
        else
            cow(i,j,:) = 255;
            gcow(i,j) = 1;
            kk(i,j,:)=temp(i,j,:);
        end
        if(j<y1 || j>y2 || i<x1 || i>x2)
            cow(i,j,:) = 255;
            gcow(i,j) = 1;
            kk(i,j,:) = temp(i,j,:);
        end
    end
end
figure, imshow(kk)
figure, imshow(cow)
imwrite(kk,'C:\Users\steven\Desktop\bg.jpg','jpg')
%cow = rgb2gray(cow);
imwrite(cow,'C:\Users\steven\Desktop\cow.jpg','jpg')
imwrite(gcow,'C:\Users\steven\Desktop\cowth.jpg','jpg')