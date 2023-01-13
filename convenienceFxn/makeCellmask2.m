function [cmask] = makeCellmask2(transIm)
%%% Changed closing element in m4 to 10    -- WHC 08.02.21
%%% Changed m2 threshold to be zero if m-s/1.5 is negative  -- WHC 08.16.21

G2 = fspecial('gaussian', [30 30], 2);
f = ones(10,10);
box = ones(size(transIm));
edgebox = box-bwmorph(box,'erode',50);
edgebox = imcomplement(edgebox);

i1 = transIm;
%i2 = imclearborder(i1);
%i2 = imfilter(i1,G2);
i2 = wiener2(i1,[10 10]);
%i4 = imfilter(double(edge(i2,'canny')),f);
%i5 = immultiply(im2bw(i4),edgebox);

i3 = mat2gray(i2);
a = i3(:);
m = mean(a);
s = std(double(a));
upThresh = m+s/1.5;
lowThresh = (abs(m-s/1.5)+(m-s/1.5))/2;
m1 = im2bw(i3,upThresh);
m2 = imcomplement(im2bw(i3,lowThresh));
m3 = m1+m2;
m3 = imclearborder(im2bw(m3));
%figure
%imshow(m3)

m4 = imclose(m3,strel('disk',10));
m5 = imclearborder(im2bw(m4));
m6 = bwareaopen(m5,40000);
m7 = imclose(m6,strel('disk',50));
m8 = imfill(m7,'holes');
%figure
%imshow(m8)

%i7 = imclose(i6,strel('disk',50));
%i8 = imfill(i7,'holes');


cmask = m8;
end