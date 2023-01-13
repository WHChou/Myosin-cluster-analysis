function [cmask] = makeCellmask(transIm)

G2 = fspecial('gaussian', [50 50], 2);
f = ones(10,10);
box = ones(size(transIm));
edgebox = box-bwmorph(box,'erode',50);
edgebox = imcomplement(edgebox);

i1 = transIm;
%i2 = imclearborder(i1);
i2 = wiener2(i1,[10 10]);

i3 = imfilter(i2,G2);
i3c = imcrop(i3,[1,1,250,250]);
%i3 = i3 - median(median(i3c))*1.1;

i4 = imfilter(double(edge(i3,'canny')),f);
i5 = immultiply(im2bw(i4),edgebox);


i6 = bwareaopen(i5,40000);

% figure
% imshow(i6)

i7 = imclose(i6,strel('disk',50));
i7 = imclose(i7,strel('disk',10));
i8 = imfill(i7,'holes');

 figure
 imshowpair(i6,i8)
 pause(2);
 close;

cmask = i8;
end