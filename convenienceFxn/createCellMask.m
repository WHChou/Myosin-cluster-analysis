function mask = createCellMask(im, method)
%%%%%%%%%%
% Create cell mask from cytoskeletal images
% Automatic method uses Patrick's code in the tfmcode folder
% Input:
%         im - 2D matrix representing an image
%         method - 'manual' or 'auto'
% Output:
%         mask
% WHC 2021.08.11
%%%%%%%%%%
if strcmpi(method, 'manual')
    imshow(im, []);
    roi = drawpolygon();      % requires Matlab R2018b or later
    mask = createMask(roi);
    close;
elseif strcmpi(method, 'auto1')
    mask = makeCellmask(im);
elseif strcmpi(method, 'auto2')
    mask = makeCellmask2(im);
else
    error('Unknown mask creation method!');
end