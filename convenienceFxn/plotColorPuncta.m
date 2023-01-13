function [] = plotColorPuncta(im, points, baseInt, inverted, showcolorbar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates a map of myosin puncta color-coded by intensity, 
%   overlaid onto the origin myosin image
% Inputs:
%   im: original myosin image 
%   points: lists of puncta information (output from feature2D)
%   baseInt: intensity intervals at which color changes
%   inverted: if 1, the myosin image will be inverted
% Output:
%   A Matlab figure with color-coded myosin puncta overlaid onto myosin
%   image
% Wen-hung Chou, 2022/04/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 5
        showcolorbar = 0;
    end
    
    if inverted
        figure; imshow(imcomplement(im), []); hold on
    else
        figure; imshow(im, []); hold on
    end
    
    co = viridis(10);
    for ii = 1:10
        flag = min(floor(points(:,6)/baseInt)+ 1, 10) == ii;
        plot(points(flag,1), points(flag,2), '.', 'Color', co(ii, :), 'MarkerSize', 8); hold on
    end
     
    if showcolorbar
        figure;
        colormap viridis
        colorbar('Ticks', linspace(0,1,10), 'TickLabels', ...
            {['<', int2str(baseInt)], int2str(baseInt*2), int2str(baseInt*3), int2str(baseInt*4), ...
            int2str(baseInt*5), int2str(baseInt*6), int2str(baseInt*7), int2str(baseInt*8), int2str(baseInt*9), ...
            ['>', int2str(baseInt*10)]});
        set(gca, 'FontSize', 20)
    end
end