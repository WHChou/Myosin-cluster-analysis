function [pts, Erej] = filtClosePoints(pts, thresh)
%%%%
% Filter out particles that are too close to each other
% Use Delaunay triangulation, find triangle edges that are too short
% Input:
%        pts - A N x M matrix with (x,y) info of the points, where the
%                 first two columns are x and y
%        thresh - Threshold for how close between the points to be filtered
% Output:
%        pts1 - Filtered points array
% WHC 2021.08.11
%%%%%%%%

DT = delaunayTriangulation(pts(:, 1:2));
vert = DT.Points;
E = edges(DT);
Erej = [];
for jj = 1:length(E)      
    % Loop through all edges, find short edges
    x = vert(E(jj,:), 1); 
    y = vert(E(jj,:), 2); 
    dist = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    if dist < thresh
        rejID1 = find(pts(:,1)==x(1)&pts(:,2)==y(1));
        rejID2 = find(pts(:,1)==x(2)&pts(:,2)==y(2));
        Erej = [Erej; pts(rejID1, :); pts(rejID2, :)];
        pts(max([rejID1 rejID2]), :) = [];
        pts(min([rejID1 rejID2]), :) = [];
    end
end
