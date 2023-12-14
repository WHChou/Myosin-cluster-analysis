function [onFiber, fiberLength] = isOnFiber(points1, fibers, imActin)
    onFiber = [];  % Placeholder for puncta info
    fiber_id = unique(fibers(:,3));
    fiberLength = zeros(1, size(fiber_id, 1));
    for id = 1:length(fiber_id)   % For each fiber
        oneFiber = fibers(fibers(:,3)==fiber_id(id), :);  % One fiber
        L=0;  % Placeholder for whole fiber length
        for ct = 1:size(oneFiber, 1)-1  % For each fiber segment in a single fiber
            %ct
            xi = [oneFiber(ct, 1), oneFiber(ct+1, 1)];  %% Two end points of a filament segment
            yi = [oneFiber(ct, 2), oneFiber(ct+1, 2)]; 
            % Calculate length of fiber segment
            par_L = sqrt((xi(2)-xi(1))^2+(yi(2)-yi(1))^2);
            L = L+par_L;
            % See if any puncta belongs to that segment
            for ii = 1:length(points1)
                xp = points1(ii,1);
                yp = points1(ii,2);
                flag = 1;
                if abs((xi(2)-xi(1))*(yi(1)-yp)-(xi(1)-xp)*(yi(2)-yi(1)))/sqrt((xi(2)-xi(1))^2+(yi(2)-yi(1))^2) > 3
                    % distance from point to line segment > (3 pixels+linewidth/2)
                    flag = 0;
                end
                if (xi(2)-xi(1))*(xp-xi(1))+(yi(2)-yi(1))*(yp-yi(1))<0
                    % If point is outside line segment
                    flag = 0;
                end
                if ((xp-xi(1))^2+(yp-yi(1))^2)>par_L^2
                    % If point is outside line segment
                    flag = 0;
                end
                if flag
                    r = feature2Dcore(imActin, 1, 3, points1(ii,1), points1(ii,2));
                    acInt = r(:,6);
                    %imActin(round(points1(ii,2)), round(points1(ii,1)));
                    points = [points1(ii, :), acInt, fiber_id(id)];
                    onFiber = [onFiber; points];
                    %points1(ii, 1:2)
                end
            end
        end
        fiberLength(id) = L;
    end
end