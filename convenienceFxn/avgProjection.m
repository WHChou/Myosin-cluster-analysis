function im2d = avgProjection(im3d)
    [~, bp] = max(sum(sum(im3d)));
    if bp == 1
        im2d = mean(im3d(:,:,bp:bp+2), 3);
    elseif bp == size(im3d, 3)
        im2d = mean(im3d(:,:,bp-2:bp), 3);
    else
        im2d = mean(im3d(:,:,bp-2:bp+2), 3);
    end
end
