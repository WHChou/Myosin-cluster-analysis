function plotPuncta(img, acceptedPuncta, rejectedPuncta, closePuncta, markerSz)
    figure; imshow(img, []); hold on;
    plot(acceptedPuncta(:,1), acceptedPuncta(:,2), 'r.', 'MarkerSize', markerSz)
    if ~isempty(rejectedPuncta)
        plot(rejectedPuncta(:,1), rejectedPuncta(:,2), 'g.', 'MarkerSize', markerSz)
    end
    if ~isempty(closePuncta)
        plot(closePuncta(:,1), closePuncta(:,2), 'c.', 'MarkerSize', markerSz)
    end
end
