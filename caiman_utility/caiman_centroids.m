function [centroids, area, max_pixel] = caiman_centroids(A, dims)
%% extract the peak location and centroid of the binary segment fotprint
[~, nsegs] = size(A);
centroids = NaN(nsegs, 2);
% weighted centroid, [x, y]
area = NaN(nsegs, 1);
% weighted centroid, [x, y]
max_pixel = NaN(nsegs, 2);
% average max pixel, [x, y]
for i = 1:nsegs
    suba = A(:,i);
    if any(suba)
        
    ind = find(max(suba)==suba);
    [x,y] = ind2sub(dims, ind);
    max_pixel(i,:) = round([mean(x), mean(y)]);
    ima = reshape(suba, dims);
%     temp = regionprops(ima>0, 'centroid');
%     ima = uint8(255*ima./max(ima(:)));
    temp = regionprops(ima>0, ima, 'WeightedCentroid', 'Area');
    centroids(i,:) = [temp.WeightedCentroid(2), temp.WeightedCentroid(1)];
    area(i) = [temp.Area];
    end
end