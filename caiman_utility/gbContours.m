function contours = gbContours(A, dims, segs, spatial_threshold, resize_scale)
%% Input
% A = spatial contour matrix of size dims(1)*dims(2) x segs
% dims = height and width of the image
% segs = segments to make contours for
% spatial_threshold = cutoff for normalized pixel values , below are set to zero
%% Output
% contours = segs by dim(1) by dim(2)

if nargin < 5
    resize_scale = 1;
end
if ~exist('spatial_threshold', 'var') || isempty(spatial_threshold)
    spatial_threshold = 0;
end

if ~exist('segs', 'var') || isempty(segs)
    segs = 1:size(A, 2);
end
dims = double(dims);
newdims = ceil(dims*resize_scale);
nsegs = length(segs);
contours = NaN([max(segs),newdims(1),newdims(2)]);
% spatial_threshold = .5;
for j = 1:nsegs
    seg = segs(j);
    c = reshape(full(A(:,seg)), dims);
    c = imresize(c, resize_scale);
    if length(spatial_threshold)>1
        c(c<=(spatial_threshold(1)*max(c(:)))) = 0;
        c(c>=(spatial_threshold(2)*max(c(:)))) = 0;
    else
        c(c<(spatial_threshold*max(c(:)))) = 0;
    end
    if any(isnan(c(:)))
        c(isnan(c)) = 0;
    end
    contours(seg,:,:) = c;
end