function [rgb_im] = color_contours_im(contours, cm, bg, c_scale)
%%
[nsegs, h, w] = size(contours);
if nargin<2 || isempty(cm)
    cm = rand(nsegs, 3);
    cm = (cm+.5)./2;
    [~, ord] = sort(rand(nsegs,1));
    cm = cm(ord,:);
end
if nargin<3
    bg = zeros(size(contours,2), size(contours,3));
end
if nargin<4
    c_scale = 1;
end
if ~isvarname('cm') || isempty(cm)
    cm = rand(nsegs, 3);
    cm = (cm+.5)./2;
    [~, ord] = sort(rand(nsegs,1));
    cm = cm(ord,:);
end
bg_im = cat(3, bg, bg, bg);
rgb_im = zeros(h, w, 3) + bg_im;

for i = 1:nsegs
    c = squeeze(contours(i,:,:));
%     c = 300.*c./max(c(:));
    rgb_im(:,:,1) = rgb_im(:,:,1) + c.*cm(i,1).*c_scale;
    rgb_im(:,:,2) = rgb_im(:,:,2) + c.*cm(i,2).*c_scale;
    rgb_im(:,:,3) = rgb_im(:,:,3) + c.*cm(i,3).*c_scale;
end
rgb_im(rgb_im>1) = 1;
