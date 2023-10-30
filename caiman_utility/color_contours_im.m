function [rgb_im] = color_contours_im(contours, cm, bg)
%%
if nargin<3
    bg = zeros(size(contours,2), size(contours,3));
end
[nsegs, h, w] = size(contours);
if ~isvarname('cm') || isempty(cm)
%     cm = jet(nsegs);
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
    rgb_im(:,:,1) = rgb_im(:,:,1) + c.*cm(i,1);
    rgb_im(:,:,2) = rgb_im(:,:,2) + c.*cm(i,2);
    rgb_im(:,:,3) = rgb_im(:,:,3) + c.*cm(i,3);
end
