function [rgb_im] = make_color_im(bg, r, g, b)
%%
if isempty(bg)
    bg = zeros(size(r));
end
% rgb_im = cat(3, bg*0, bg*0, bg*0);
rgb_im = cat(3, bg, bg, bg);
if exist('r', 'var') && ~isempty(r)
    targ_diffs = size(bg) - size(r);
    targ_inds = [1+floor(targ_diffs(1)/2), 1+floor(targ_diffs(2)/2)];
    rgb_im(targ_inds(1):(targ_inds(1)+size(r,1)-1), targ_inds(2):(targ_inds(2)+size(r,2)-1), 1) = rgb_im(targ_inds(1):(targ_inds(1)+size(r,1)-1), targ_inds(2):(targ_inds(2)+size(r,2)-1), 1)+ r;
end

if exist('g', 'var') && ~isempty(g)
    targ_diffs = size(bg) - size(g);
    targ_inds = [1+floor(targ_diffs(1)/2), 1+floor(targ_diffs(2)/2)];
    rgb_im(targ_inds(1):(targ_inds(1)+size(g,1)-1), targ_inds(2):(targ_inds(2)+size(g,2)-1), 2) = rgb_im(targ_inds(1):(targ_inds(1)+size(g,1)-1), targ_inds(2):(targ_inds(2)+size(g,2)-1), 2) + g;
end

if exist('b', 'var') && ~isempty(b)
    targ_diffs = size(bg) - size(b);
    targ_inds = [1+floor(targ_diffs(1)/2), 1+floor(targ_diffs(2)/2)];
    rgb_im(targ_inds(1):(targ_inds(1)+size(b,1)-1), targ_inds(2):(targ_inds(2)+size(b,2)-1), 3) = rgb_im(targ_inds(1):(targ_inds(1)+size(b,1)-1), targ_inds(2):(targ_inds(2)+size(b,2)-1), 3) + b;
end