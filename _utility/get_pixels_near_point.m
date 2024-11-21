function [valid_inds, dists] = get_pixels_near_point(h, w, x, y, radius)
im = true(h,w);

idx = find(im);
[col, row] = ind2sub([h,w], idx);

dists = sqrt((x-col).^2 + (y-row).^2);
valid_inds = dists < radius;

end