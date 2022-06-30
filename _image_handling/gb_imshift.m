function imOut = gb_imshift(imIn, r, c)

[h, w, d3] = size(imIn);

hv = 1:h;
wv = 1:w;

hshift = circshift(hv, r);
wshift = circshift(wv, c);

imOut = imIn*NaN;
for dim = 1:d3
    imOut(:,:,dim) = (imIn(hshift, wshift, dim));
end

