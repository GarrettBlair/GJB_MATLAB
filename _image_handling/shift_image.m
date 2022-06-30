function [im] = shift_image(im, h, w)
% Shifts input image matrix by height and width
if h>0
    im = [im(h+1:end, :); zeros(h, size(im,2))];
elseif h<0
    im = [zeros(abs(h), size(im,2)); im(1:end+h, :)];
end
if w>0
    im = [im(:, w+1:end) zeros(size(im,1), w)];
elseif w<0
    im = [zeros( size(im,1), abs(w)) im(:, 1:end+w)];
end
% if w>0
%     im = [ zeros(size(im,1), w) im(:, w+1:end)];
% elseif w<0
%     im = [ im(:, 1:end+w) zeros( size(im,1), abs(w))];
% end

end