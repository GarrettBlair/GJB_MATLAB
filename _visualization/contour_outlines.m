function [conts, conts_outline] = contour_outlines(A_contour_matix, min_spatial, min_outline, dims)
if nargin == 1
    % decent defaults
    min_spatial = .35;
    min_outline = .5;
end
conts = normalize_cols(A_contour_matix); 
conts(conts<min_spatial) = 0;
conts_outline = conts;
conts_outline(conts_outline>=min_outline) = 0;

if nargin>3
    conts = reshape(conts, [dims size(A_contour_matix, 2)]);
    conts = squeeze(sum(conts,3))>0;
    conts_outline = reshape(conts_outline, [dims size(A_contour_matix, 2)]);
    conts_outline = squeeze(sum(conts_outline,3))>0;
end


end