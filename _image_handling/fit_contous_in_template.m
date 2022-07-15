function [contours_out] = fit_contous_in_template(contours, contour_bounds, template_size)
%%
% contours = gbContours(ms.neuron.fullA, ms.neuron.dims, [], .6);
% contour_bounds = temp_params.params.crop_params.cropROI;
% % [x1 y1; x2 y2]
% template_size = [camData.ROI.height, camData.ROI.width];

[nsegs, h, w] = size(contours);

c = contour_bounds;
contours_out = zeros([nsegs, template_size]);
contours_out(:, c(1,2):c(2,2)-1, c(1,1):c(2,1)-1) = contours;
