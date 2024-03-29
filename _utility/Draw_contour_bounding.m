function [good_inds, bad_inds, contours, bounds] = Draw_contour_bounding(A, dims, background, idx_components, skip_flag)
%%
quit_radius = 20;
nsegs = size(A,2);
segID = 1:size(A,2);
% % clear
% % load('G:\Miniscope Data\Hipp6\2020_03_05\16_55_16_linear1\Miniscope\caiman_cnmfe_out3.mat')
% % background = maxFrame./255;

if issparse(A)
    A = full(A);
end
A = normalize_cols(A);
contour_draw_thresh = .6;
contours = gbContours(A, dims, [], contour_draw_thresh);

%
if ~exist('background', 'var') || isempty(background)
    background = zeros(dims)+.3;
else
    background = background - min(background(:));
    background = background/max(background(:));
end
%%
if exist('idx_components', 'var') && ~isempty(idx_components)
%     mc = squeeze(sum(contours(idx_components, :, :), 1));
    bmap = viridis(length(idx_components)*4);
    bmap = bmap(randperm(length(idx_components)*4), :);
    bmap(:,1) = 0;
    bmap = ones(length(idx_components),1)*[0 1 .3];
    mc = color_contours_im(contours(idx_components, :, :), bmap);
    idx_components_bad = setdiff(segID, idx_components);
%     bmap = plasma(length(idx_components_bad)*4);
%     bmap = bmap(randperm(length(idx_components_bad)*4), :);
%     bmap(:,2) = 0;
    bmap = ones(length(idx_components_bad),1)*[1 0 .3];
    bc = color_contours_im(contours(idx_components_bad, :, :), bmap );
else
    if ~any(background(:))
        mc = squeeze(sum(contours, 1))>0;
        bc = zeros(dims);
    else
        mc = squeeze(sum(contours, 1));
        bc = zeros(dims);
    end
end

if ~exist('skip_flag', 'var') || isempty(skip_flag)
    skip_flag = false;
end
bc = bc*3; mc = mc*3;
%
im = zeros(dims(1), dims(2), 3);
im = background*.8 + bc*.1 + mc*.1;
% im(:,:,1) = background*.9 + bc;
% im(:,:,2) = background*.9 + mc;
% im(:,:,3) = background*.9;
%
h = figure(375); clf; 
set(gcf, 'Position', [224 183 1332 769])
subplot_tight(1,2,1)
image(im);
axis image
title('Draw polygon! Click near first point to end, q=quit, r=restart')
hold on
ind = 1;
opx = [];
opy = [];
%%
if ~skip_flag
    [opx(ind, 1), opy(ind, 1), button] = ginput(1);
else
    button = 113;
end
if button~=113 % if not q get points
    plot(opx(ind), opy(ind), 'mo')
    rectangle('Position', [opx(ind)-quit_radius/2, opy(ind)-quit_radius/2, quit_radius, quit_radius],'Curvature',[1 1], 'EdgeColor', 'r')
    pd = 1000;
    button = 1;
    %%
    while pd >=quit_radius/2 && button~=113 % q
        ind = ind+1;
        [opx(ind), opy(ind), button] = ginput(1);
        if button == 114 % r
            % restart
            figure(375); clf; subplot_tight(1,2,1)
            image(im);
            axis image
            title('Draw polygon! Click near first point to end, q=quit, r=restart')
            hold on
            ind = 1;
            opx = [];
            opy = [];
            [opx(ind), opy(ind)] = ginput(1);
            plot(opx(ind), opy(ind), 'mo')
            rectangle('Position', [opx(ind)-quit_radius/2, opy(ind)-quit_radius/2, quit_radius, quit_radius],'Curvature',[1 1], 'EdgeColor', 'r')
            pd = 1000;
            ind = ind+1;
            [opx(ind), opy(ind), button] = ginput(1);
        elseif button==113 && length(opx)==1
            opx = double([1    1       dims(2)     dims(2), 1]);
            opy = double([1    dims(1) dims(1)     1,       1]);
            plot(opx, opy, 'mo-', 'LineWidth', 3)
        end
        plot(opx(ind-1:ind), opy(ind-1:ind), 'mo-')
        pd = sqrt( (opx(1) - opx(ind)).^2 + (opy(1) - opy(ind)).^2 );
    end
else % else use whole FOV, no contours excluded
    opx = double([1    1       dims(2)     dims(2), 1]);
    opy = double([1    dims(1) dims(1)     1,       1]);
    plot(opx, opy, 'mo-', 'LineWidth', 3)
end
good_flag = true(nsegs,1);
for j = 1:nsegs
    a = (reshape(A(:,j), [dims]))>contour_draw_thresh;%0;
    
    [yy, xx] = ind2sub(size(a), find(a));
    isgood = inpolygon(xx, yy, opx, opy);
%     if ~any(a(:)) % if no pixles are left, it's bad
%         good_flag(j) = false;
%     end
    prop_in_poly = sum(isgood)/length(isgood);
    if prop_in_poly < .5
        good_flag(j) = false;
    end
end
%
edge = 20;
mc = squeeze(sum(contours(good_flag, :, :), 1));
bounds.x = opx;
bounds.y = opy;
bounds.mh1 = max(find(sum(mc,2), 1, 'first')-edge, 1);
bounds.mh2 = min(find(sum(mc,2), 1, 'last')+edge, size(mc,1));
bounds.mw1 = max(find(sum(mc,1), 1, 'first')-edge, 1);
bounds.mw2 = min(find(sum(mc,1), 1, 'last')+edge, size(mc,2));

% bc = squeeze(sum(contours(~good_flag, :, :), 1));
idx_components = find(good_flag);
bmap = viridis(length(idx_components)*4);
bmap = bmap(randperm(length(idx_components)*4), :);
bmap(:,1) = 0;
mc = color_contours_im(contours(idx_components, :, :), bmap);
idx_components_bad = setdiff(segID, idx_components);

if any(idx_components_bad)
%     bmap = plasma(length(idx_components_bad)*4);
%     bmap = bmap(randperm(length(idx_components_bad)*4), :);
%     bmap(:,2) = 0;
    bmap = ones(length(idx_components_bad),1)*[1 0 .3];
    bc = color_contours_im(contours(idx_components_bad, :, :), bmap );
else
    bc = 0;
end

figure(h);
subplot_tight(1,2,2)
bc = bc*3; mc = mc*3;
%
im = zeros(dims(1), dims(2), 3);
im = background*.8 + bc*.1 + mc*.1;

% im(:,:,1) = background*.9 + bc;
% im(:,:,2) = background*.9 + mc;
% im(:,:,3) = background*.9;
% im = im(bounds.mh1:bounds.mh2, bounds.mw1:bounds.mw2, :);

image(im); hold on
% plot([opx opx(1)]-bounds.mw1, [opy opy(1)]-bounds.mh1, 'mo-')
title(sprintf('Inside: %d  ,   Outside: %d', sum(good_flag), sum(~good_flag)))
rectangle('Position', [ bounds.mw1 bounds.mh1 bounds.mw2-bounds.mw1 bounds.mh2-bounds.mh1 ], 'EdgeColor', 'g')
% axis([ 0 bounds.mw2-bounds.mw1 0 bounds.mh2-bounds.mh1])
axis image
drawnow;

good_inds = find(good_flag);
bad_inds = find(~good_flag);
poly_vals = [opx, opy];











