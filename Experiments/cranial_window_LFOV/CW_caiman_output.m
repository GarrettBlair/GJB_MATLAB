ms.parentDir = 'D:\Sample Data\cortical_window_1P\megha_cw_mouse\mouse_cw2_rec1\';
% ms.parentDir = 'D:\Sample Data\cortical_window_1P\megha_cw_mouse\mouse_cw1_rec1\';
ms.tiffName = [ms.parentDir 'dff.tiff'];
% ms.manROIs = readmatrix([ms.parentDir 'manual_points.csv']);
% ms.manROIs = ms.manROIs(:,6:7);
% ms.midline = readmatrix([ms.parentDir 'midline_vals.csv']);
% ms.midline = ms.midline(:,2);
caimanFilename = sprintf('%scaiman_cnmfe_out.mat', ms.parentDir);
sminSweepFilename = sprintf('%sdeconv_sweep.mat', ms.parentDir);
caiman_data = load(caimanFilename);




%% m = movmedian(ms.midline, 5); m = m-median(m);
params.smin_vals = [-30:5:-5];
params.remove_bad_caiman_segs = true;
params.reuse_contour_crop = 'bounding_box.mat';
% Cull contours
% correct for python indexing start at 0
% caiman_data.idx_components_bad = caiman_data.idx_components_bad+1;
% caiman_data.idx_components = caiman_data.idx_components+1;
[nsegs,nframes] = size(caiman_data.C);

[smat, smat_weighted, good_idx, ~] = deconv_sweep_read(sminSweepFilename, params.smin_vals);
if params.remove_bad_caiman_segs
    % all_good_idx = find(sum(good_idx,1)==size(good_idx,1));
    all_good_idx = find(sum(good_idx,1)>0);
    bad_idx = setdiff(1:size(caiman_data.C,1), all_good_idx);
    caiman_data.idx_components = all_good_idx;
    caiman_data.idx_components_bad = bad_idx;
else
    caiman_data.idx_components = 1:size(caiman_data.C,1);
    caiman_data.idx_components_bad = [];
end

temp = sum(smat, 1);
caiman_data.S_mat = reshape(temp, [nsegs, nframes]);
temp = sum(smat_weighted, 1);
caiman_data.S_matw = reshape(temp, [nsegs, nframes]);

caiman_data.fullA = single(zeros(caiman_data.dims(1)*caiman_data.dims(2), size(caiman_data.C,1)));
for aloop = 0:100:floor(size(caiman_data.C,1)/100)*100
    temp = eval(sprintf('caiman_data.fullA%d', aloop));
    ns = size(temp,2);
    i1=aloop+1; i2 = aloop+ns;
    caiman_data.fullA(:, i1:i2) = single(temp);
    caiman_data = rmfield(caiman_data, sprintf('fullA%d', aloop));
end
clearvars temp
%

if ~isempty(params.reuse_contour_crop)
    tempCropName = sprintf('%s/%s', ms.parentDir, params.reuse_contour_crop);
%     tempCropName = sprintf('%s/MiniLFOV/%s', ms.parentDir, params.reuse_contour_crop);
    load(tempCropName, 'valid_contour_bounds');
    if exist('valid_contour_bounds', 'var')
        draw_bounds = false;
        nsegs = size(caiman_data.C,1);
        good_flag = true(nsegs,1);
        for j = 1:nsegs
            a = (reshape(caiman_data.A(:,j), [caiman_data.dims]))>0;
            [yy, xx] = ind2sub(size(a), find(a));
            isgood = inpolygon(xx, yy, valid_contour_bounds.x, valid_contour_bounds.y);
            prop_in_poly = sum(isgood)/length(isgood);
            if prop_in_poly < .5
                good_flag(j) = false;
            end
        end
        bad_inds = find(~good_flag);
        allbad = unique([caiman_data.idx_components_bad, bad_inds']);

    else
        draw_bounds = true;
    end
else
%     tempCropName = sprintf('%s/MiniLFOV/%s', ms.parentDir, 'bounding_box.mat');
    tempCropName = sprintf('%s/%s', ms.parentDir, 'bounding_box.mat');
    draw_bounds = true;
end
if draw_bounds
    [~, bad_inds, ~, valid_contour_bounds] = Draw_contour_bounding(caiman_data.fullA, ...
        caiman_data.dims, caiman_data.maxFrame, caiman_data.idx_components, false);
    save(tempCropName, 'valid_contour_bounds')
    allbad = unique([caiman_data.idx_components_bad, bad_inds']);
end
fprintf('\nRemoving %d bad components\n', length(allbad))
neuron = remove_segments(caiman_data, allbad, false);

% S2 = S2(ms.neuron.idx_components,:);


% neuron = caiman_data;
ms.neuron = neuron;
ms.valid_contour_bounds = valid_contour_bounds;

%%

xseg = seg_peak(:,2);
yseg = seg_peak(:,1);
xm = round(ms.manROIs(:,1));
ym = round(ms.manROIs(:,2));

Y = imread_big(ms.tiffName);
dff_t = zeros(size(ms.manROIs,1),size(Y,3));
craw_t = zeros(size(ms.manROIs,1),size(Y,3));
for j = 1:size(ms.manROIs,1)
    x = xm(j);
    y = ym(j);
    dd = sqrt((xseg-xm(j)).^2 + (yseg-ym(j)).^2 );
    [~, nearestord] = sort(dd, 'ascend');
    xs = x-2:x+2;
    xs = xs(xs>0 & xs<size(Y,2));
    ys = y-2:y+2;
    ys = ys(ys>0 & ys<size(Y,1));
    dff_t(j,:) = squeeze(mean(mean(Y(ys,xs,:))))';
    craw_t(j,:) = craw(nearestord(1),:);
end
%%
im = neuron.maxFrame./255/2;
im2 = neuron.minFrame./255/2;

spks = normalize_rows(ms.neuron.S_matw);
craw = normalize_rows(neuron.C+neuron.YrA);

rng('default')
nsegs = size(neuron.fullA, 2);
segs = randi(nsegs, [20, 1]);

[seg_cent, seg_area, seg_peak] = caiman_centroids(neuron.fullA, neuron.dims);
neuron.seg_center = seg_cent;
neuron.seg_area = seg_area;
neuron.seg_peak = seg_peak;
neuron = rmfield(neuron, 'mc_xshifts');
neuron = rmfield(neuron, 'mc_yshifts');
neuron = rmfield(neuron, 'dataIso');
neuron.A = neuron.fullA;
neuron = rmfield(neuron, 'fullA');
neuron = rmfield(neuron, 'S_mat');
neuron = rmfield(neuron, 'S_matw');
[contours] = gbContours(neuron.fullA, neuron.dims, [], .5);
bad = seg_area>=800 | seg_area<=150;
% cm = ones(nsegs,3);
% cm(bad, [2 3]) = 0;
% cm(~bad, [1 3]) = 0;
% r = squeeze(sum(contours(bad,:,:)));
% g = squeeze(sum(contours(~bad,:,:)));
% rgb_im = cat(3, r, g, g*0);
% % % y1 = 640; y2 = y1+200;
% % % x1 = 327; x2 = x1+200;
y1 = 250; y2 = y1+200;
x1 = 1050; x2 = x1+200;
cm =  plasma(nsegs);
[~,ord] = sort(rand(nsegs,1));
[rgb_im] = color_contours_im(contours,cm(ord,:));
%%
mic_per_px = 3.5;
half_mm = round(500/mic_per_px);


subrgb = rgb_im(y1:y2, x1:x2,:);
figure(1); clf
subplot_tight(3,2,1)
image(rgb_im*3 + im); axis image
hold on
rectangle('Position', [x1 y1 x2-x1 y2-y1], 'EdgeColor', 'r')
plot([1050, 1050+half_mm], [900 900], 'w', 'LineWidth', 5);
subplot_tight(3,2,2)
bg = im(y1:y2, x1:x2);%-im2(640:640+100, 327:327+100);
bg = normalize_matrix(bg);
image(subrgb*3 + bg); axis image
hold on
plot([50, 50+half_mm], [190 190], 'w', 'LineWidth', 5);
rectangle('Position', [1 1 x2-x1 y2-y1], 'EdgeColor', 'r', 'LineWidth', 3)

subplot_tight(3,1,2:3, [.05 .1])
imagesc(craw); colormap hot
% stacked_traces(craw(segs,:), .8, {'k'});
axis tight
colorbar
fr = 1/11;
t = 0:fr:size(neuron.C,2);
mt = mod(t, 150);
% segs2plot = [0:300:size(neuron.C,1)-100, size(neuron.C,1)]+.5;
segs2plot = [0:10:50]+.5;
set(gca, 'XTick', find(mt==0), 'XTickLabel', t(mt==0));
xlabel('Seconds')
set(gca, 'YTick', segs2plot, 'YTickLabel', segs2plot-.5);
ylabel('Cell ID')
% figure(2); clf
% imagesc(craw); colormap hot
% stacked_traces(craw, 5, {'k'});
% axis tight
