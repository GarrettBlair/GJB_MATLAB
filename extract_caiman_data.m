function [ms] = extract_caiman_data(ms, params, cameraName)

caimanFilename = sprintf('%s/%s/caiman_cnmfe_out.mat', ms.parentDir, cameraName);
sminSweepFilename = sprintf('%s/%s/deconv_sweep.mat', ms.parentDir, cameraName);
caiman_data = load(caimanFilename);

%% Cull contours
% correct for python indexing start at 0
% caiman_data.idx_components_bad = caiman_data.idx_components_bad+1;
% caiman_data.idx_components = caiman_data.idx_components+1;
[nsegs,nframes] = size(caiman_data.C);


if isfile(sminSweepFilename)
    [smat, smat_weighted, good_idx, ~] = deconv_sweep_read(sminSweepFilename, params.smin_vals);
    all_good_idx = find(sum(good_idx,1)>0);
    bad_idx = setdiff(1:size(caiman_data.C,1), all_good_idx);
    caiman_data.good_idx_smat = good_idx;
    caiman_data.idx_components = all_good_idx;
    caiman_data.idx_components_bad = bad_idx;
    
    temp = sum(smat, 1);
    caiman_data.S_mat = reshape(temp, [nsegs, nframes]);
    temp = sum(smat_weighted, 1);
    caiman_data.S_matw = reshape(temp, [nsegs, nframes]);
else
    caiman_data.good_idx_smat = caiman_data.idx_components+1;
    caiman_data.idx_components = caiman_data.idx_components+1;
    caiman_data.idx_components_bad = caiman_data.idx_components_bad+1;
    caiman_data.S_mat = false(nsegs, nframes);
    caiman_data.S_matw = false(nsegs, nframes);
end

if length(ms.timestamps) ~= nframes || length(ms.frameNum) ~= nframes
    warning('!~!~! Frame number discrepancy found between ms and caiman files!')
    disp([ms.fileName])
    last_ts = length(ms.timestamps);
    if last_ts < nframes && isfield(ms.warnings, 'TrackerCrash')
        disp([ms.warnings.TrackerCrash])
        fprintf('~~~Tracker crashed,\n\tconcatenating caiman_data at %d index\n', last_ts);
        caiman_data.C               = caiman_data.C(:, 1:last_ts);
        caiman_data.S               = caiman_data.S(:, 1:last_ts);
        caiman_data.YrA             = caiman_data.YrA(:, 1:last_ts);
        caiman_data.mc_xshifts      = caiman_data.mc_xshifts(1:last_ts,:);
        caiman_data.mc_yshifts      = caiman_data.mc_yshifts(1:last_ts,:);
%         caiman_data.dataIso         = caiman_data.dataIso(1:last_ts, :);
        caiman_data.S_mat           = caiman_data.S_mat(:, 1:last_ts);
        caiman_data.S_matw          = caiman_data.S_matw(:, 1:last_ts);
    else
        error('Uknown solution')
    end
end
%% Cropping contours manually if needed
if ~isempty(params.reuse_contour_crop)
    tempCropName = sprintf('%s%s/%s', ms.parentDir, cameraName, params.reuse_contour_crop);
    if isfile(tempCropName)
%     tempCropName = sprintf('%s/MiniLFOV/%s', ms.parentDir, params.reuse_contour_crop);
    prev_judgement = load(tempCropName, 'valid_contour_bounds');
    else
        prev_judgement = [];
    end
%     prev_judgement = [];
    if isfield(prev_judgement, 'valid_contour_bounds')
        valid_contour_bounds = prev_judgement.valid_contour_bounds;
        draw_bounds = false;
        nsegs = size(caiman_data.C,1);
        good_flag = true(nsegs,1);
        for j = 1:nsegs
            a = (reshape(caiman_data.fullA(:,j), [caiman_data.dims]))>0;
            [yy, xx] = ind2sub(size(a), find(a));
            isgood = inpolygon(xx, yy, valid_contour_bounds.x, valid_contour_bounds.y);
            prop_in_poly = sum(isgood)/length(isgood);
            if prop_in_poly < .5
                good_flag(j) = false;
            end
        end
        bad_inds = find(~good_flag);
        allbad = unique([caiman_data.idx_components_bad, bad_inds']);
        
        contours = gbContours(caiman_data.fullA, caiman_data.dims, [], .5);
        background = caiman_data.minFrame./255;
        bc = squeeze(sum(contours(~good_flag, :, :), 1));
        mc = squeeze(sum(contours(good_flag, :, :), 1));
        figure;
        x0 = 100;
        set(gcf, 'Name', ms.fileName, 'Position', [x0 x0 x0+2*caiman_data.dims(2), x0+2*caiman_data.dims(1)])
%         subplot_tight(1,2,2)
        im = zeros(caiman_data.dims(1), caiman_data.dims(2), 3);
        im(:,:,1) = background*.9 + bc;
        im(:,:,2) = background*.9 + mc;
        im(:,:,3) = background*.9;
%         im = im(bounds.mh1:bounds.mh2, bounds.mw1:bounds.mw2, :);
        
        image(im); hold on
        plot(valid_contour_bounds.x, valid_contour_bounds.y, 'mo-')
        % plot([opx opx(1)]-bounds.mw1, [opy opy(1)]-bounds.mh1, 'mo-')
        title(sprintf('Inside: %d  ,   Outside: %d', sum(good_flag), sum(~good_flag)))
        axis image
        drawnow;
        
    else
        draw_bounds = true;
    end
else
    tempCropName = sprintf('%s/%s/%s', ms.parentDir, cameraName, 'bounding_box.mat');
    draw_bounds = true;
end
if draw_bounds
    [~, bad_inds, ~, valid_contour_bounds] = Draw_contour_bounding(caiman_data.fullA, ...
        caiman_data.dims, caiman_data.maxFrame, caiman_data.idx_components, params.skip_contour_bounding);
    save(tempCropName, 'valid_contour_bounds')
end



if params.remove_all_bad_caiman_segs == true
    allbad = unique([caiman_data.idx_components_bad, bad_inds']);
    fprintf('\nRemoving %d bad components and outside components\n', length(allbad))
    neuron = remove_segments(caiman_data, allbad, false);
elseif params.remove_segs_outside == true
    allbad = bad_inds;
    fprintf('\nRemoving %d components outside drawing\n', length(allbad))
    neuron = remove_segments(caiman_data, allbad, false);
else
    neuron = caiman_data;
end

if isfield(ms ,'goodFrames') && (length(ms.frameNum) - length(ms.goodFrames) > 0)
    fprintf('\nRemoving %d bad frames\n', length(ms.frameNum) - length(ms.goodFrames))
    [neuron] = caiman_downsample(neuron, ms.goodFrames, false);
else
    fprintf('\nNo bad frames found/removed\n')
end
ms.neuron = neuron;

ms.valid_contour_bounds = valid_contour_bounds;
