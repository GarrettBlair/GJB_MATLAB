function [ms] = APA_generate_placemaps(ms, params)
caimanFilename = sprintf('%s/MiniLFOV/caiman_cnmfe_out.mat', ms.parentDir);
sminSweepFilename = sprintf('%s/MiniLFOV/deconv_sweep.mat', ms.parentDir);
caiman_data = load(caimanFilename);

%% Cull contours
% correct for python indexing start at 0
% caiman_data.idx_components_bad = caiman_data.idx_components_bad+1;
% caiman_data.idx_components = caiman_data.idx_components+1;
[nsegs,nframes] = size(caiman_data.C);



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
        caiman_data.dataIso         = caiman_data.dataIso(1:last_ts, :);
        caiman_data.S_mat           = caiman_data.S_mat(:, 1:last_ts);
        caiman_data.S_matw          = caiman_data.S_matw(:, 1:last_ts);
    else
        error('Uknown solution')
    end
end
%% Cropping contours manually if needed
if ~isempty(params.reuse_contour_crop)
    tempCropName = sprintf('%s/MiniLFOV/%s', ms.parentDir, params.reuse_contour_crop);
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
    tempCropName = sprintf('%s/MiniLFOV/%s', ms.parentDir, 'bounding_box.mat');
    draw_bounds = true;
end
if draw_bounds
    [~, bad_inds, ~, valid_contour_bounds] = Draw_contour_bounding(caiman_data.fullA, ...
        caiman_data.dims, caiman_data.maxFrame, caiman_data.idx_components, params.skip_contour_bounding);
    save(tempCropName, 'valid_contour_bounds')
    allbad = unique([caiman_data.idx_components_bad, bad_inds']);
end

if params.remove_bad_caiman_segs
    fprintf('\nRemoving %d bad components\n', length(allbad))
    neuron = remove_segments(caiman_data, allbad, false);
else
    neuron = caiman_data;
end

ms.neuron = neuron;
ms.valid_contour_bounds = valid_contour_bounds;
spks = normalize_rows(ms.neuron.S_matw);
% temp_ts = cat(1, ms.timestamps, ms.timestamps(end))./1000;
% dt = abs(diff(temp_ts));
if isfield(ms, 'arena')
[~, speed_epochs] = get_speed_epochs(ms.arena.speed_smooth, params);


ms.speed_epochs = speed_epochs;
is_moving       = speed_epochs;
ms.head_ori     = [];
[ms.room]       = construct_place_maps_2D(ms.room,  ms.room.x(is_moving),  ms.room.y(is_moving),  ms.dt(is_moving), spks(:, is_moving), params.pos_bins, params);
[ms.arena]      = construct_place_maps_2D(ms.arena, ms.arena.x(is_moving), ms.arena.y(is_moving), ms.dt(is_moving), spks(:, is_moving), params.pos_bins, params);
[ms.head_ori]   = construct_place_maps_1D(ms.head_ori,   ms.ori.yaw(is_moving), ms.dt(is_moving), spks(:, is_moving), params.yaw_bin, params);

% %%
% shock_zone_center = pi/2; % typical room shock configuration
% shock_zone_size = pi/6; % size in rad from center to edge
% distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure
% 
% x = ms.room.x; y = ms.room.y; t = ms.timestamps./1000;
% [th, rth] = cart2pol(x,y);
% % figure(1); clf; polarhistogram(th,24);
% 
% d = shock_zone_center - th;
% d = abs(mod(d + pi, 2*pi) - pi);
% 
% 
% % [temp]    = construct_place_maps_1D(ms.room,   th(is_moving), ms.dt(is_moving), spks(:, is_moving), [-pi:pi/6:pi], params);
% [temp]    = construct_place_maps_1D(ms.room,   d(is_moving), ms.dt(is_moving), spks(:, is_moving), [0:pi/32:pi], params);
% qqq = (normalize_rows(temp.pfields_smooth));
% [~, mp] = max(qqq, [], 2);
% [~, ord] = sort(mp);
% qqq = qqq(ord,:);
% 
% s = ms.room.shockTimes;
% shockVec = false(length(ms.timestamps),1);
% e = ms.room.entranceTimes;
% entrVec = false(length(ms.timestamps),1);
% for i = 1:length(s)
%     ind = find(min(abs(ms.timestamps - s(i))) == abs(ms.timestamps - s(i)));
%     shockVec(ind) = true;
% end
% for i = 1:length(e)
%     ind = find(min(abs(ms.timestamps - e(i))) == abs(ms.timestamps - e(i)));
%     entrVec(ind) = true;
% end
% peth_mean = zeros(size(spks,1), 2*22*2 + 1);
% for i = 1:size(spks,1)
% [peth_mean(i,:)] = gb_PETH(spks(i,:), entrVec, 22*2, 22*2);
% % [peth_mean(i,:)] = gb_PETH(spks(i,:), shockVec, 22*2, 22*2);
% end
% [~, mp] = max(peth_mean, [], 2);
% [~, ord] = sort(mp);
% peth_mean = peth_mean(ord,:);

% theta spike maps and distance maps, also PETH of shocks
%%
%     imagesc(ms.ori.pfields_smooth)

% figure(9); clf; hold on
% spks = normalize_rows(ms.neuron.C);
% spks2 = normalize_rows(ms.neuron.S);
% scale = .8;
% for i = 1:size(spks,1)
%     s = spks2(i,:)*scale;
%     plot(s + i - 1, 'k')
%     s = spks(i,:)*scale;
%     plot(s + i - 1, 'r')
% end
[ms.arena.pcell_stats] = place_cell_stats(spks, ms.arena.pfields, ms.arena.spkmap, ms.arena.vmap);
[ms.room.pcell_stats] = place_cell_stats(spks, ms.room.pfields,  ms.room.spkmap,  ms.room.vmap);

%% seting up optimal density subplot

ms.arena.pfield_alpha = ~isnan(ms.arena.vmap);
ms.room.pfield_alpha = ~isnan(ms.room.vmap);
if params.plotting
    nsegs = size(spks,1);
    ss = get(0,'screensize');
    ar = mean([1, ss(3)/ss(4)])-1;
    n = ceil(sqrt(nsegs));
    nc = round(n*(1-ar))+1;
    nr = round(n*(1+ar));
    
    figure(10); clf;
    set(gcf, 'Name', 'ARENA FRAME')
    colormap(plasma)
    figure(11); clf
    set(gcf, 'Name', 'ROOM FRAME')
    colormap(viridis)
    figure(10);
    for i = 1:nsegs
        subplot_tight(nc, nr, i, [.01 .01])
        p = squeeze(ms.arena.pfields_smooth(i,:,:));
        %     p = squeeze(ms.arena.pfields(i,:,:));
        p = p./max(p(:));
        imagesc(p, 'AlphaData', ms.arena.pfield_alpha);
        axis square off
    end
    figure(11);
    for i = 1:nsegs
        subplot_tight(nc, nr, i, [.01 .01])
        p = squeeze(ms.room.pfields_smooth(i,:,:));
        %     p = squeeze(ms.room.pfields(i,:,:));
        p = p./max(p(:));
        imagesc(p, 'AlphaData', ms.room.pfield_alpha);
        axis square off
    end
end
else
    warning('No position data within MS, skipped place maps!')
end % END PLACE CELL PROCESSING
%%
if false
    spks = normalize_rows(ms.neuron.S)>0;
    % spks = normalize_rows(ms.neuron.C+ms.neuron.YrA);
    spks(isnan(spks)) = 0;
    % [pco_tau, pco_prob, sumspks] = Fenton_pco(spks, 11*5, false, 'Kendall');
    samplerate = round(1/median(ms.dt));
    [pco_tau, pco_prob, sumspks] = Fenton_pco(spks, samplerate*10, false, 'Kendall');
    if params.plotting
        figure(4); clf;
        subplot(1,3,1)
        imagesc(normalize_rows(sumspks));
        subplot(1,3,2)
        imagesc(pco_tau, [0 .3]);
        subplot(1,3,3)
        imagesc(pco_prob);
    end
end
end