function [ms] = APA_generate_placmaps(ms, params)
caimanFilename = sprintf('%s/MiniLFOV/caiman_cnmfe_out.mat', ms.parentDir);
caiman_data = load(caimanFilename);
% correct for python indexing start at 0
caiman_data.idx_components_bad = caiman_data.idx_components_bad+1;
caiman_data.idx_components = caiman_data.idx_components+1;

%% Cull contours
[good_inds, bad_inds, ~, bounds] = Draw_contour_bounding(caiman_data.fullA, ...
    caiman_data.dims, caiman_data.maxFrame, caiman_data.idx_components, false);
allbad = unique([caiman_data.idx_components_bad, bad_inds']);
fprintf('\nRemoving %d bad components\n', length(allbad))
neuron = remove_segments(caiman_data, allbad, false);

% neuron = caiman_data;
ms.neuron = neuron;
ms.contour_bounds = bounds;
spks = normalize_rows(ms.neuron.S);
% temp_ts = cat(1, ms.timestamps, ms.timestamps(end))./1000;
% dt = abs(diff(temp_ts));
[partition_ID, speed_epochs] = get_speed_epochs(ms.arena.speed_smooth, params);

% ms.is_moving    = speed_epochs;
ms.speed_epochs    = speed_epochs;
is_moving       = speed_epochs;
ms.head_ori = [];
[ms.room]   = construct_place_maps_2D(ms.room,  ms.room.x(is_moving),  ms.room.y(is_moving),  ms.dt(is_moving), spks(:, is_moving), params.pos_bins, params);
[ms.arena]  = construct_place_maps_2D(ms.arena, ms.arena.x(is_moving), ms.arena.y(is_moving), ms.dt(is_moving), spks(:, is_moving), params.pos_bins, params);
[ms.head_ori]    = construct_place_maps_1D(ms.head_ori,   ms.ori.yaw(is_moving), ms.dt(is_moving), spks(:, is_moving), params.yaw_bin, params);
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
[a] = place_cell_stats(spks, ms.arena.pfields, ms.arena.spkmap, ms.arena.vmap);
[r] = place_cell_stats(spks, ms.room.pfields,  ms.room.spkmap,  ms.room.vmap);

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