function [ms] = APA_generate_placemaps(ms, params)

spks = ms.spks; %normalize_rows(ms.neuron.S_matw);

% temp_ts = cat(1, ms.timestamps, ms.timestamps(end))./1000;
% dt = abs(diff(temp_ts));
if isfield(ms, 'arena')
[~, speed_epochs] = get_speed_epochs(ms.arena.speed_smooth, params);

[nsegs, nframes] = size(ms.neuron.C);

ms.speed_epochs = speed_epochs;
is_moving       = speed_epochs;
ms.is_moving       = speed_epochs;
ms.head_ori     = [];

ms.room.split_vec = [];
ms.arena.split_vec = [];
[ms.room]       = construct_place_maps_2D(ms.room,  ms.room.x(is_moving),  ms.room.y(is_moving),  ms.dt(is_moving), spks(:, is_moving), params.pos_bins, params.pos_bins, params);
[ms.arena]      = construct_place_maps_2D(ms.arena, ms.arena.x(is_moving), ms.arena.y(is_moving), ms.dt(is_moving), spks(:, is_moving), params.pos_bins, params.pos_bins, params);

full_splits = NaN(nframes,1);
full_splits(find(is_moving)) = ms.room.pfield_split_vec;
full_splits(find(~is_moving)) = uint8(interp1(find(is_moving), double(ms.room.pfield_split_vec), find(~is_moving), 'nearest'));
ms.room.split_vec = full_splits;
full_splits = NaN(nframes,1);
full_splits(find(is_moving)) = ms.arena.pfield_split_vec;
full_splits(find(~is_moving)) = uint8(interp1(find(is_moving), double(ms.arena.pfield_split_vec), find(~is_moving), 'nearest'));
ms.arena.split_vec = full_splits;

[ms.room.pfield_decode]  = pfield_split_bayesian_decoding(ms, ms.room.x,  ms.room.y,  spks, params.pos_bins, params.pos_bins, 1:size(spks,1), params.num_random_shuffle_decode, params.ipos_int_time);
[ms.arena.pfield_decode] = pfield_split_bayesian_decoding(ms, ms.arena.x, ms.arena.y, spks, params.pos_bins, params.pos_bins, 1:size(spks,1), params.num_random_shuffle_decode, params.ipos_int_time);


[ms.head_ori]   = construct_place_maps_1D(ms.head_ori,   ms.ori.yaw(is_moving), ms.dt(is_moving), spks(:, is_moving), params.yaw_bins, params);

[ms.arena.pcell_stats] = place_cell_stats(spks, ms.arena.pfields_smooth, ms.arena.spkmap_smooth, ms.arena.vmap);
[ms.room.pcell_stats]  = place_cell_stats(spks, ms.room.pfields_smooth,  ms.room.spkmap_smooth,  ms.room.vmap);

%%
if isfield(params, 'num_random_shuffle_pcell')
    nrand = params.num_random_shuffle_pcell;
    ms.arena.pcell_stats.infoPerSpike_rand = NaN(nsegs, nrand);
    ms.arena.pcell_stats.sparsity_rand     = NaN(nsegs, nrand);
    ms.arena.pcell_stats.splitcorr_rand    = NaN(nsegs, nrand);
    ms.room.pcell_stats.infoPerSpike_rand  = NaN(nsegs, nrand);
    ms.room.pcell_stats.sparsity_rand      = NaN(nsegs, nrand);
    ms.room.pcell_stats.splitcorr_rand     = NaN(nsegs, nrand);
    a_infoPerSpike_rand = NaN(nsegs, nrand);
    a_sparsity_rand     = NaN(nsegs, nrand);
    a_coh_rand          = NaN(nsegs, nrand);
    a_corr_rand         = NaN(nsegs, nrand);
    r_infoPerSpike_rand = NaN(nsegs, nrand);
    r_sparsity_rand     = NaN(nsegs, nrand);
    r_coh_rand          = NaN(nsegs, nrand);
    r_corr_rand         = NaN(nsegs, nrand);
    
    fprintf('\t\tPeforming %d shuffles for mutual info, sparsity, & split corr... ', nrand)
    
    roomxm = ms.room.x(is_moving);    roomym = ms.room.y(is_moving);
    arenaxm = ms.arena.x(is_moving);  arenaym = ms.arena.y(is_moving);
    dt = ms.dt(is_moving);
    
    minshift = floor(length(roomxm)*.1);
    maxshift = floor(length(roomxm)*.9);
    shiftval = randi([minshift maxshift], nrand, 4);
    randdir = 1*(rand(nrand, 4)>=.5);
    randdir(randdir==0) = -1;
    shiftval = shiftval.*randdir;
    spks_mov = spks(:, is_moving);
    pos_bins = params.pos_bins;
    h = tic;
    reps = linspace(0, nrand, 11);
    parfor_progbar = params.parfor_progbar;
    ppm = 0;
    if parfor_progbar == true
        ppm = ParforProgressbar(nrand,'parpool', {'local'}, 'showWorkerProgress',true,...
            'progressBarUpdatePeriod',5,'title','Pfield info randLoop');
    end
    parfor randLoop = 1:nrand
        rx = circshift(roomxm,  shiftval(randLoop,1));
        ry = circshift(roomym,  shiftval(randLoop,2));
        ax = circshift(arenaxm, shiftval(randLoop,3));
        ay = circshift(arenaym, shiftval(randLoop,4));
        room_temp = [];
        arena_temp = [];
        [room_temp]                         = construct_place_maps_2D(room_temp,  rx,  ry,  dt, spks_mov, pos_bins, pos_bins, params);
        [room_rand_pcell_stats]             = place_cell_stats(spks, room_temp.pfields_smooth,  room_temp.spkmap_smooth,  room_temp.vmap);
        r_infoPerSpike_rand(:, randLoop)    = room_rand_pcell_stats.infoPerSpike;
        r_sparsity_rand(:, randLoop)        = room_rand_pcell_stats.sparsity;
        r_coh_rand(:, randLoop)             = room_rand_pcell_stats.coherence;
        r_corr_rand(:, randLoop)            = room_temp.split_corr;
        
        [arena_temp]                        = construct_place_maps_2D(arena_temp, ax,  ay,  dt, spks_mov, pos_bins, pos_bins, params);
        [arena_rand_pcell_stats]            = place_cell_stats(spks, arena_temp.pfields_smooth, arena_temp.spkmap_smooth, arena_temp.vmap);
        a_infoPerSpike_rand(:, randLoop)    = arena_rand_pcell_stats.infoPerSpike;
        a_sparsity_rand(:, randLoop)        = arena_rand_pcell_stats.sparsity;
        a_coh_rand(:, randLoop)             = arena_rand_pcell_stats.coherence;
        a_corr_rand(:, randLoop)            = arena_temp.split_corr;
        if parfor_progbar == true
            pause(.001)
            ppm.increment();
        end
    end
    fprintf(' Done! %.2f seconds\n', toc(h))
    ms.arena.pcell_stats.infoPerSpike_rand = a_infoPerSpike_rand;
    ms.room.pcell_stats.infoPerSpike_rand  = r_infoPerSpike_rand;
    ms.arena.pcell_stats.sparsity_rand = a_sparsity_rand;
    ms.room.pcell_stats.sparsity_rand  = r_sparsity_rand;
    ms.arena.pcell_stats.coherence_rand = a_coh_rand;
    ms.room.pcell_stats.coherence_rand  = r_coh_rand;
    ms.arena.pcell_stats.splitcorr_rand = a_corr_rand;
    ms.room.pcell_stats.splitcorr_rand  = r_corr_rand;
    
    temp = ms.arena.pcell_stats.infoPerSpike*ones(1, nrand);
    ms.arena.pcell_stats.infoProb = sum(temp<ms.arena.pcell_stats.infoPerSpike_rand, 2)./nrand;
    temp = ms.arena.pcell_stats.sparsity*ones(1, nrand);
    ms.arena.pcell_stats.sparsityProb = sum(temp<ms.arena.pcell_stats.sparsity_rand, 2)./nrand;
    temp = ms.arena.split_corr*ones(1, nrand);
    ms.arena.pcell_stats.splitcorrProb = sum(temp<ms.arena.pcell_stats.splitcorr_rand, 2)./nrand;
    temp = ms.arena.pcell_stats.coherence*ones(1, nrand);
    ms.arena.pcell_stats.coherenceProb = sum(temp<ms.arena.pcell_stats.coherence_rand, 2)./nrand;
    
    temp = ms.room.pcell_stats.infoPerSpike*ones(1, nrand);
    ms.room.pcell_stats.infoProb = sum(temp<ms.room.pcell_stats.infoPerSpike_rand, 2)./nrand;
    temp = ms.room.pcell_stats.sparsity*ones(1, nrand);
    ms.room.pcell_stats.sparsityProb = sum(temp<ms.room.pcell_stats.sparsity_rand, 2)./nrand;
    temp = ms.room.split_corr*ones(1, nrand);
    ms.room.pcell_stats.splitcorrProb = sum(temp<ms.room.pcell_stats.splitcorr_rand, 2)./nrand;
    temp = ms.room.pcell_stats.coherence*ones(1, nrand);
    ms.room.pcell_stats.coherenceProb = sum(temp<ms.room.pcell_stats.coherence_rand, 2)./nrand;
end

%% seting up optimal density subplot

ms.arena.pfield_alpha = ~isnan(ms.arena.vmap);
ms.room.pfield_alpha = ~isnan(ms.room.vmap);
if false % params.plotting
    nsegs = size(spks,1);
    ss = get(0,'screensize');
    ar = mean([1, ss(3)/ss(4)])-1;
    n = ceil(sqrt(nsegs));
    nc = round(n*(1-ar))+1;
    nr = round(n*(1+ar));
    
%     figure(10); clf;
%     set(gcf, 'Name', 'ARENA FRAME')
%     colormap(plasma)
%     figure(11); clf
%     set(gcf, 'Name', 'ROOM FRAME')
    nfigs = ceil(nsegs/400);
    for figloop=1:nfigs
    subsegs = 400*(figloop-1)+1:400*(figloop);
    subsegs = subsegs(subsegs<=nsegs);
    for i = 1:length(subsegs)
    figure(10+figloop);
        subplot_tight(20, 20, i, [.01 .01])
        p = squeeze(ms.arena.pfields_smooth(subsegs(i),:,:));
        %     p = squeeze(ms.arena.pfields(i,:,:));
        pr = p./max(p(:));
        imagesc(p, 'AlphaData', ms.arena.pfield_alpha);
        axis square off
%     end
%     for i = 1:nsegs
    figure(100+figloop);
        subplot_tight(20, 20, i, [.01 .01])
        p = squeeze(ms.room.pfields_smooth(subsegs(i),:,:));
        %     p = squeeze(ms.room.pfields(i,:,:));
        pa = p./max(p(:));
        imagesc([pr pa], 'AlphaData', [ ms.room.pfield_alpha  ms.arena.pfield_alpha]);
        axis square off
    end
    colormap(viridis)
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