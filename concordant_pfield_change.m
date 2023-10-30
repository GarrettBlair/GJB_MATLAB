function concordant_pfield_change(ms, params)
%%
% clear

% load('D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_07_06_H17_16_26_CON22_@placecells_HPC_miniscope1.mat')
% load('D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_07_06_H17_16_26_CON22_@placecells_ACC_miniscope2.mat')
% load('D:\GarrettBlair\APA\Hipp18240\processed_files\2022_09_16_H17_03_07_TR9_@placecells.mat')
params.occupancy_thresh = 0;

% spks = ms.spks>0;
% [spks_bin, group] = bin_spks_time(spks, params.ipos_int_time, ms.timestamps./1000, false);
spks_bin = ms.room.svm_decoding.spks_bin;
rx = ms.room.svm_decoding.x;
ry = ms.room.svm_decoding.y;
dt = ms.room.svm_decoding.dt;
ax = ms.arena.svm_decoding.x;
ay = ms.arena.svm_decoding.y;


%%
sig2 = nanmean(abs(ms.room.momentary_pos_info) - abs(ms.arena.momentary_pos_info));
% room_surprise = -1*sig2';%abs(ms.room.conjoint_ipos_av) - abs(ms.arena.conjoint_ipos_av);
% room_surprise = abs(ms.room.conjoint_ipos_av) - abs(ms.arena.conjoint_ipos_av);
room_surprise = ms.arena.conjoint_ipos_av>ms.room.conjoint_ipos_av;
% room_inds = find( room_surprise<0 ); % greater 'surprise' in arena frame==room more likely
% arena_inds = find( room_surprise>0 ); % greater 'surprise' in room frame==arena more likely
room_inds = find(ms.arena.conjoint_ipos_av>ms.room.conjoint_ipos_av); % greater 'surprise' in arena frame==room more likely
arena_inds = find(ms.arena.conjoint_ipos_av<ms.room.conjoint_ipos_av); % greater 'surprise' in room frame==arena more likely

fitmodel= false;
if fitmodel==true
r_cone = ms.room.conjoint_ipos_av;
a_cone = ms.arena.conjoint_ipos_av;
vmin = min(min([r_cone, a_cone]));
r_cone = r_cone - vmin;
a_cone = a_cone - vmin;
vmax = max(max([r_cone, a_cone]));
r_cone = r_cone./vmax;
a_cone = a_cone./vmax;
% r_cone = normalize_matrix(r_cone);
% a_cone = normalize_matrix(a_cone);
goods = ~isnan(r_cone) & ~isnan(a_cone);
r_cone = r_cone(goods) + .001;
a_cone = a_cone(goods) + .001;

ebins = [0:0.001:1.001];
[r_model] = gb_conjoint_ensemble_model_fitting(r_cone, nanmedian(r_cone)+2*nanstd(r_cone), ebins, true);
[a_model] = gb_conjoint_ensemble_model_fitting(a_cone, nanmedian(a_cone)+2*nanstd(a_cone), ebins, true);
end
% room_surprise = abs(ms.room.conjoint_ipos_min) - abs(ms.arena.conjoint_ipos_min);
% sig1 = ms.room.conjoint_ipos_av - ms.arena.conjoint_ipos_av;
% room_inds = find( ms.arena.conjoint_ipos_av > ms.room.conjoint_ipos_av ); % greater 'surprise' in arena frame==room more likely
% arena_inds = find( ms.room.conjoint_ipos_av > ms.arena.conjoint_ipos_av ); % greater 'surprise' in room frame==arena more likely

nsamples = min([length(room_inds), length(arena_inds)]);
[~, ords] = sort(rand(nsamples, 2));
[~, ord] = sort(rand(length(room_surprise), 1));
room_inds  = room_inds(ords(:,1));
arena_inds = arena_inds(ords(:,1));

randord = [room_inds(1:2:end); arena_inds(1:2:end)];


ms_roomsig = [];
% ms_roomsig = ms;
[ms_roomsig.room]       = construct_place_maps_2D(ms_roomsig,  rx(room_inds),  ry(room_inds),  dt(room_inds), spks_bin(:, room_inds), params.pos_bins, params.pos_bins, params);
[ms_roomsig.arena]      = construct_place_maps_2D(ms_roomsig, ax(room_inds), ay(room_inds), dt(room_inds), spks_bin(:, room_inds), params.pos_bins, params.pos_bins, params);
[ms_roomsig.null]       = construct_place_maps_2D(ms_roomsig,  rx(randord),  ry(randord),  dt(randord), spks_bin(:, randord), params.pos_bins, params.pos_bins, params);
% [ms_roomsig.null]       = construct_place_maps_2D(ms_roomsig,  rx,  ry,  dt, spks_bin(:, :), params.pos_bins, params.pos_bins, params);

ms_arenasig = [];
% ms_arenasig = ms;
[ms_arenasig.room]       = construct_place_maps_2D(ms_arenasig,  rx(arena_inds),  ry(arena_inds),  dt(arena_inds), spks_bin(:, arena_inds), params.pos_bins, params.pos_bins, params);
[ms_arenasig.arena]      = construct_place_maps_2D(ms_arenasig, ax(arena_inds), ay(arena_inds), dt(arena_inds), spks_bin(:, arena_inds), params.pos_bins, params.pos_bins, params);
[ms_arenasig.null]       = construct_place_maps_2D(ms_arenasig,  ax(randord),  ay(randord),  dt(randord), spks_bin(:, randord), params.pos_bins, params.pos_bins, params);
% [ms_arenasig.null]       = construct_place_maps_2D(ms_arenasig,  ax,  ay,  dt, spks_bin(:, :), params.pos_bins, params.pos_bins, params);
%%
figure(944); clf; 
scatter(ms_arenasig.arena.split_corr, ms_roomsig.room.split_corr, 'o', 'MarkerFaceAlpha', .2, 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'none')
figure(943); clf; 
scatter(ms_roomsig.arena.split_corr, ms_arenasig.room.split_corr, 'o', 'MarkerFaceAlpha', .2, 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'none')
figure(942); clf; 
scatter(ms_arenasig.arena.split_corr - ms_roomsig.arena.split_corr, ms_roomsig.room.split_corr - ms_arenasig.room.split_corr, 'o', 'MarkerFaceAlpha', .2, 'MarkerFaceColor', 'm', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'none')
axis equal
% scatter(ms_arenasig.arena.pcell_stats.infoPerSpike, ms_roomsig.room.pcell_stats.infoPerSpike, 'o', 'MarkerFaceAlpha', .2, 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'none')


figure(11); clf; bins = [-.5:.025:1];
subplot(2,1,1); hold on; 
histogram(ms_roomsig.room.split_corr, bins)
histogram(ms_arenasig.room.split_corr, bins)
% histogram(ms_roomsig.null.split_corr, bins)
subplot(2,1,2); hold on
histogram(ms_arenasig.arena.split_corr, bins)
histogram(ms_roomsig.arena.split_corr, bins)
% histogram(ms_arenasig.null.split_corr, bins)

room_pfield_corr_change = ms_roomsig.room.split_corr - ms_arenasig.room.split_corr;
arena_pfield_corr_change = ms_arenasig.arena.split_corr - ms_roomsig.arena.split_corr;

figure(12); clf; bins = [-1.5:.05:1.5];
subplot(3,1,1); hold on; 
histogram(ms_roomsig.room.split_corr - ms_arenasig.room.split_corr, bins)
subplot(3,1,2); hold on
histogram(ms_arenasig.arena.split_corr - ms_roomsig.arena.split_corr, bins)
subplot(3,1,3); hold on
histogram([arena_pfield_corr_change], bins)
histogram([room_pfield_corr_change], bins)

figure(945); clf; 
scatter(arena_pfield_corr_change, room_pfield_corr_change, 'o', 'MarkerFaceAlpha', .2, 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'none')

[ms_roomsig.arena.pcell_stats] = place_cell_stats(spks_bin, ms_roomsig.arena.pfields_smooth, ms_roomsig.arena.spkmap_smooth, ms_roomsig.arena.vmap);
[ms_roomsig.room.pcell_stats]  = place_cell_stats(spks_bin, ms_roomsig.room.pfields_smooth,  ms_roomsig.room.spkmap_smooth,  ms_roomsig.room.vmap);
[ms_arenasig.arena.pcell_stats] = place_cell_stats(spks_bin, ms_arenasig.arena.pfields_smooth, ms_arenasig.arena.spkmap_smooth, ms_arenasig.arena.vmap);
[ms_arenasig.room.pcell_stats]  = place_cell_stats(spks_bin, ms_arenasig.room.pfields_smooth,  ms_arenasig.room.spkmap_smooth,  ms_arenasig.room.vmap);

%%
figure(13); clf; bins = [0:.025:2];
subplot(2,1,1); hold on; 
histogram(ms_roomsig.room.pcell_stats.infoPerSpike, bins)
histogram(ms_arenasig.room.pcell_stats.infoPerSpike, bins)
% histogram(ms_roomsig.null.split_corr, bins)
subplot(2,1,2); hold on
histogram(ms_arenasig.arena.pcell_stats.infoPerSpike, bins)
histogram(ms_roomsig.arena.pcell_stats.infoPerSpike, bins)
% histogram(ms_arenasig.null.split_corr, bins)

room_pfield_info_change = ms_roomsig.room.pcell_stats.infoPerSpike - ms_arenasig.room.pcell_stats.infoPerSpike;
arena_pfield_info_change = ms_arenasig.arena.pcell_stats.infoPerSpike - ms_roomsig.arena.pcell_stats.infoPerSpike;

figure(14); clf; ;
subplot(3,1,1); hold on; 
% histogram(ms_roomsig.room.pcell_stats.infoPerSpike - ms_arenasig.room.pcell_stats.infoPerSpike, bins)
histogram2(ms_roomsig.room.pcell_stats.infoPerSpike, ms_arenasig.room.pcell_stats.infoPerSpike, bins, bins)
subplot(3,1,2); hold on
histogram2(ms_roomsig.arena.pcell_stats.infoPerSpike, ms_arenasig.arena.pcell_stats.infoPerSpike, bins, bins)
% histogram(ms_arenasig.arena.pcell_stats.infoPerSpike - ms_roomsig.arena.pcell_stats.infoPerSpike, bins)
subplot(3,1,3); hold on; bins = [-1:.05:1];
histogram([arena_pfield_info_change], bins)
histogram([room_pfield_info_change], bins)
% hc = histcounts2(arena_pfield_info_change, room_pfield_info_change, bins, bins);
% pcolor(log2(hc))
figure(945); clf; 
scatter(arena_pfield_info_change, room_pfield_info_change, 'o', 'MarkerFaceAlpha', .2, 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'none')
% scatter(arena_pfield_info_change, room_pfield_info_change)
%%



% rdr = nanmean(ms.room.svm_decoding.pred_err(room_inds)) / nanmean(ms.room.svm_decoding.rand_dist_ex(room_inds))
% rda = nanmean(ms.room.svm_decoding.pred_err(arena_inds)) / nanmean(ms.room.svm_decoding.rand_dist_ex(arena_inds))
% ada = nanmean(ms.arena.svm_decoding.pred_err(arena_inds)) / nanmean(ms.arena.svm_decoding.rand_dist_ex(arena_inds))
% adr = nanmean(ms.arena.svm_decoding.pred_err(room_inds)) / nanmean(ms.arena.svm_decoding.rand_dist_ex(room_inds))
rdr = nanmean(ms.room.svm_decoding.pred_err(room_inds)) / nanmean(ms.room.svm_decoding.rand_err_median)
rda = nanmean(ms.room.svm_decoding.pred_err(arena_inds)) / nanmean(ms.room.svm_decoding.rand_err_median)
ada = nanmean(ms.arena.svm_decoding.pred_err(arena_inds)) / nanmean(ms.arena.svm_decoding.rand_err_median)
adr = nanmean(ms.arena.svm_decoding.pred_err(room_inds)) / nanmean(ms.arena.svm_decoding.rand_err_median)


rdr = nanmean(ms.room.pfield_decode.decode_dist(room_inds)) / nanmean(ms.room.pfield_decode.decode_dist_shuffle_median)
rda = nanmean(ms.room.pfield_decode.decode_dist(arena_inds)) / nanmean(ms.room.pfield_decode.decode_dist_shuffle_median)
ada = nanmean(ms.arena.pfield_decode.decode_dist(arena_inds)) / nanmean(ms.room.pfield_decode.decode_dist_shuffle_median)
adr = nanmean(ms.arena.pfield_decode.decode_dist(room_inds)) / nanmean(ms.room.pfield_decode.decode_dist_shuffle_median)
drawnow







%%














