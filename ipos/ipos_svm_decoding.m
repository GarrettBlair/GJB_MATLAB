function ipos_svm_decoding(ipos_file, dataSaveName, figSaveName, x_fold_training)
% load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\ipos_files\Acc20832_TR7_int_0.25sec.mat');
% load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\2023_01_13_H19_26_00_TR7_@placecells.mat');
% load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\ipos_files\Acc19947_TR6_int_0.25sec.mat') 
% load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\2023_01_13_H14_38_46_TR6_@placecells.mat');
%%
load(ipos_file, 'ms_room_temp', 'room_momentary_pos', 'ms_arena_temp', 'arena_momentary_pos',...
    'room_ensemble_prob', 'arena_ensemble_prob') 
% x_fold_training = 5;
% vars2save = {'ipos_file' 'x_fold_training' 'room_angle*' 'room_pos*'...
%     'arena_angle*' 'arena_pos*', 'dataSaveName' , 'figSaveName'};

%% ROOM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shock_zone_center = pi/2; % typical room shock configuration
% shock_zone_size = pi/6; % size in rad from center to edge
% distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure
rng(1)

ref_string = {'room' 'arena'};
% ref_string = {'room'};
% decode_var = {'pos'}; is2D = [true];
decode_var = {'angle'}; is2D = [false];
% decode_var = {'angle', 'pos'}; is2D = [false, true];

% var2use = {'ipos_mean' 'ipos_ensemble'};
% var2use = {'spks' 'ipos' };
var2use = {'spks'};
% var2use = {'spks' 'ipos' 'ipos_mean' 'ipos_ensemble'};
% var2use = {'ipos_mean' 'ipos_ensemble'};

vars2save = {'ipos_file' 'x_fold_training', 'dataSaveName' , 'figSaveName'};
for i = 1:length(ref_string)
    for j = 1:length(decode_var)
        vars2save = cat(2, vars2save, [ref_string{i} '_' decode_var{j} '*']);
    end
end



pos_bins = [-45, -36:4:36, 45];
occ_thresh = 1; % minimum occupancy in time (sec) to be valid fo mapping error

spks = 1*ms_room_temp.spks_bin;%ms.neuron.S_matw>0;
% ipos = 1*ipos_ra;
ipos = (abs(room_momentary_pos) - abs(arena_momentary_pos));
ipos_mean = nanmean(abs(room_momentary_pos) - abs(arena_momentary_pos));
ipos_ensemble = abs(room_ensemble_prob) - abs(arena_ensemble_prob);
ipos_ensemble = ipos_ensemble';
for refLoop = 1:length(ref_string)
    eval(sprintf('ref_struct = ms_%s_temp;', ref_string{refLoop})) 
    x = ref_struct.x; y = ref_struct.y; occ_map = ref_struct.vmap; t = ms_room_temp.t;
    [pos_angle, rth] = cart2pol(x,y);
    for decodeVarLoop = 1:length(decode_var)
        strName = [ref_string{refLoop} '_' decode_var{decodeVarLoop}];
        %         eval(sprintf('%s = [];', strName))
        temp = [];
%         temp.abins          = -pi:pi/12:pi; % range for binning angular position
%         temp.rbins          = 0:5:45; % range fo binning polar distance
        temp.abins          = -pi:pi/24:pi; % range for binning angular position
        temp.rbins          = 0:5:45; % range fo binning polar distance
        temp.abin_center    = temp.abins(2:end) - abs(diff(temp.abins))/2;
        temp.rbin_center    = temp.rbins(2:end) - abs(diff(temp.rbins))/2;
        temp.theta          = pos_angle;
        temp.r              = rth;
        for varLoop = 1:length(var2use)
            %%
            eval(sprintf('predVar = %s;', var2use{varLoop}))
            fprintf('Decoding  [%s_%s]  xfold=%d\n', strName, var2use{varLoop}, x_fold_training)
            room_bias  = ipos_mean>0;
            arena_bias  = ipos_mean<=0;
            if is2D(decodeVarLoop)
                varStruct = svm_decode_sub(temp, {pos_angle, rth}, {temp.abins temp.rbins}, predVar, x_fold_training, 'eucl_dist_real');
                varStruct2 = svm_decode_sub_iposbias(varStruct, ipos_mean, {pos_angle, rth}, {temp.abins temp.rbins}, predVar, x_fold_training, 'eucl_dist_real');
%                 rStruct = svm_decode_sub(temp, {pos_angle(room_bias), rth(room_bias)}, {temp.abins temp.rbins}, predVar(:, room_bias), x_fold_training, 'polar_dist_real');
%                 aStruct = svm_decode_sub(temp, {pos_angle(arena_bias), rth(arena_bias)}, {temp.abins temp.rbins}, predVar(:, arena_bias), x_fold_training, 'polar_dist_real');
            else
                varStruct = svm_decode_sub(temp, pos_angle, temp.abins, predVar, x_fold_training, 'polar_dist_real');
                varStruct2 = svm_decode_sub_iposbias(varStruct, ipos_mean, pos_angle, temp.abins, predVar, x_fold_training, 'polar_dist_real');
%                 aStruct = svm_decode_sub(temp, pos_angle(arena_bias), temp.abins, predVar(:, arena_bias), x_fold_training, 'polar_dist_real');
            end     
            %%%%% CHANGE DIFF TO ERR
            varStruct2.predictor     = predVar;
            varStruct2.predictor_var = var2use{varLoop};
            varStruct2.ipos_mean     = ipos_mean;
            
            [sum_all, counts_all]       = make_occupancymap_2D(x,  y,  varStruct.pred_err, pos_bins, pos_bins);
            av_err_map = sum_all./counts_all; av_err_map(occ_map <= occ_thresh) = NaN;
            varStruct2.error_map     = av_err_map;
            
            % coherent reference
            [sum_r, counts_r]       = make_occupancymap_2D(x,  y,  varStruct2.pred_err_same, pos_bins, pos_bins);
            av_err_map = sum_r./counts_r; av_err_map(occ_map <= occ_thresh) = NaN;
            varStruct2.error_map_same    = av_err_map;
            % incoherent reference
            [sum_a, counts_a]       = make_occupancymap_2D(x,  y,  varStruct2.pred_err_diff, pos_bins, pos_bins);
            av_err_map = sum_a./counts_a; av_err_map(occ_map <= occ_thresh) = NaN;
            varStruct2.error_map_diff    = av_err_map;

            figure; 
            set(gcf, 'Name', sprintf('%s_%s', strName, var2use{varLoop}))
            subplot(2,2,1)
            imagesc(varStruct2.error_map); axis square; colorbar
            subplot(2,2,2)
            imagesc(varStruct2.error_map_same); axis square; colorbar
            subplot(2,2,3)
            imagesc(varStruct2.error_map_diff); axis square; colorbar
            subplot(2,2,4)
            imagesc(occ_map); axis square; colorbar
            drawnow
%             eval(sprintf('%s = temp;', strName))
%             eval(sprintf('%s_%s_single = varStruct;', strName, var2use{varLoop}))
            eval(sprintf('%s_%s = varStruct2;', strName, var2use{varLoop}))
        end
    end
end

%%
% [room_angle_iposr]  = svm_decode_sub(room_angle, room_pos_angle, room_angle.abins, spks_ipos_room,   x_fold_training, 'polar_dist_real');
% [room_angle_iposa]  = svm_decode_sub(room_angle, room_pos_angle, room_angle.abins, spks_ipos_arena,  x_fold_training, 'polar_dist_real');

% [room_pos_spks]   = svm_decode_sub(room_angle, {room_pos_angle, rth}, {room_angle.abins room_angle.rbins}, spks,           x_fold_training, 'polar_dist_real');
% [room_pos_ipos]   = svm_decode_sub(room_angle, {room_pos_angle, rth}, {room_angle.abins room_angle.rbins}, spks_ipos,      x_fold_training, 'polar_dist_real');
% [room_pos_iposr]  = svm_decode_sub(room_angle, {room_pos_angle, rth}, {room_angle.abins room_angle.rbins}, spks_ipos_room, x_fold_training, 'polar_dist_real');
% [room_pos_iposa]  = svm_decode_sub(room_angle, {room_pos_angle, rth}, {room_angle.abins room_angle.rbins}, spks_ipos_arena, x_fold_training, 'polar_dist_real');

% [room_dist_spks]   = svm_decode_sub(room_dist, sz_dist, room_dist.bins, spks,             x_fold_training, 'eucl_dist_real');
% [room_dist_ipos]   = svm_decode_sub(room_dist, sz_dist, room_dist.bins, spks_ipos,        x_fold_training, 'eucl_dist_real');
% [room_dist_iposr]  = svm_decode_sub(room_dist, sz_dist, room_dist.bins, spks_ipos_room,   x_fold_training, 'eucl_dist_real');
% [room_dist_iposa]  = svm_decode_sub(room_dist, sz_dist, room_dist.bins, spks_ipos_arena,  x_fold_training, 'eucl_dist_real');
%%
if false
fig = figure;
set(fig, 'Position', [216 351 931 454], 'Name', ipos_file)
clrstyle = [.1 .1 .1; .5 .2 .8; .8 .2 .2; .2 .2 .7];

subplot(2,2,1); cla; hold on;
if exist('room_angle_spks', 'var')
room_angle_means = summary_figs(room_angle_spks, room_angle_ipos, room_angle_ipos_mean, room_angle_ipos_ensemble, clrstyle);
title('Room Frame Angle'); set(gca, 'XTickLabel', {'spks' 'ipos' 'ipos_mean' 'ipos_ensemble'},'XTickLabelRotation', 12)
end
% subplot(2,1,2); hold on;
% dist_means = summary_figs(room_dist_spks, room_dist_ipos, room_dist_iposr, room_dist_iposa);
% title('Room Frame Shock Dist.'); set(gca, 'XTickLabel', {'spks' 'ipos(room)' 'ipos(arena)' 'ipos(room-arena)'})
subplot(2,2,3); cla; hold on;
if exist('room_pos_spks', 'var')
room_pos_means = summary_figs(room_pos_spks, room_pos_ipos, room_pos_ipos_mean, room_pos_ipos_ensemble, clrstyle);
ylim([-5 85])
title('Room Frame Position'); set(gca, 'XTickLabel', {'spks' 'ipos' 'ipos_mean' 'ipos_ensemble'},'XTickLabelRotation', 12)
end

subplot(2,2,2); cla; hold on;
if exist('arena_angle_spks', 'var')
arena_angle_means = summary_figs(arena_angle_spks, arena_angle_ipos, arena_angle_ipos_mean, arena_angle_ipos_ensemble, clrstyle);
ylim([-.2 3.5])
title('Arena Frame Angle'); set(gca, 'XTickLabel', {'spks' 'ipos' 'ipos_mean' 'ipos_ensemble'},'XTickLabelRotation', 12)
end
% subplot(2,1,2); hold on;
% dist_means = summary_figs(arena_dist_spks, arena_dist_ipos, arena_dist_iposr, arena_dist_iposa);
% title('Room Frame Shock Dist.'); set(gca, 'XTickLabel', {'spks' 'ipos(room)' 'ipos(arena)' 'ipos(room-arena)'})

subplot(2,2,4); cla; hold on;
if exist('arena_pos_spks', 'var')
arena_pos_means = summary_figs(arena_pos_spks, arena_pos_ipos, arena_pos_ipos_mean, arena_pos_ipos_ensemble, clrstyle);
ylim([-5 85])
title('Arena Frame Position'); set(gca, 'XTickLabel', {'spks' 'ipos' 'ipos_mean' 'ipos_ensemble'},'XTickLabelRotation', 12)
end
drawnow
saveas(fig, figSaveName)
end
save(dataSaveName, vars2save{:});

% figure;
% set(gcf, 'Name', sprintf('%s   -   %s_%s', ipos_file, strName, var2use{varLoop}))
% subplot(2,2,1)
% imagesc(arena_pos_ipos_ensemble.error_map)
% subplot(2,2,2)
% imagesc(arena_pos_ipos_ensemble.error_mapA)
% subplot(2,2,3)
% imagesc(arena_pos_ipos_ensemble.error_mapR)
% subplot(2,2,4)
% imagesc(occ_map)
% drawnow


end
%%
function means = summary_figs(room_angle_spks, room_angle_ipos, room_angle_iposr, room_angle_iposa, clrstyle)
angdiff         = room_angle_spks.pred_err';
angdiff_r       = room_angle_spks.rand_err';
angdiff_ip      = room_angle_ipos.pred_err';
angdiff_r_ip    = room_angle_ipos.rand_err';
angdiff_ipr     = room_angle_iposr.pred_err';
angdiff_r_ipr   = room_angle_iposr.rand_err';
angdiff_ipa     = room_angle_iposa.pred_err';
angdiff_r_ipa   = room_angle_iposa.rand_err';
means = [];
means.diff_spks             = mean(angdiff);
means.diff_ipos             = mean(angdiff_ip);
means.diff_iposroom         = mean(angdiff_ipr);
means.diff_iposarena        = mean(angdiff_ipa);
means.rand_err_spks        = mean(angdiff_r);
means.rand_err_ipos        = mean(angdiff_r_ip);
means.rand_err_iposroom    = mean(angdiff_r_ipr);
means.rand_err_iposarena   = mean(angdiff_r_ipa);

% figure; hold on;
% plot(t, pos_angle)
% plot(t, room_angle_iposr.bin_center(room_angle_iposr.pred_real))
% plot(t, room_angle_spks.bin_center(room_angle_ipos.pred_real))

% figure(8); clf;
% hold on
% v = violinplot([angdiff angdiff_ip angdiff_ipr angdiff_ipa], [1 2 3 4]);%, 'Bandwidth', pi/12) ;
cla
v = violinplot([angdiff angdiff_ip angdiff_ipr angdiff_ipa], [1 2 3 4]);%, 'Bandwidth', pi/12) ;
for i = 1:4
v(i).EdgeColor = clrstyle(i,:)*0;
v(i).ShowData = 0;
v(i).ShowMean = 0;
v(i).ShowNotches = 0;
v(i).ViolinColor = clrstyle(i,:);
v(i).ViolinAlpha = .3;
rectangle('Position', [i-.45 0 .425 max(max([angdiff angdiff_ip angdiff_ipr angdiff_ipa]))], 'FaceColor', 'w', 'EdgeColor', 'w')
end
rms = .35;
plot([1+rms 1-rms], [nanmedian(angdiff_r)      nanmedian(angdiff_r)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([2+rms 2-rms], [nanmedian(angdiff_r_ip)     nanmedian(angdiff_r_ip)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([3+rms 3-rms], [nanmedian(angdiff_r_ipr)       nanmedian(angdiff_r_ipr)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([4+rms 4-rms], [nanmedian(angdiff_r_ipa)    nanmedian(angdiff_r_ipa)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
boxplot([angdiff angdiff_ip angdiff_ipr angdiff_ipa], 'Notch', 'off', 'Widths', .2,...
    'plotstyle', 'traditional', 'Jitter', .35, 'Symbol', '', 'Color', clrstyle, 'Positions', [1 2 3 4]-.2)
axis([0.25 4.75 -.2 3.45 ])


drawnow
end
%%

% end