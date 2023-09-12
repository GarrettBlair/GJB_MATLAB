%%
% datfolder = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\Stress_platform\DAT_files';
datfolder = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\OCAM imaging\DAT_files';
animals = {'Hipp18239' 'Hipp18240' 'Hipp17861' 'Hipp17862' 'Hipp17863' 'Hipp17864'};
is_stressed = [0 1 0 1 0 1]==1;
TR_nums = [0:3; 24:27; 24:27; 0:3; 0:3; 0:3];

animals = {'Hipp18239' 'Hipp18240' 'Hipp17862' 'Hipp17863' 'Hipp17864'};
is_stressed = [0 1 1 0 1]==1;
TR_nums = [0:3; 24:27; 0:3; 0:3; 0:3];

params = [];
params.arena_radius             = 40; % in cm
params.arena_center             = [ 127.5, 127.5]; % pixel center of behav cam, [x,y]
params.pixpercm                 = 3.1220; % pixel radius of behav cam env
params.behav_fps                = 30;
params.behav_smoothing_interval = .25; % in seconds, length of smoothing kernel

params.pos_bins                 = -40:4:40; % in cm, x and y
params.yaw_bin                  = -pi:pi/8:pi;
params.occupancy_thresh         = 1; % in seconds, minimum time in each bin to be counted for place map
params.nan_interp               = true;
params.correct_dt               = false; % correct for large jumps in timestamp file when constructing vmap
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
%


[numA, numTR] = size(TR_nums);

num_entr = NaN(numA, numTR);
num_shocks = NaN(numA, numTR);
total_time = NaN(numA, numTR);
for aLoop = 1:numA
    for sessLoop = 1:numTR
        room_tracking_fname    = sprintf('%s\\%s_TR%d_Room.dat', datfolder, animals{aLoop}, TR_nums(aLoop, sessLoop));
        arena_tracking_fname   = sprintf('%s\\%s_TR%d_Arena.dat', datfolder, animals{aLoop}, TR_nums(aLoop, sessLoop));
        if ~isfile(room_tracking_fname) || ~isfile(arena_tracking_fname) 
            error('DAT FILE NOT FOUND')
        end
        [room, arena, params_sub] = behavior_DAT_tracking_eval(room_tracking_fname, arena_tracking_fname, params);
        tot_min = (room.timestamps(end)/1000)/60;
        num_entr(aLoop, sessLoop) = room.num_entrances;%/tot_min;
        num_shocks(aLoop, sessLoop) = room.num_shocks;%/tot_min;
        total_time(aLoop, sessLoop) = tot_min;
        figure(8);
        subplot(numA, numTR, numTR*(aLoop-1) + sessLoop)
        plot3(room.x, room.y, room.timestamps./1000)
        figure(9);
        subplot(numA, numTR, numTR*(aLoop-1) + sessLoop)
        room_vmap = histcounts(room.pol_theta, [-pi:pi/12:pi]);
%         room_vmap = make_occupancymap_2D(room.x, room.y, room.timestamps./1000, params.pos_bins, params.pos_bins);
        room_vmap = room_vmap-min(room_vmap(:));
        room_vmap = room_vmap./max(room_vmap(:));
        plot(room_vmap)
        
    end
end
%%
% num_entr_tot    = [num_entr(:, 1)+num_entr(:, 2), num_entr(:, 3)+num_entr(:, 4)];
% num_shocks_tot  = [num_shocks(:, 1)+num_shocks(:, 2), num_shocks(:, 3)+num_shocks(:, 4)];
num_entr_tot    = [(num_entr(:, 1)+num_entr(:, 2))./(total_time(:, 1)+total_time(:, 2)),...
    (num_entr(:, 3)+num_entr(:, 4))./(total_time(:, 3)+total_time(:, 4))];
num_shocks_tot    = [(num_shocks(:, 1)+num_shocks(:, 2))./(total_time(:, 1)+total_time(:, 2)),...
    (num_shocks(:, 3)+num_shocks(:, 4))./(total_time(:, 3)+total_time(:, 4))];
% num_shocks_tot  = [num_shocks_min(:, 1)+num_shocks_min(:, 2), num_shocks_min(:, 3)+num_shocks_min(:, 4)];
stress_ent_mean = nanmean(num_entr_tot(is_stressed, :), 1);
cont_ent_mean   = nanmean(num_entr_tot(~is_stressed, :), 1);
stress_shk_mean = nanmean(num_shocks_tot(is_stressed, :), 1);
cont_shk_mean   = nanmean(num_shocks_tot(~is_stressed, :), 1);

figure(963); clf; 
% plotx = [1 2 4 5];
plotx = [1 3];
d1 = plotx(1)-.5;
d2 = plotx(end)+.5;
stresscolor = [.5 .1 .1];
controlcolor = [.3 .3 .3];

subplot(121); hold on
title('Entrances')
for j = 1:numA
    if is_stressed(j)==1
        c = stresscolor;
    else
        c = controlcolor;
    end
    plot(plotx, num_entr_tot(j,:), '-', 'Color', c, 'LineWidth', 1.25); 
end
plot(plotx, cont_ent_mean, 'Color', controlcolor.*2, 'LineWidth', 3); 
plot(plotx, stress_ent_mean, 'Color', stresscolor.*2, 'LineWidth', 3); 
ylabel('Number of entrances/min')
xlabel('Session')
% axis([d1 d2 -10 125])
axis([d1 d2 -1 5])
% set(gca, 'XTick', plotx, 'YTick', [0:20:100])
set(gca, 'XTick', plotx, 'YTick', [0:1:5], 'XTickLabel', {'pre' 'post'})

subplot(122); hold on
title('# of Shocks')
for j = 1:numA
    if is_stressed(j)==1
        c = stresscolor;
    else
        c = controlcolor;
    end
    plot(plotx, num_shocks_tot(j,:), '-', 'Color', c, 'LineWidth', 1.25); 
end
plot(plotx, cont_shk_mean, 'Color', controlcolor.*2, 'LineWidth', 3); 
plot(plotx, stress_shk_mean, 'Color', stresscolor.*2, 'LineWidth', 3); 
% axis([d1 d2 -10 280])
axis([d1 d2 -1 10])
ylabel('Number of shocks/min')
xlabel('Session')
% set(gca, 'XTick', plotx, 'YTick', [0:50:250])
set(gca, 'XTick', plotx, 'YTick', [0:2:10], 'XTickLabel', {'pre' 'post'})
%%
names = {'animal', 'stressed', 'Day1Entr', 'Day3Entr', 'Day1Shk', 'Day3Shk'};
T = table(animals', is_stressed', num_entr_tot(:,1), num_entr_tot(:,2), num_shocks_tot(:,1), num_shocks_tot(:,2), 'VariableNames', names(:));
fname = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\Stress_platform\BehavEffects.csv';
% writetable(T, fname)
%%
matchingfile = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\Stress_platform\Hipp18240\matching_contours\matching_matrix.mat';
filenames = ...
   {'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\Stress_platform\Hipp18240\processed_files\2022_09_28_H17_02_49_TR24_@placecells.mat';...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\Stress_platform\Hipp18240\processed_files\2022_09_28_H17_31_07_TR25_@placecells.mat';...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\Stress_platform\Hipp18240\processed_files\2022_09_29_H19_23_01_TR99_@placecells.mat';...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\Stress_platform\Hipp18240\processed_files\2022_09_30_H18_18_17_TR26_@placecells.mat';...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\Stress_platform\Hipp18240\processed_files\2022_09_30_H18_50_38_TR27_@placecells.mat'};
nsess = 5;
pop_sec_res = 60;
cell_sec_res = 1;

load(matchingfile, 'cellmap')
sess_corr = NaN(nsess);
sess_absdiff = NaN(nsess);

num_shared = [cellmap'>0] * [cellmap>0];
min_shared = min(num_shared(:));
numcells2use = min_shared;%[];%
for s1 = 1:nsess-1
    for s2 = s1+1:nsess
        %%
        temp = load(filenames{s1}, 'ms');
        ms1 = temp.ms;
        temp = load(filenames{s2}, 'ms');
        ms2 = temp.ms;
        
        
        matched = cellmap(:,s1)>0 & cellmap(:,s2)>0;
        cells_s1 = cellmap(matched,s1);
        cells_s2 = cellmap(matched,s2);
        spks1 = ms1.neuron.S_matw(cells_s1, :);
        spks2 = ms2.neuron.S_matw(cells_s2, :);
        
        spks1 = normalize_rows(spks1);
        spks2 = normalize_rows(spks2);
        nsamples = 6943;%%%%%%%%%%%%%% min(size(spks1,2), size(spks2,2)); % minimum samles in stress recording
        spks1 = spks1(:, 1:nsamples);
        spks2 = spks2(:, 1:nsamples);
        
        pop_dt1 = round(1/median(ms1.dt))*pop_sec_res;
        pop_dt2 = round(1/median(ms2.dt))*pop_sec_res;
        pcell_dt1 = round(1/median(ms1.dt))*cell_sec_res;
        pcell_dt2 = round(1/median(ms2.dt))*cell_sec_res;
        
        if ~isempty(numcells2use)
            [~, randord] = sort(rand(sum(matched), 1));
            spks1 = spks1(randord(1:numcells2use), :);
            spks2 = spks2(randord(1:numcells2use), :);
        end
        [timecorr1, ~, ~] = Fenton_pop_stability(spks1, pop_sec_res, ms1.timestamps(1:nsamples)./1000, false);
        [timecorr2, ~, ~] = Fenton_pop_stability(spks2, pop_sec_res, ms2.timestamps(1:nsamples)./1000, false);
        joint_t = ms1.timestamps(1:nsamples)./1000;
        joint_t = [joint_t; joint_t(end) + cell_sec_res + ms2.timestamps(1:nsamples)./1000];
        [timecorr3, ~, ~] = Fenton_pop_stability([spks1, spks2], pop_sec_res, joint_t, false);
        %             [popcorr3, ~, cellcorr3, ~, ~, ~] = Fenton_pop_stability(cat(2, spks1, spks2), pop_dt1, false);
        
        [cellcorr1, ~, ~] = Fenton_cell_corr(spks1, cell_sec_res, ms1.timestamps(1:nsamples)./1000, false);
        [cellcorr2, ~, ~] = Fenton_cell_corr(spks2, cell_sec_res, ms2.timestamps(1:nsamples)./1000, false);
        
        c1 = triu(cellcorr1,1);
        c2 = triu(cellcorr2,1);
        lowerinds = c1==0 & c2==0;
        corrs1 = c1(~lowerinds);
        corrs2 = c2(~lowerinds);
        sess_corr(s1,s2) = corr(corrs1, corrs2);
        sess_corr(s2,s1) = sess_corr(s1,s2);
        
        pw_diff = corrs1 - corrs2;
        sess_absdiff(s1,s2) = mean(abs(pw_diff));
        sess_absdiff(s2,s1) = sess_absdiff(s1,s2);
        
        [corrs1_sort, ord] = sort(corrs1, 'descend');
        corrs2_sort = corrs2(ord);
        tops = floor(length(ord)/4);
        corrs2_sort = corrs2_sort(1:tops);
        corrs1_sort = corrs1_sort(1:tops);
    end
end
figure; subplot(121); imagesc(sess_corr, [-.2 1]); subplot(122); imagesc(sess_absdiff, [0 .15]);
figure; subplot(131); imagesc(timecorr1, [-.2 .6]); subplot(132); imagesc(timecorr2, [-.2 .6]); subplot(133); imagesc(timecorr3, [-.2 .6]);
figure; subplot(131); imagesc(cellcorr1, [-.3 .3]); subplot(132); imagesc(cellcorr2, [-.3 .3]); subplot(133); imagesc(cellcorr1-cellcorr2, [-.3 .3]);
figure; hold on; plot(corrs1_sort); plot(corrs2_sort)



















