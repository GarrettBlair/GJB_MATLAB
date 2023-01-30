dataDir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\';
allFiles = {'2022_09_08_H16_52_24_TR0_@placecells.mat';...
'2022_09_10_H15_45_23_TR1_@placecells.mat';...    
'2022_09_10_H16_29_57_TR2_@placecells.mat';...    
'2022_09_12_H15_51_50_TR3_@placecells.mat';...    
'2022_09_12_H16_35_51_TR4_@placecells.mat';...    
'2022_09_13_H17_04_46_TR5_@placecells.mat';...
'2022_09_13_H17_52_09_TR6_@placecells.mat';...
'2022_09_14_H17_07_45_TR7_@placecells.mat';...
'2022_09_14_H17_52_22_WTR8_@placecells.mat';...
'2022_09_16_H17_03_07_TR9_@placecells.mat';...
'2022_09_16_H17_41_09_WTR10_@placecells.mat';...
'2022_09_16_H18_19_23_TR11_@placecells.mat';...
'2022_09_18_H15_22_23_TR12_@placecells.mat';...
'2022_09_18_H16_00_44_TR13_@placecells.mat';...
'2022_09_19_H16_37_23_TR14_@placecells.mat';...
'2022_09_19_H17_07_47_DRK15_@placecells.mat';...
'2022_09_19_H17_39_13_DRK16_@placecells.mat';...
'2022_09_21_H17_49_55_TR17_@placecells.mat';...
'2022_09_21_H18_22_12_DRK18_@placecells.mat';...
'2022_09_24_H16_11_01_WTR19_@placecells.mat';...
'2022_09_24_H16_46_19_TR20_@placecells.mat';...
'2022_09_27_H17_23_57_WTR21_@placecells.mat';...
'2022_09_27_H17_58_42_TR22_@placecells.mat';...
'2022_09_27_H18_26_28_TR23_@placecells.mat';...
'2022_09_28_H17_02_49_TR24_@placecells.mat';...
'2022_09_28_H17_31_07_TR25_@placecells.mat';...
'2022_09_29_H19_23_01_TR99_@placecells.mat';... 
'2022_09_30_H18_18_17_TR26_@placecells.mat';...  
'2022_09_30_H18_50_38_TR27_@placecells.mat'};
for i = 1:length(allFiles)
    allFiles{i} = sprintf('%s%s', dataDir, allFiles{i});
end

HAB_files = {sprintf('%s%s', dataDir, '2022_09_08_H16_09_56_HAB_@placecells.mat')};
% TR_ind = [0 1 2 3 4 5 6 7 8 9 11 12 13 14 17 20 22 23 24 25]+1;
TR_ind =  [7 9  17 22]+1;
WTR_ind = [8 10 19 21]+1;
DRK_ind = [15 16 18]+1;
TR_files = allFiles(TR_ind);
WTR_files = allFiles(WTR_ind);
DRK_files = allFiles(DRK_ind);

% [hab_runs_z] = pfield_dispersion(HAB_files{1}, true);
tr_runs_z = []; wtr_runs_z = []; drk_runs_z = []; 
tr_roompref = []; wtr_roompref = []; drk_roompref = []; 
for i = 1:length(TR_files)
[tr_runs_z(i), tr_roompref(i)] = pfield_dispersion(TR_files{i}, true);
end
for i = 1:length(WTR_files)
[wtr_runs_z(i), wtr_roompref(i)] = pfield_dispersion(WTR_files{i}, true);
end
for i = 1:length(DRK_files)
[drk_runs_z(i), drk_roompref(i)] = pfield_dispersion(DRK_files{i}, true);
end
%% Behavior
animals = {'Hipp18240'};

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

datfolder = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\DAT_files';

TR_nums = TR_ind-1;
[numA, numTR] = size(TR_nums);

TR_entr = NaN(numA, numTR);
TR_shocks = NaN(numA, numTR);
for aLoop = 1:numA
    for sessLoop = 1:numTR
        room_tracking_fname    = sprintf('%s\\%s_TR%d_Room.dat', datfolder, animals{aLoop}, TR_nums(aLoop, sessLoop));
        arena_tracking_fname   = sprintf('%s\\%s_TR%d_Arena.dat', datfolder, animals{aLoop}, TR_nums(aLoop, sessLoop));
        if ~isfile(room_tracking_fname) || ~isfile(arena_tracking_fname)
            error('DAT FILE NOT FOUND')
        end
        [room, arena, params_sub] = behavior_DAT_tracking_eval(room_tracking_fname, arena_tracking_fname, params);
        tot_min = (room.timestamps(end)/1000)/60;
        if isnan(room.first_entrance)
            warning('arena shocks used')
            TR_entr(aLoop, sessLoop) = arena.num_entrances;%/tot_min;
            TR_shocks(aLoop, sessLoop) = arena.num_shocks;%/tot_min;
        else
            TR_entr(aLoop, sessLoop) = room.num_entrances;%/tot_min;
            TR_shocks(aLoop, sessLoop) = room.num_shocks;%/tot_min;
        end
    end
end

WTR_nums = WTR_ind-1;
[numA, numWTR] = size(WTR_nums);

WTR_entr = NaN(numA, numWTR);
WTR_shocks = NaN(numA, numWTR);
for aLoop = 1:numA
    for sessLoop = 1:numWTR
        room_tracking_fname    = sprintf('%s\\%s_WTR%d_Room.dat', datfolder, animals{aLoop}, WTR_nums(aLoop, sessLoop));
        arena_tracking_fname   = sprintf('%s\\%s_WTR%d_Arena.dat', datfolder, animals{aLoop}, WTR_nums(aLoop, sessLoop));
        if ~isfile(room_tracking_fname) || ~isfile(arena_tracking_fname) 
            error('DAT FILE NOT FOUND')
        end
        [room, arena, params_sub] = behavior_DAT_tracking_eval(room_tracking_fname, arena_tracking_fname, params);
        tot_min = (room.timestamps(end)/1000)/60;
        if isnan(room.first_entrance)
            warning('arena shocks used')
            WTR_entr(aLoop, sessLoop) = arena.num_entrances;%/tot_min;
            WTR_shocks(aLoop, sessLoop) = arena.num_shocks;%/tot_min;
        else
            WTR_entr(aLoop, sessLoop) = room.num_entrances;%/tot_min;
            WTR_shocks(aLoop, sessLoop) = room.num_shocks;%/tot_min;
        end
    end
end

DRK_nums = DRK_ind-1;
[numA, numDRK] = size(DRK_nums);

DRK_entr = NaN(numA, numDRK);
DRK_shocks = NaN(numA, numDRK);
for aLoop = 1:numA
    for sessLoop = 1:numDRK
        room_tracking_fname    = sprintf('%s\\%s_DRK%d_Room.dat', datfolder, animals{aLoop}, DRK_nums(aLoop, sessLoop));
        arena_tracking_fname   = sprintf('%s\\%s_DRK%d_Arena.dat', datfolder, animals{aLoop}, DRK_nums(aLoop, sessLoop));
        if ~isfile(room_tracking_fname) || ~isfile(arena_tracking_fname)
            error('DAT FILE NOT FOUND')
        end
        [room, arena, params_sub] = behavior_DAT_tracking_eval(room_tracking_fname, arena_tracking_fname, params);
        tot_min = (room.timestamps(end)/1000)/60;
        if isnan(room.first_entrance)
            warning('arena shocks used')
            DRK_entr(aLoop, sessLoop) = arena.num_entrances;%/tot_min;
            DRK_shocks(aLoop, sessLoop) = arena.num_shocks;%/tot_min;
        else
            DRK_entr(aLoop, sessLoop) = room.num_entrances;%/tot_min;
            DRK_shocks(aLoop, sessLoop) = room.num_shocks;%/tot_min;
        end
    end
end



%%
figure(826); clf; hold on
yyaxis('left')
scatter(1+tr_runs_z*0, tr_runs_z)
scatter(2+wtr_runs_z*0, wtr_runs_z)
scatter(3+drk_runs_z*0, drk_runs_z)
xlim([0 4])
ylim([-6 0])
yyaxis('right')
scatter(1.2+TR_shocks*0, TR_shocks,'k')
scatter(2.2+WTR_shocks*0, WTR_shocks, 'k')
scatter(3.2+DRK_shocks*0, DRK_shocks, 'k')
xlim([0 4])
ylim([0 120])

figure(827); clf; hold on
yyaxis('left')
scatter(1+tr_roompref*0, tr_roompref)
scatter(2+wtr_roompref*0, wtr_roompref)
scatter(3+drk_roompref*0, drk_roompref)
xlim([0 4])
ylim([-.15 .15])
yyaxis('right')
scatter(1.2+TR_shocks*0, TR_shocks,'k')
scatter(2.2+WTR_shocks*0, WTR_shocks, 'k')
scatter(3.2+DRK_shocks*0, DRK_shocks, 'k')
xlim([0 4])
ylim([0 120])
% scatter(TR_ind, tr_roompref)
% scatter(WTR_ind, wtr_roompref)
% scatter(DRK_ind, drk_roompref)

%%
load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\2022_09_16_H17_03_07_TR9_@placecells.mat')
% load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\2022_09_14_H17_52_22_WTR8_@placecells.mat')
craw = normalize_rows(ms.neuron.C+ms.neuron.YrA);
% dff = movmedian(dff, 11, 2);
% dff = conv2(dff, [1 1 1 0 0], 2);
% dff = (diff(ms.neuron.YrA, 1, 2));
% % yra = ms.neuron.YrA;
% % yra =  yra - median(yra,2)*ones(1,size(yra,2));% + std(ms.neuron.YrA,[],2)*2;
% % yra =  yra./(  std(ms.neuron.YrA,[],2)*ones(1,size(yra,2))  );
% % dff(yra<=2) = 0;

% [theta_ang, dist] = cart2pol(ms.room.x, ms.room.y);
%%
[theta_ang, dist] = cart2pol(ms.arena.x, ms.arena.y);
ang_bin = [-pi:pi/8:pi];
[~, ~, thetabin] = histcounts(theta_ang, ang_bin);
nanind = isnan(thetabin);
thetabin(nanind) = interp1(find(~nanind), thetabin(~nanind), find(nanind), 'nearest');
figure(3); clf; plot(thetabin, 'k.-')

%
craw = normalize_rows(ms.neuron.C+ms.neuron.YrA);
[craw_bin] = bin_spks_time(craw, .5, ms.timestamps./1000, false);
craw_bin = normalize_matrix(craw_bin);
[theta_ds_bin] = average_spks_time(thetabin', .5, ms.timestamps./1000, false, 'median');
% a1 = craw(1,:);
% [spmat, spt, spc] = CellsortFindspikes(craw, 2.5, 1/11, 50, true);
% N   = 120;
% Fs  = 11;
% Fst  = 2;
% Fp  = 1.8;
% Ap  = 0.01;
% Ast = 80;
% Rp  = (10^(Ap/20) - 1)/(10^(Ap/20) + 1);
% Rst = 10^(-Ast/20);
% NUM = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge');
% % a = fvtool(NUM,'Fs',Fs)
% LP_FIR = dsp.LowpassFilter('SampleRate',Fs,...
%     'DesignForMinimumOrder',false,'FilterOrder',N,...
%     'PassbandFrequency',Fp,'PassbandRipple',Ap,'StopbandAttenuation',Ast);


% craw = normalize_rows(ms.neuron.C);
[ReducedDataArray, legitimacy_vec] = Ziv_LaplacianEigenVectors(craw, true);
v2=ReducedDataArray{3};

[topo_out] = Ziv_TopoClustering(v2, thetabin(legitimacy_vec));

% options.kappa = 0; 
% options.n_t = 500; % Default: n_t=500. If number of samples per manifold M>500, then set at M. 
% options.flag_NbyM = false;
% [output] = manifold_analysis_outputanchor({ReducedDataArray{3}}, options);
%%
x = ms.arena.x;
y = ms.arena.y;
[theta_ang, dist] = cart2pol(x, y);
ang_bins = [-pi:pi/8:pi];
[~, ~, thetabin] = histcounts(theta_ang, ang_bins);
nanind = isnan(thetabin);
thetabin(nanind) = interp1(find(~nanind), thetabin(~nanind), find(nanind), 'nearest');

dist_bins = [0, 22, 45];
[~, ~, distbin] = histcounts(dist, dist_bins);
nanind = isnan(distbin);
distbin(nanind) = interp1(find(~nanind), distbin(~nanind), find(nanind), 'nearest');
ang_by_dist_bin = thetabin.*distbin;
figure(3); clf; plot(ang_by_dist_bin, 'k.-')

% figure(34); clf; 
% hold on
% cm = viridis(length(ang_bin)-1);
% for i = 1:length(ang_bin)-1
%     scatter3(v2(thetabin==i, 2), v2(thetabin==i, 3), v2(thetabin==i, 4), 20, '.', 'MarkerEdgeColor', cm(i,:));
% end

c_label = topo_out.global.ind_augmented3;
n_cluster = max(c_label);
% cm = jet(n_cluster);
cm = viridis(n_cluster);
[~, ord] = sort(rand(n_cluster,1));
cm = cm(ord,:);
figure(35); clf; 
hold on
for i = 1:n_cluster
    scatter3(v2(c_label==i, 2), v2(c_label==i, 3), v2(c_label==i, 4), 20, '.', 'MarkerEdgeColor', cm(i,:));
end



% figure(36); clf; 
% plot(thetabin,'k')
% hold on
% for i = 1:n_cluster
%     scatter(find(c_label==i), thetabin(c_label==i), 50, '.', 'MarkerEdgeColor', cm(i,:));
% end

figure(37); clf; colormap gray
nsub = ceil(sqrt(n_cluster));
x = x-min(x); x = 20*x./max(x);
y = y-min(y); y = 20*y./max(y);
for i = 1:n_cluster
    subplot_tight(nsub, nsub, i); cla; hold on
%     plot(x, y, 'Color', [.4 .4 .4]);
    imagesc(ms.arena.vmap);
    scatter(x(c_label==i), y(c_label==i), 20, '.', 'MarkerEdgeColor', cm(i,:));
    axis image off
end


x = ms.room.x;
y = ms.room.y;
[theta_ang_room, dist] = cart2pol(x, y);
ang_bin = [-pi:pi/8:pi];
[~, ~, thetabin] = histcounts(theta_ang_room, ang_bin);
nanind = isnan(thetabin);
thetabin(nanind) = interp1(find(~nanind), thetabin(~nanind), find(nanind), 'nearest');

% figure(38); clf; 
% plot(thetabin,'k')
% hold on
% for i = 1:n_cluster
%     scatter(find(c_label==i), thetabin(c_label==i), 50, '.', 'MarkerEdgeColor', cm(i,:));
% end

figure(39); clf; colormap gray
nsub = ceil(sqrt(n_cluster));
x = x-min(x); x = 20*x./max(x);
y = y-min(y); y = 20*y./max(y);
for i = 1:n_cluster
    subplot_tight(nsub, nsub, i); cla; hold on
%     plot(x, y, 'Color', [.4 .4 .4]);
    imagesc(ms.room.vmap);
    scatter(x(c_label==i), y(c_label==i), 20, '.', 'MarkerEdgeColor', cm(i,:));
    axis image off
end


% figure(40); clf; hold on
% plot(c_label)








