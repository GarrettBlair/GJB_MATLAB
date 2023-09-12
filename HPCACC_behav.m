params = [];
if isempty(params)
    params.arena_radius             = 40; % in cm
    params.arena_center             = [ 127.5, 127.5]; % pixel center of behav cam, [x,y]
    params.pixpercm                 = 3.1220; % pixel radius of behav cam env
    params.behav_fps                = 30;
    params.behav_smoothing_interval = .25; % in seconds, length of smoothing kernel
    
    params.pos_bins                 = [-45, -36:4:36, 45]; % in cm, x and y
    params.yaw_bin                  = -pi:pi/8:pi;
    params.occupancy_thresh         = 1; % in seconds, minimum time in each bin to be counted for place map
    params.pfield_kernel_radius     = 3; % kernel ends up being [n*2 + 1] in bins
    params.smin_vals                = -50:5:-10; % smin values used to create the 'deconv_sweep.mat'
    % params.speed_thresh             = 5; % speed thresh in cm/sec
    params.num_partitions           = 2;
    % params.max_spd_thresh           = 100;
    % params.min_spd_thresh           = 5;
    params.max_spd_thresh           = 100;
    params.min_spd_thresh           = 0;
    params.min_samples              = 10;
    params.ipos_int_time             = .25;%.2; % seconds, binning time for computing momentary spatial information
    
    % params.rotate_behav             = true;
    params.nan_interp               = true;
    params.remove_bad_caiman_segs   = false;
    params.correct_dt               = true; % correct for large jumps in timestamp file when constructing vmap
    params.plotting                 = false;
end
%%
ddir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA';
animals = {'HPCACC24500' 'HPCACC24501' 'HPCACC24502' 'HPCACC24503'};
trials = {'HAB0', 'TR0',  'TR1',  'TR2',  'TR3',  'TR4', 'TR6',  'TR7',  'TR8',  'TR9', 'TR10',  'TR11'...
    'TR12', 'TR13', 'TR14', 'TR15', 'TR16', };
% no TR5
na = length(animals);
nt = length(trials);

numShocks      = NaN(nt, na);
numEntr        = NaN(nt, na);
occ_diff_room  = NaN(nt, na);
occ_diff_arena = NaN(nt, na);
sessMinutes    = NaN(nt, na);

for anum = 1:na
    for trloop = 1:nt
        roomfname  = [ddir '\DAT_files\' animals{anum} '_' trials{trloop} '_Room.dat'];
        arenafname = [ddir '\DAT_files\' animals{anum} '_' trials{trloop} '_Arena.dat'];
        if isfile(roomfname)
            warning('off', 'MATLAB:structOnObject')
            [room, arena, params_sub] = behavior_DAT_tracking_eval(roomfname, arenafname, params);
            sessMinutes( trloop, anum ) = (room.timestamps(end) - room.timestamps(1))/(60000);
            numShocks( trloop, anum ) = room.num_shocks;
            numEntr( trloop, anum ) = room.num_entrances;
            hr = histcounts(room.pol_theta,  [-pi:pi/16:pi], 'Normalization', 'probability');
            ha = histcounts(arena.pol_theta, [-pi:pi/16:pi], 'Normalization', 'probability');
            hr = sort(hr);
            ha = sort(ha);
            null_h = ones(size(hr))./length(hr);
%             figure; hold on; plot(ha); plot(hr)
%             plot(abs(hr-null_h))
    %         occ_diff( trloop, anum ) = sum(abs(ha-hr));
            occ_diff_room( trloop, anum ) = sum(abs(hr-null_h));
            occ_diff_arena( trloop, anum ) = sum(abs(ha-null_h));
        end
    end
end
%%
figure(1023); clf
set(gcf, 'Position', [128 646 1045 297])
subplot(1,3,1);
plot(numEntr./sessMinutes)
set(gca, 'XTick', [1:nt], 'XTickLabel', trials, 'XTickLabelRotation', -90)
subplot(1,3,2);
plot(occ_diff_room)
set(gca, 'XTick', [1:nt], 'XTickLabel', trials, 'XTickLabelRotation', -90)
subplot(1,3,3);
plot(occ_diff_arena)
legend(animals)   
set(gca, 'XTick', [1:nt], 'XTickLabel', trials, 'XTickLabelRotation', -90)
% save('D:\Sample Data\HPCACC_behav.mat', 'animals', 'trials', 'numShocks', 'numEntr', 'params')