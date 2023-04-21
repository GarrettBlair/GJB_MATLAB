function [ms, room, arena, params_sub] = behavior_DAT_tracking_file(ms, room_tracking_fname, arena_tracking_fname, params)
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

params_sub = params;
ksize = round(params_sub.behav_fps*params_sub.behav_smoothing_interval);
kern = ones(ksize, 1); kern = kern./sum(kern(:));

[room, arena, params_sub] = behavior_DAT_tracking_eval(room_tracking_fname, arena_tracking_fname, params_sub);
%%
if abs( ms.timestamps(end) - room.timestamps(end) )> 1000  || abs( ms.timestamps(end) - arena.timestamps(end) )> 1000
    warning('!~!~! Large timestamp discrepancy found between ms and dat files!')
    disp([ms.fileName])
%     return
    if ms.timestamps(end) > room.timestamps(end)
        diff_ts = find(ms.timestamps > room.timestamps(end), 1)-1;
        fprintf('~~~Tracker crashed before ms software ended, \n\ttimestamps cuttof at %d index\n', diff_ts);
        ms.timestamps       = ms.timestamps(1:diff_ts);
        ms.frameNum         = ms.frameNum(1:diff_ts);
        ms.dt               = ms.dt(1:diff_ts);
        ms.dt_corrected     = ms.dt_corrected(1:diff_ts);
        ms.ori              = ms.ori(1:diff_ts, :);
        ms.warnings.TrackerCrash = ...
            sprintf('Tracker crashed before ms software ended, timestamps cuttof at %d index', diff_ts);
    else
        error('Uknown solution')
    end
end

ms.room.x = interp1(room.timestamps, room.x, ms.timestamps, 'linear', 'extrap');
ms.room.y = interp1(room.timestamps, room.y, ms.timestamps, 'linear', 'extrap');
ms.arena.x = interp1(arena.timestamps, arena.x, ms.timestamps, 'linear', 'extrap');
ms.arena.y = interp1(arena.timestamps, arena.y, ms.timestamps, 'linear', 'extrap');
ms.room.speed   = interp1(room.timestamps, room.speed, ms.timestamps, 'linear', 'extrap');
ms.arena.speed  = interp1(arena.timestamps, arena.speed, ms.timestamps, 'linear', 'extrap');

% Extract the shock times, entrances, and approaches
ms.room.entranceTimes = room.timestamps(room.entrance_start_idx);
ms.room.shockTimes = room.timestamps(room.shock_start_idx);
ms.arena.entranceTimes = arena.timestamps(arena.entrance_start_idx);
ms.arena.shockTimes = arena.timestamps(arena.shock_start_idx);

% % dt = ms.dt_corrected;
% % ms.room.speed  = sqrt(diff([ms.room.x(1); ms.room.x]).^2   + diff([ms.room.y(1); ms.room.y]).^2)./dt;
% % ms.arena.speed = sqrt(diff([ms.arena.x(1); ms.arena.x]).^2 + diff([ms.arena.y(1); ms.arena.y]).^2)./dt;
ms.room.speed_smooth  = conv(ms.room.speed,  kern, 'same');
ms.arena.speed_smooth = conv(ms.arena.speed, kern, 'same');
is_moving = ms.arena.speed_smooth>params_sub.min_spd_thresh;

[~, speed_epochs] = get_speed_epochs(ms.arena.speed_smooth, params_sub);
ms.speed_epochs = speed_epochs;
ms.is_moving = is_moving;


if params_sub.plotting
    figure(3); clf;
    subplot(2,2,1)
    hold on;
    scatter3(behav.arena_x, behav.arena_y, ts, 1, 'k.')
    scatter3(behav.x, behav.y, ts, 1, 'b.')
    title('Extracted')
    set(gca, 'View', [-60 60])
    axis square
    
    subplot(2,2,2)
    hold on;
    plot3(ax, ay, ts, 'k-')
    plot3(x, y, ts, 'b-')
    title('Interp & Scaled')
    set(gca, 'View', [-60 60])
    axis square
    
    subplot(2,2,3)
    hold on;
    plot3(x, y, ts, 'b-')
    plot3(xx, yy, ts, 'r-')
    title('Room (b) vs. Arena (r)')
    set(gca, 'View', [-60 60])
    axis square
    
    subplot(4,2,6)
    % hold on;
    plot(ts, behav.spd_roomframe, 'b-')
    title('Room speed (b)')
    axis tight
    subplot(4,2,8)
    plot(ts, behav.spd_arenaframe, 'r-')
    title('Arena speed (r)')
    % set(gca, 'View', [-60 60])
    axis tight
end
end