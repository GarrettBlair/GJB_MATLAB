function [ms, room, arena, params_sub] = behavior_DAT_tracking_file(ms, room_tracking_fname, arena_tracking_fname, params)
if isempty(params)
APA_rat_imaging_params_current
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
    elseif ms.timestamps(end) < room.timestamps(end)
        diff_ts = find(room.timestamps > ms.timestamps(end), 1)-1;
        fprintf('~~~MS software ended before TRACKER software ended, \n\t .dat timestamps cuttof at %d index\n', diff_ts);
        fields = fieldnames(room);
        length2cut = length(room.timestamps);
        for jj = 1:length(fields)
            temp = eval(sprintf('room.%s;', fields{jj}));
            temp_info = whos('temp');
            if any(temp_info.size == length2cut) && ~strcmp(temp_info.class, 'char')
                temp = temp(1:diff_ts);
                eval(sprintf('room.%s = temp;', fields{jj}));
            end
        end
        room.entrance_start_idx = room.entrance_start_idx(room.entrance_start_idx<= length(room.timestamps));
        room.shock_start_idx    = room.shock_start_idx(room.shock_start_idx<= length(room.timestamps));
        ms.warnings.TrackerTS_diffs = ...
            sprintf('MS software ended before TRACKER software ended, .dat timestamps cuttof at %d index', diff_ts);
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
if isfield(room, 'Angle')
ms.room.headangle   = interp1(room.timestamps, room.Angle, ms.timestamps, 'linear', 'extrap');
ms.arena.headangle  = interp1(arena.timestamps, arena.Angle, ms.timestamps, 'linear', 'extrap');
end

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
    figure; clf;
    subplot(2,2,1)
    hold on;
%     scatter3(room.x, room.y, room.timestamps./1000, 1, 'k.')
%     scatter3(ms.room.x, ms.room.y, ms.timestamps./1000, 1, 'r.')
    plot(ms.room.x, ms.room.y, 'r-')
    title('Extracted Room')
%     set(gca, 'View', [-50 50])
    axis square
    
    subplot(2,2,2)
    hold on;
%     scatter3(arena.x, arena.y, arena.timestamps./1000, 1, 'k.')
    plot(ms.arena.x, ms.arena.y, 'b-')
    title('Extracted Arena')
%     set(gca, 'View', [-50 50])
    axis square
    
    subplot(2,2,3:4)
    hold on;
    plot(ms.timestamps./1000, ms.arena.speed, 'k', 'LineWidth', 2)
    
    plot(ms.timestamps(ms.is_moving)./1000, ms.arena.speed_smooth(ms.is_moving), 'g.')
    title('Arena Speed')
    axis tight
end
end