function [SR, SL, LR, LL, zones, turn_arounds] = linear_segmentation_ScaryMaze_LFOV( x, y, speed, dt, always_check_return, plotting)
% %% Draw what points are valid
% 
% left.bounds = ...
%     [-120.01 -30;
%     -120.01 15;
%     -200 15;
%     -200 -30]; % valid points to include in left goal, Rectanglular shock maze
% right.bounds = ...
%     [120.01 -30;
%     120.01 15;
%     200 15;
%     200 -30]; % valid points to include in right goal, Rectanglular shock maze
% short.bounds  = ...
%     [-120 -30;
%     -120 15;
%     120 15;
%     120 -30]; % valid points to include in long path, Rectanglular shock maze
% long.bounds  = ...
%     [-200   15.01;
%     200   15.01;
%     200   200;
%     -200   200]; % valid points to include in short path, Rectanglular shock maze
% 
% 
% right.inside = inpolygon(x, y, right.bounds(:,1), right.bounds(:,2));
% left.inside = inpolygon(x, y, left.bounds(:,1), left.bounds(:,2));
% short.inside  = inpolygon(x, y,  short.bounds(:,1),  short.bounds(:,2));
% long.inside  = inpolygon(x, y,  long.bounds(:,1),  long.bounds(:,2));
[short, long, left, right] = ScaryMaze_zones(x,y);

right.inside_indx    = find(right.inside);
left.inside_indx    = find(left.inside);
short.inside_indx     = find(short.inside);
long.inside_indx     = find(long.inside);
%%
if isempty(long.inside_indx)
    long.inside_indx(1) = NaN;
end
if isempty(short.inside_indx)
    short.inside_indx(1) = NaN;
end
if isempty(left.inside_indx)
    left.inside_indx(1) = NaN;
end
if isempty(right.inside_indx)
    right.inside_indx(1) = NaN;
end

seg_i = 0;
% seg_start   = min([left.inside_indx(1) right.inside_indx(1)]);
seg_start   = min([left.inside_indx(1) right.inside_indx(1) short.inside_indx(1) long.inside_indx(1)]);
if right.inside(seg_start) == 1
    on_rightward_journey = false; % journey is from right to left
elseif left.inside(seg_start) == 1
    on_rightward_journey = true; % journey is from left to right
elseif short.inside(seg_start) == 1
    if (left.inside_indx(1) < right.inside_indx(1)) || isnan(right.inside_indx(1))
        on_rightward_journey = true; % journey is from left to right
    elseif (left.inside_indx(1) > right.inside_indx(1)) || isnan(left.inside_indx(1))
        on_rightward_journey = false; % journey is from left to right
    end
    seg_end = min([left.inside_indx(1), right.inside_indx(1)]);
    seg_i = seg_i+1;
    all_seg_dir_rightwards(seg_i) = on_rightward_journey;
    all_seg_start(seg_i) = seg_start;
    all_seg_end(seg_i) = seg_end;
    
    seg_start = seg_end + find((short.inside(seg_end:end) | long.inside(seg_end:end)), 1);
elseif long.inside(seg_start) == 1
    if (left.inside_indx(1) < right.inside_indx(1)) || isnan(right.inside_indx(1))
        on_rightward_journey = true; % journey is from left to right
    elseif (left.inside_indx(1) > right.inside_indx(1)) || isnan(left.inside_indx(1))
        on_rightward_journey = false; % journey is from left to right
    end
    seg_end = min([left.inside_indx(1), right.inside_indx(1)]);
    seg_i = seg_i+1;
    all_seg_dir_rightwards(seg_i) = on_rightward_journey;
    all_seg_start(seg_i) = seg_start;
    all_seg_end(seg_i) = seg_end;

    seg_start = seg_end + find((short.inside(seg_end:end) | long.inside(seg_end:end)), 1);
end
if plotting
    figure; 
    subplot(121)
    hold on
    plot(x,y)
    plot(left.bounds(:,1), left.bounds(:,2), ':', 'LineWidth', 2)
    plot(right.bounds(:,1), right.bounds(:,2), ':', 'LineWidth', 2)
    plot(long.bounds(:,1), long.bounds(:,2), 'LineWidth', 2)
    plot(short.bounds(:,1), short.bounds(:,2), 'LineWidth', 2)
end

%% While loop to step through data

path_template.index         = false(size(x));
path_template.master_ind    = zeros(size(x));
path_template.pass_ind      = NaN(1,2);
path_template.pass_time     = NaN(1,1);
path_template.pass_length   = NaN(1,1);
path_template.min_spd       = NaN(1,1);
path_template.mean_spd      = NaN(1,1);


SR = path_template;
SL = path_template;
LR = path_template;
LL = path_template;
zones.right_goal    = right;
zones.left_goal     = left;
zones.short_path    = short;
zones.long_path     = long;


%
ending      = false;

min_turnaround_length = 3/dt;
turn_around_counter = 0;
turn_start = NaN(100,1);
turn_end = NaN(100,1);
turn_peak = NaN(100,1);
turn_type = NaN(100,1);
turn_key = {'SR' 'LR' 'SL' 'LL'};
if ~isnan(seg_start)
    while ending == false
        switch on_rightward_journey % if rightwards, goal is right zone, return is left zone
            case false
                if always_check_return
                    check_return = true;
                    while check_return == true
                        next_return = seg_start + find(right.inside(seg_start:end), 1);
                        next_goal = seg_start + find(left.inside(seg_start:end), 1);
                        if next_goal < next_return
                            check_return = false;
                        elseif isempty(next_return)
                            check_return = false;
                        else
                            % Find turn arounds and split segment into two
                            tempseg = seg_start:next_return;
                            if length(tempseg)>=min_turnaround_length
                                turn_around_counter = turn_around_counter + 1;
                                maxy_change = max(y(tempseg) - y(seg_start));
                                maxx_change = max(x(tempseg) - x(seg_start));
                                if maxy_change > maxx_change
                                    turn_index = seg_start + find(y(tempseg)-y(seg_start) == maxy_change, 1, 'last')-1;
                                else
                                    turn_index = seg_start + find(x(tempseg)-x(seg_start) == maxx_change, 1, 'last')-1;
                                end
                                turn_start(turn_around_counter) = seg_start;
                                
                                turn_peak(turn_around_counter) = turn_index;
                                if short.inside(turn_index) == 1
                                    turn_type(turn_around_counter) = 1;
                                elseif long.inside(turn_index) == 1
                                    turn_type(turn_around_counter) = 2;
                                end
                                turn_end(turn_around_counter) = next_return;
                            end
                            seg_start = next_return + min([find(short.inside(next_return:end), 1) find(long.inside(next_return:end), 1)])-1;
                        end
                    end
                    seg_end = seg_start + find((left.inside(seg_start:end)), 1) - 1;
                else
                    next_return = seg_start + find(right.inside(seg_start:end), 1);
                    next_goal = seg_start + find(left.inside(seg_start:end), 1);
                    seg_end = min([next_return next_goal]) - 1;
                    %
                    next_path = min([find(short.inside(seg_end:end), 1, 'first') find(long.inside(seg_end:end), 1, 'first')]) - 1;
                    seg_end = seg_end + next_path;
                end
            case true
                if always_check_return
                    check_return = true;
                    while check_return == true
                        next_return = seg_start + find(left.inside(seg_start:end), 1);
                        next_goal = seg_start + find(right.inside(seg_start:end), 1);
                        if next_goal < next_return
                            check_return = false;
                        elseif isempty(next_return)
                            check_return = false;
                        else
                            % Find turn arounds and split segment into two
                            tempseg = seg_start:next_return;
                            if length(tempseg)>=min_turnaround_length
                                turn_around_counter = turn_around_counter + 1;
                                maxy_change = max(y(tempseg) - y(seg_start));
                                maxx_change = max(x(tempseg) - x(seg_start));
                                if maxy_change > maxx_change
                                    turn_index = seg_start + find(y(tempseg)-y(seg_start) == maxy_change, 1, 'last')-1;
                                else
                                    turn_index = seg_start + find(x(tempseg)-x(seg_start) == maxx_change, 1, 'last')-1;
                                end
                                turn_start(turn_around_counter) = seg_start;
                                turn_end(turn_around_counter) = next_return;
                                turn_peak(turn_around_counter) = turn_index;
                                if short.inside(turn_index) == 1
                                    turn_type(turn_around_counter) = 1;
                                elseif long.inside(turn_index) == 1
                                    turn_type(turn_around_counter) = 2;
                                end
                            end
                            seg_start = next_return + min([find(short.inside(next_return:end), 1) find(long.inside(next_return:end), 1)])-1;
                        end
                    end
                    seg_end = seg_start + find((right.inside(seg_start:end)), 1) - 1;
                else
                    next_return = seg_start + find(left.inside(seg_start:end), 1);
                    next_goal = seg_start + find(right.inside(seg_start:end), 1);
                    seg_end = min([next_return next_goal]) - 1;
                    %
                    next_path = min([find(short.inside(seg_end:end), 1, 'first') find(long.inside(seg_end:end), 1, 'first')]) - 1;
                    seg_end = seg_end + next_path;
                end
        end
        
        if ~isempty(seg_end)
            seg_i = seg_i+1;
            all_seg_dir_rightwards(seg_i) = on_rightward_journey;
            all_seg_start(seg_i) = seg_start;
            all_seg_end(seg_i) = seg_end;
        elseif isempty(seg_end) && ~isempty(seg_start)
            if ~always_check_return
            seg_end = length(x);
            seg_i = seg_i+1;
            all_seg_dir_rightwards(seg_i) = on_rightward_journey;
            all_seg_start(seg_i) = seg_start;
            all_seg_end(seg_i) = seg_end;
            else % not complete don't count
%             seg_end = length(x);
%             seg_i = seg_i+1;
%             all_seg_dir_rightwards(seg_i) = on_rightward_journey;
%             all_seg_start(seg_i) = seg_start;
%             all_seg_end(seg_i) = seg_end;
            end
            break
        else
            break
        end
        seg_start = seg_end + find((short.inside(seg_end:end) | long.inside(seg_end:end)), 1)-1;
        last_z = max([find(right.inside(1:seg_start)==1, 1,'last') find(left.inside(1:seg_start)==1, 1,'last')]);
        if right.inside(last_z) == 1
            on_rightward_journey = false; % journey is from right to left
        elseif left.inside(last_z) == 1
            on_rightward_journey = true; % journey is from left to right
        end
        
    end
else
    fprintf('No place to start! Check zones and position data');
end

if turn_around_counter > 0
    turn_arounds.index = turn_peak(1:turn_around_counter);
    turn_arounds.type = turn_type(1:turn_around_counter);
    turn_arounds.key = turn_key;
    turn_arounds.start = turn_start;
    turn_arounds.end = turn_end;
else
    turn_arounds.index = [];
    turn_arounds.type = [];
    turn_arounds.key = turn_key;
    turn_arounds.start = [];
    turn_arounds.end = [];
end

sr_ind = 0;
sl_ind = 0;
lr_ind = 0;
ll_ind = 0;
%%
spd_win = 4;
xmid = x<=120 & x>=-120;
ymid = y<=120 & y>=20;
for i = 1:seg_i
    seg_start = all_seg_start(i);
    seg_end = all_seg_end(i);
    seg = seg_start:seg_end;
    
    seg_mid = intersect(seg_start:seg_end, find(xmid));
    spd = speed(seg_mid);% cm spd along the central x portion
%     spd = speed(seg_start:seg_end);% cm spd 
    if isempty(spd)
        spd = NaN;
    end
    switch all_seg_dir_rightwards(i)
        case false
                if ~isempty(long.inside(seg_start:seg_end))
                    long_prop = sum(long.inside(seg_start:seg_end))/length(long.inside(seg_start:seg_end));
                else
                    long_prop = 0;
                end
                if ~isempty(short.inside(seg_start:seg_end))
                    short_prop = sum(short.inside(seg_start:seg_end))/length(short.inside(seg_start:seg_end));
                else
                    short_prop = 0;
                end
                if short_prop > long_prop                % short path
                    sr_ind = sr_ind + 1;
                    SR.index(seg) = true;
                    SR.pass_time(sr_ind) = (seg_end - seg_start)*dt; % total time for pass
                    SR.pass_ind(sr_ind,:) = [seg_start, seg_end];
                    SR.pass_length(sr_ind,:) = seg_end - seg_start;
                    SR.min_spd(sr_ind,:) = min(spd);
                    SR.mean_spd(sr_ind,:) = mean(spd);
                    run_id = 1;
                else % if short_prop < long_prop
                    % long path
                    lr_ind = lr_ind + 1;
                    LR.index(seg) = true;
                    LR.pass_time(lr_ind) = (seg_end - seg_start)*dt; % total time for pass
                    LR.pass_ind(lr_ind,:) = [seg_start, seg_end];
                    LR.pass_length(lr_ind,:) = seg_end - seg_start;
                    
                    long_mid = intersect(seg_start:seg_end, find(xmid));
                    long_spd = speed(long_mid);% cm spd along the central x portion
                    left_mid = intersect(seg_start:seg_end, find(ymid & x<0));
                    left_spd = speed(left_mid);% cm spd along the central x portion
                    right_mid = intersect(seg_start:seg_end, find(ymid & x>0));
                    right_spd = speed(right_mid);% cm spd along the central x portion
                    spd = [left_spd; long_spd; right_spd];
                    if isempty(spd)
                        spd = NaN;
                    end
                    
                    LR.min_spd(lr_ind,:) = min(spd);
                    LR.mean_spd(lr_ind,:) = mean(spd);
                    run_id = 2;
                end
        case true
                if ~isempty(long.inside(seg_start:seg_end))
                    long_prop = sum(long.inside(seg_start:seg_end))/length(long.inside(seg_start:seg_end));
                else
                    long_prop = 0;
                end
                if ~isempty(short.inside(seg_start:seg_end))
                    short_prop = sum(short.inside(seg_start:seg_end))/length(short.inside(seg_start:seg_end));
                else
                    short_prop = 0;
                end
                if short_prop > long_prop
                    % short path
                    sl_ind = sl_ind + 1;
                    SL.index(seg) = true;
                    SL.pass_time(sl_ind) = (seg_end - seg_start)*dt; % total time for pass
                    SL.pass_ind(sl_ind,:) = [seg_start, seg_end];
                    SL.pass_length(sl_ind,:) = seg_end - seg_start;
                    SL.min_spd(sl_ind,:) = min(spd);
                    SL.mean_spd(sl_ind,:) = mean(spd);
                    run_id = 3;
                else % if short_prop < long_prop
                    % long path
                    ll_ind = ll_ind + 1;
                    LL.index(seg) = true;
                    LL.pass_time(ll_ind) = (seg_end - seg_start)*dt; % total time for pass
                    LL.pass_ind(ll_ind,:) = [seg_start, seg_end];
                    LL.pass_length(ll_ind,:) = seg_end - seg_start;
                    
                    long_mid = intersect(seg_start:seg_end, find(xmid));
                    long_spd = speed(long_mid);% cm spd along the central x portion
                    left_mid = intersect(seg_start:seg_end, find(ymid & x<0));
                    left_spd = speed(left_mid);% cm spd along the central x portion
                    right_mid = intersect(seg_start:seg_end, find(ymid & x>0));
                    right_spd = speed(right_mid);% cm spd along the central x portion
                    spd = [left_spd; long_spd; right_spd];
                    if isempty(spd)
                        spd = NaN;
                    end
                    LL.min_spd(ll_ind,:) = min(spd);
                    LL.mean_spd(ll_ind,:) = mean(spd);
                    run_id = 4;
            end
    end

    SL.master_ind(seg) = run_id;
    LL.master_ind(seg) = run_id;
    SR.master_ind(seg) = run_id;
    LR.master_ind(seg) = run_id;
end

% plot3(x,y,t)
if plotting
    subplot(122)
    hold on
    t = 1:length(x);
    plot3(x,y,t)
    i = SR.index;
    plot3(x(i),y(i),t(i), 'b.')
    i = SL.index;
    plot3(x(i),y(i),t(i), 'r.')
    i = LL.index;
    plot3(x(i),y(i),t(i), 'm.')
    i = LR.index;
    plot3(x(i),y(i),t(i), 'c.')
    set(gca, 'View', [-18 51.5])
end
%

% SR.pass_ind = SR.pass_ind(~isnan(SR.pass_ind(:,1)),:);
% SL.pass_ind = SL.pass_ind(~isnan(SL.pass_ind(:,1)),:);
% LR.pass_ind = LR.pass_ind(~isnan(LR.pass_ind(:,1)),:);
% LL.pass_ind = LL.pass_ind(~isnan(LL.pass_ind(:,1)),:);
SL.total_time = sum(SL.index)*dt;
SR.total_time = sum(SR.index)*dt;
LL.total_time = sum(LL.index)*dt;
LR.total_time = sum(LR.index)*dt;
% % % % % % % % % % % % % % Linear_plots