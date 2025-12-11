function [ZoneStruct] = DUAL8MAZE_performance_eval(room, arena, params, debug)
debug = false;
angle_cut = 15; % cutoff angle to define choice and arm sections
rho_cut   = 30;%25; % cutoff distance between choice and center


%% define trajectories along the arena
% x=arena.rawx; y=arena.rawy; ts=arena.timestamps/1000
if params.during_rotation_only == true;
    arena_angdiff = movmean(unwrap(arena.angular_offset-arena.angular_offset(1)), 30*.5);
    arena_rotate_start = find( abs(arena_angdiff)> pi/4 , 1);
    arena.x(1:arena_rotate_start) = NaN;
    arena.y(1:arena_rotate_start) = NaN;
end


[room.pol_theta, room.pol_rho] = cart2pol(room.x, room.y);
[arena.pol_theta, arena.pol_rho] = cart2pol(arena.x, arena.y);

% x=arena.x; y=arena.y; ts=arena.timestamps/1000;
% x=room.x; y=room.y; ts=room.timestamps/1000;
t = arena.pol_theta -pi/2;
t(t<-pi)=t(t<-pi)+2*pi;
t = (t)*180/pi;
rho = arena.pol_rho;

center_zone = rho <= rho_cut;
top_choice_zone = (rho >rho_cut) & (abs(t) <= angle_cut);
bot_choice_zone = (rho >rho_cut) & (abs(t) >= 180-angle_cut);
either_choice_zone = top_choice_zone | bot_choice_zone;
right_zone  = (rho > rho_cut) & (t < -angle_cut) & (t > -180+angle_cut);
left_zone = (rho > rho_cut) & (t >  angle_cut) & (t <  180-angle_cut);
either_arm_zone = left_zone | right_zone;
strnames = {'room', 'arena'};

% room.entrance_start_idx  = -1; % initialized in case no reinforcement
% arena.entrance_start_idx = -1; % initialized in case no reinforcement

for strLoop = 1:length(strnames)
    %%
    % input structure
    eval(sprintf('behav = %s;', strnames{strLoop}));
    if ~isfield(behav, 'DAT_fileinfo')
        behav.DAT_fileinfo = false;
        behav.entrance_start_idx = -1;
        eval(sprintf('%s = behav;', strnames{strLoop}));
    end
    if isfield(behav.DAT_fileinfo, 'ReinforcedSector')
        %%
        % find entries into the reinforcement zone
%         z1 = mod(behav.DAT_fileinfo.ReinforcedSector(1)-behav.DAT_fileinfo.ReinforcedSector(2)/2, 360);
%         z2 = mod(behav.DAT_fileinfo.ReinforcedSector(1)+behav.DAT_fileinfo.ReinforcedSector(2)/2, 360);
%         d = mod(rad2deg(behav.pol_theta+pi)+behav.DAT_fileinfo.ReinforcedSector(1), 360);
%         d = mod(rad2deg(behav.pol_theta+pi)+180, 360);
        
        d = mod(rad2deg(behav.pol_theta) - behav.DAT_fileinfo.ReinforcedSector(1)+180, 360);
        [~, in_zone_theta] = angular_distance(d, 180, behav.DAT_fileinfo.ReinforcedSector(2)/2);
%         in_zone_theta =  ( d >= z1 )& ( d <= z2 );
        r1 = rho_cut*behav.DAT_fileinfo.TrackerResolution_PixPerCM; % behav.DAT_fileinfo.ReinforcedSector(3);
        % r2 = room.DAT_fileinfo.RoomTrackReinforcedSector(4);
        r = behav.pol_rho*behav.DAT_fileinfo.TrackerResolution_PixPerCM;
        in_zone_rho =  ( r >= r1 );% & ( r <=  r2);
        
        in_real_zone = in_zone_theta & in_zone_rho;
        
        % find entries into the a fake zone on the opposite area
%         z1 = mod(behav.DAT_fileinfo.ReinforcedSector(1)-behav.DAT_fileinfo.ReinforcedSector(2)/2, 360);
%         z2 = mod(behav.DAT_fileinfo.ReinforcedSector(1)+behav.DAT_fileinfo.ReinforcedSector(2)/2, 360);
%         faked = mod(rad2deg(behav.pol_theta+pi+pi)+180, 360);
%         in_fake_zone =  ( faked >= z1 )& ( faked <= z2 );
        
%         faked = mod(rad2deg(behav.pol_theta+pi+pi)+180, 360);
        faked = mod(rad2deg(behav.pol_theta) - behav.DAT_fileinfo.ReinforcedSector(1), 360);
        [~, in_zone_opposite] = angular_distance(faked, 180, behav.DAT_fileinfo.ReinforcedSector(2)/2);;
        in_fake_zone = in_zone_opposite & in_zone_rho;
        
        
        % second_conv = round(.5/median(room.dt));
        % in_real_zone = conv(in_real_zone, ones(second_conv,1)./second_conv, 'same')>0;
        % in_fake_zone = conv(in_fake_zone, ones(second_conv,1)./second_conv, 'same')>0;
        
        enter_real = in_real_zone<0;
        enter_real(2:end) = in_real_zone(2:end)==1 & in_real_zone(1:end-1)==0;
        enter_fake = in_fake_zone<0;
        enter_fake(2:end) = in_fake_zone(2:end)==1 & in_fake_zone(1:end-1)==0;
        
        %%%%%%
        pref = sum(enter_real) / sum(enter_real|enter_fake);
        entrance_real_idx = find(enter_real);
        entrance_fake_idx = find(enter_fake);
        
        time_to_retrieve = NaN(length(entrance_real_idx),1);
        retrieve_idx = NaN(length(entrance_real_idx),1);
        time_to_retrieve_fake = NaN(length(entrance_fake_idx),1);
        fake_retrieve_idx = NaN(length(entrance_fake_idx),1);
        ts = behav.timestamps./1000;
        for ind = 1:length(entrance_real_idx)%length(room.entrance_start_idx)
            e = entrance_real_idx(ind);
            r = find(center_zone(e:end), 1, 'first')+e;
            if ~isempty(r)
                time_to_retrieve(ind) =  ts(r) - ts(e);
                retrieve_idx(ind) = r;
            end
        end
        for ind = 1:length(entrance_fake_idx)%length(room.entrance_start_idx)
            e = entrance_fake_idx(ind);
            r = find(center_zone(e:end), 1, 'first')+e;
            if ~isempty(r)
                time_to_retrieve_fake(ind) = ts(r) - ts(e);
                fake_retrieve_idx(ind) = r;
            end
        end
        bads = isnan(entrance_fake_idx) | isnan(fake_retrieve_idx);
        entrance_fake_idx = entrance_fake_idx(~bads);
        fake_retrieve_idx = fake_retrieve_idx(~bads);
        bads = isnan(entrance_real_idx) | isnan(retrieve_idx);
        entrance_real_idx = entrance_real_idx(~bads);
        retrieve_idx      = retrieve_idx(~bads);
        
%         figure; 
%         plot(behav.x, behav.y)
%         hold on
%         plot(behav.x(in_real_zone), behav.y(in_real_zone), 'ro')
%         plot(behav.x(in_fake_zone), behav.y(in_fake_zone), 'ks')
%         plot(behav.x(behav.entrance_start_idx), behav.y(behav.entrance_start_idx), 'go')
        
        behav.reinforced_pref          = pref;
        behav.entrance_real_idx         = entrance_real_idx;
        behav.entrance_fake_idx         = entrance_fake_idx;
        behav.time_to_retrieve          = time_to_retrieve;
        behav.time_to_retrieve_fake     = time_to_retrieve_fake;
        
        if debug %%%%%%%%%%%% define entries into the reinforced zone and a fake control zone
            %     close all
            figure;
            subplot(2,1,1)
            hold on; plot(d)
            plot(find(in_real_zone), d(in_real_zone), 'g.');
            plot(find(enter_real), d(enter_real), 'gx');
            plot(find(in_fake_zone), d(in_fake_zone), 'm.');
            plot(find(enter_fake), d(enter_fake), 'mx');
            subplot(2,1,2)
            hold on; plot(room.x, room.y)
            plot(room.x(enter_real), room.y(enter_real), 'gx');
            plot(room.x(enter_fake), room.y(enter_fake), 'mo');
            
            figure;
            subplot(121); hold on;
            x = arena.x; y = arena.y;
            plot3(x,y, ts/1000, 'k-')
            plot3(x(center_zone), y(center_zone), ts(center_zone)/1000, 'g.','MarkerSize', 5)
            plot3(x(either_choice_zone), y(either_choice_zone), ts(either_choice_zone)/1000, 'm.','MarkerSize', 5)
            plot3(x(left_zone), y(left_zone), ts(left_zone)/1000, 'r.','MarkerSize', 5)
            plot3(x(right_zone), y(right_zone), ts(right_zone)/1000, 'b.','MarkerSize', 5)
            drawnow
            subplot(122); hold on
            x = room.x; y = room.y;
            plot3(x,y, ts/1000, 'k-')
            plot3(x(center_zone), y(center_zone), ts(center_zone)/1000, 'g.','MarkerSize', 5)
            plot3(x(either_choice_zone), y(either_choice_zone), ts(either_choice_zone)/1000, 'm.','MarkerSize', 5)
            plot3(x(left_zone), y(left_zone), ts(left_zone)/1000, 'r.','MarkerSize', 5)
            plot3(x(right_zone), y(right_zone), ts(right_zone)/1000, 'b.','MarkerSize', 5)
            drawnow
            
            figure
            subplot(121); hold on
            plot3(room.x,room.y, room.timestamps/1000, 'k-')
            scatter3(room.x(room.entrance_start_idx),room.y(room.entrance_start_idx), room.timestamps(room.entrance_start_idx)/1000, 500, 'r.')
            ylabel('Y-pos'); xlabel('X-pos'); zlabel('Session time (min)')
            set(gca, 'View', [30.3407   62.8791], 'ZTick', [0:300:1800], 'ZTickLabel', [0:5:30], 'FontWeight', 'bold', 'FontName', 'Arial bold')
            subplot(122); hold on
            plot3(arena.x,arena.y, arena.timestamps/1000, 'k-')
            scatter3(arena.x(room.entrance_start_idx), arena.y(room.entrance_start_idx), arena.timestamps(room.entrance_start_idx)/1000, 500, 'r.')
            ylabel('Y-pos'); xlabel('X-pos'); zlabel('Session time (min)')
            set(gca, 'View', [30.3407   62.8791], 'ZTick', [0:300:1800], 'ZTickLabel', [0:5:30], 'FontWeight', 'bold', 'FontName', 'Arial bold')
        end
        % output structure
        eval(sprintf('%s = behav;', strnames{strLoop}));
    end
end
%%    
    
zone_sequence = arena.x*0;
zone_change = arena.x*0;
zone_change_end = arena.x*0;
% % version for partial left and rights
% zone_sequence(center_zone) = 1;
% zone_sequence(either_choice_zone) = 2;
% zone_sequence(right_zone) = 3;
% zone_sequence(left_zone) = 4;
% % version for partial left and rights
zone_sequence(center_zone)      =  1;
zone_sequence(top_choice_zone)  =  2;
zone_sequence(bot_choice_zone)	= -2;
zone_sequence(left_zone)        =  3;
zone_sequence(right_zone)       = -3;

zone_change(2:end)          = zone_sequence(2:end)     ~= zone_sequence(1:end-1);
zone_change_end(1:end-1)    = zone_sequence(1:end-1)   ~= zone_sequence(2:end);

zs = zone_sequence(zone_change==true);
zs_ind = find(zone_change==true);
zs_ind2 = find(zone_change_end==true);

%%%%% finding partial trajectories through arms, center to choice
L_seq = [1,2, 3,-2; 1,-2, 3,2]';
R_seq = [1,2,-3,-2; 1,-2,-3,2]';
[L_start, L_end, idxL] = seq_finder(L_seq, zs, zs_ind, zs_ind2);
[R_start, R_end, idxR] = seq_finder(R_seq, zs, zs_ind, zs_ind2);
% %%%%% finding full trajectories through arms center to center
% L_seq = [1,2,3,-2,1;  1,-2,3,2,1 ]';
% R_seq = [1,2,-3,-2,1; 1,-2,-3,2,1]';
% [L_start, L_end, idxL] = seq_finder(L_seq, zs, zs_ind, zs_ind2);
% [R_start, R_end, idxR] = seq_finder(R_seq, zs, zs_ind, zs_ind2);
% % setdiff(L_start_partial, L_start)
% % setdiff(full_moon_start, [L_start_partial; R_start_partial] )

[Correct_start, Correct_end, Error_start, Error_end, frame_type, is_path_rew] = ...
correct_choice_eval(L_start, L_end, R_start, R_end, room.entrance_start_idx, arena.entrance_start_idx);

L_full_seq = [1,  2,  3, -2, -3,  2, 1;
                   1, -2,  3,  2, -3, -2, 1]';
R_full_seq = [1,  2, -3, -2,  3,  2, 1;
                   1, -2, -3,  2,  3, -2, 1]';
% [full_moon_start, full_moon_end, idx_fullmoon] = seq_finder(full_moon_seq, zs, zs_ind, zs_ind2);
[L_full_start, L_full_end, ~] = seq_finder(L_full_seq, zs, zs_ind, zs_ind2);
[R_full_start, R_full_end, ~] = seq_finder(R_full_seq, zs, zs_ind, zs_ind2);

% cutoff full circles at second choice entry
for ii = 1:length(L_full_start)
%         idx = R_full_start(ii):R_full_end(ii);
        idx2 = find(zs_ind == L_full_start(ii)):find(zs_ind > L_full_end(ii), 1, 'first');
        zs_sub = abs( zone_sequence(zs_ind(idx2(1)):zs_ind(idx2(end))) );
        three_two = find(zs_sub(1:end-1)==3 & zs_sub(2:end)==2, 1, 'first');
        two_three = find(zs_sub(three_two:end-1)==2 & zs_sub(three_two+1:end)==3) + three_two;
        L_full_end(ii) = round( (two_three+ three_two)/2) + zs_ind(idx2(1));
end
for ii = 1:length(R_full_start)
%         idx = R_full_start(ii):R_full_end(ii);
        idx2 = find(zs_ind == R_full_start(ii)):find(zs_ind > R_full_end(ii), 1, 'first');
        zs_sub = abs( zone_sequence(zs_ind(idx2(1)):zs_ind(idx2(end))) );
        three_two = find(zs_sub(1:end-1)==3 & zs_sub(2:end)==2, 1, 'first');
        two_three = find(zs_sub(three_two:end-1)==2 & zs_sub(three_two+1:end)==3) + three_two;
        R_full_end(ii) = round( (two_three+ three_two)/2) + zs_ind(idx2(1));
end
[Correct_full_start, Correct_full_end, Error_full_start, Error_full_end, frame_full_type, is_path_rew_full] = ...
correct_choice_eval(L_full_start, L_full_end, R_full_start, R_full_end, room.entrance_start_idx, arena.entrance_start_idx);

if params.strict_errors_only == false
    if false
        isFull        = [Correct_full_start; Error_full_start];
        isHalf        = [Correct_start; Error_start];
        Correct_start = [Correct_start; Correct_full_start];
        Correct_end   = [Correct_end;   Correct_full_end];
        frame_type    = [frame_type;   frame_full_type];
    else % only include full circle errors?
        isFull      = [Error_full_start];
        isHalf      = [Error_start];
    end
    Error_start = [Error_start; Error_full_start];
    Error_end   = [Error_end;   Error_full_end];
    L_start     = [L_start; intersect(L_full_start, Error_full_start)];
    L_end       = [L_end;   intersect(L_full_end,   Error_full_end)];
    R_start     = [R_start; intersect(R_full_start, Error_full_start)];
    R_end       = [R_end;   intersect(R_full_end,   Error_full_end)];
else
    if false
        isFull      = [];
        isHalf      = [Correct_start; Error_start];
    else % only include full circle errors?
        isFull      = [];
        isHalf      = [Error_start];
    end
end


%%

% rews = room.entrance_start_idx;
if isfield(room.DAT_fileinfo, 'ReinforcedSector') && ~isfield(arena.DAT_fileinfo, 'ReinforcedSector')
    rr = deg2rad(room.DAT_fileinfo.ReinforcedSector(1))*ones(length(arena.angular_offset),1);
elseif ~isfield(room.DAT_fileinfo, 'ReinforcedSector') && isfield(arena.DAT_fileinfo, 'ReinforcedSector')
    rr = deg2rad(arena.DAT_fileinfo.ReinforcedSector(1)) + arena.angular_offset;
else % no reinforcement, habituation
    rr = zeros(length(arena.x),1);
%     rr = zeros(length(arena.angular_offset),1);
end
debug1=false;
if debug1; figure; end
x = room.x;  y = room.y;  ts = room.timestamps./1000;
[rew_diff_R, rew_diff_idx_R] = Tjunction_difficulty(room.x, room.y, ...
    room.timestamps./1000, rr, zone_sequence, R_start, R_end, debug1);
[rew_diff_L, rew_diff_idx_L] = Tjunction_difficulty(room.x, room.y, ...
    room.timestamps./1000, rr, zone_sequence, L_start, L_end, debug1);

rew_diff = [rew_diff_R; rew_diff_L];
rew_diff_idx = [rew_diff_idx_R; rew_diff_idx_L];

if debug1
scatter3(x(rew_diff_idx_R), y(rew_diff_idx_R), ts(rew_diff_idx_R), 20+(abs(rad2deg(rew_diff_R)))*10, 'm.')
scatter3(x(rew_diff_idx_L), y(rew_diff_idx_L), ts(rew_diff_idx_L), 20+(abs(rad2deg(rew_diff_L)))*10, 'c.')
x = arena.x;  y = arena.y;  ts = arena.timestamps./1000;
subplot(122); hold on; 
plot3(x, y, ts, 'k-')
scatter3(x(rew_diff_idx_R), y(rew_diff_idx_R), ts(rew_diff_idx_R), 20+(abs(rad2deg(rew_diff_R)))*10, 'm.')
scatter3(x(rew_diff_idx_L), y(rew_diff_idx_L), ts(rew_diff_idx_L), 20+(abs(rad2deg(rew_diff_L)))*10, 'c.')
end
% debug1=debug;
%%
isCorrect                   = zeros(length(room.x), 1);
isAngDist                   = zeros(length(room.x), 1);
isChoiceIDX                 = zeros(length(room.x), 1);
isLeft                      = zeros(length(room.x), 1);
isFullCircle                = zeros(length(room.x), 1);
isChoiceIDX(rew_diff_idx)   =  rew_diff_idx;
isAngDist(rew_diff_idx)     =  rew_diff;
isCorrect(Correct_start)    =  1;
isCorrect(Error_start)      = -1;
isLeft(L_start)             =  1;
isLeft(R_start)             = -1;
isFullCircle(isFull)        =  1;
isFullCircle(isHalf)        = -1;
isIDX                       = find(isCorrect~=0);
isCorrect                   = isCorrect(isIDX)>0;
isLeft                      = isLeft(isIDX)>0;
isAngDist                   = isAngDist(isChoiceIDX>0);
isChoiceIDX                 = isChoiceIDX(isChoiceIDX>0);
isFullCircle                = isFullCircle(isIDX)>0;
isTime                      = arena.timestamps(isIDX)./1000;

if debug1
    figure(); hold on
    x = arena.x;  y = arena.y;  ts = arena.timestamps./1000;
    plot3(x, y, ts, 'k-')
    idx = isChoiceIDX(isCorrect);
    sz  = isAngDist(isCorrect);
    scatter3(x(idx), y(idx), ts(idx), 20+(abs(rad2deg(sz)))*10, 'g.')
    idx = room.entrance_start_idx(is_path_rew);
    scatter3(x(idx), y(idx), ts(idx), 20+(abs(rad2deg(sz)))*10, 'c.')
    
    idx = isChoiceIDX(~isCorrect);
    sz  = isAngDist(~isCorrect);
    scatter3(x(idx), y(idx), ts(idx), 20+(abs(rad2deg(sz)))*10, 'r.')
end

streaks     = isCorrect*0;
streaks_err = isCorrect*0;
for ii = 1:length(isCorrect)
    idx = find(isCorrect(1:ii)==0, 1, 'last');
    if ~isempty(idx)
        streaks(ii) = sum(isCorrect(idx:ii)== true);
    end
    idx = find(isCorrect(1:ii)==1, 1, 'last');
    if ~isempty(idx)
        streaks_err(ii) = sum(isCorrect(idx:ii)== false);
    end
end

choices = [];
choices.L_start_end         = [L_start, L_end];
choices.R_start_end         = [R_start, R_end];
choices.Correct_start_end   = [Correct_start, Correct_end];
choices.Error_start_end     = [Error_start, Error_end];
choices.reinforced_frame    = frame_type;
choices.isLeft              = isLeft;
choices.choiceDist          = isAngDist;
choices.choiceIDX           = isChoiceIDX;
choices.isCorrect           = isCorrect;
choices.isIDX               = isIDX ;
choices.isFullCircle        = isFullCircle ;
choices.isTime              = isTime;
choices.streaks             = streaks;
choices.streaks_err         = streaks_err;

dxy = sqrt((diff(x).^2) + (diff(y).^2));
dxy = [dxy(1); dxy];
linear_vec = zeros(length(arena.x),1);
left_vec = false(length(arena.x),1);
right_vec = false(length(arena.x),1);
if choices.L_start_end>0
for jj = 1:size(choices.L_start_end, 1)
    inds = choices.L_start_end(jj,1):choices.L_start_end(jj,2);
    left_vec(inds) = true;
    poschange = cumsum(dxy(inds)) ./ sum(dxy(inds));
    linear_vec(inds) = poschange;
end
end
if choices.R_start_end>0
for jj = 1:size(choices.R_start_end, 1)
    inds = choices.R_start_end(jj,1):choices.R_start_end(jj,2);
    right_vec(inds) = true;
    poschange = cumsum(dxy(inds)) ./ sum(dxy(inds));
%     plot(linspace(0, 1, length(inds)), poschange)
    linear_vec(inds) = poschange;
end
end
correct_vec = false(length(arena.x),1);
error_vec = false(length(arena.x),1);
if choices.Correct_start_end>0
for jj = 1:size(choices.Correct_start_end, 1)
    inds = choices.Correct_start_end(jj,1):choices.Correct_start_end(jj,2);
    correct_vec(inds) = true;
end
end
if choices.Error_start_end>0
for jj = 1:size(choices.Error_start_end, 1)
    inds = choices.Error_start_end(jj,1):choices.Error_start_end(jj,2);
    error_vec(inds) = true;
end
end
choices.paths.left    = left_vec;
choices.paths.right   = right_vec;
choices.paths.correct = correct_vec;
choices.paths.error   = error_vec;
choices.paths.linear_vec   = linear_vec;



if debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    figure; 
    subplot(121)
    hold on;
% %     correct_runav = cumsum(isCorrect)/length(abs(isCorrect));
% %     left_runav = cumsum(isLeft)/length(abs(isLeft));
%     correct_runav = cumsum(isCorrect==1)/length(abs(isCorrect));
%     left_runav = cumsum(isLeft==1)/length(abs(isLeft));
    
    correct_runav = 0*isCorrect;
    for ii = 1:length(isCorrect)
        correct_runav(ii) = sum(isCorrect(1:ii)==1)/ii;
    end
    left_runav = 0*isLeft;
    for ii = 1:length(isLeft)
        left_runav(ii) = sum(isLeft(1:ii)==1)/ii;
    end
    
    plot(isTime/60, correct_runav, 'g-', 'Color', [.4 1 .4]/2);
    plot(isTime/60, left_runav, 'b--')
    plot([0 ts(end)]/60, [.5 .5], 'k:')
%     plot([0 ts(end)]/60, [0 1], 'k:')
%     plot([0 ts(end)]/60, [0 -1], 'k:')
    ylim([-.1 1.1])
    xlim([0 ts(end)]/60)
    ylabel('Choice performance')
    xlabel('Session time (minutes)')
    
    x=arena.x; y=arena.y; ts=arena.timestamps./1000;
    subplot(122)
    hold on;
    plot3(x,y,ts, 'k-')
    for ii = 1:length(L_start)
        idx = L_start(ii):L_end(ii);
        plot3(x(idx) ,y(idx),ts(idx), 'b-')
    end
    for ii = 1:length(R_start)
        idx = R_start(ii):R_end(ii);
        plot3(x(idx) ,y(idx),ts(idx), 'r-')
    end
%     for ii = 1:length(L_full_start)
%         idx = L_full_start(ii):L_full_end(ii);
%         plot3(x(idx) ,y(idx),ts(idx), 'b--', 'LineWidth', 2)
%     end
%     for ii = 1:length(R_full_start)
%         idx = R_full_start(ii):R_full_end(ii);
%         plot3(x(idx) ,y(idx),ts(idx), 'm-', 'LineWidth', 2)
%     end

%     for ii = 1:length(R_start_partial)
%         idx = R_start_partial(ii):R_end_partial(ii);
%         plot3(x(idx) ,y(idx),ts(idx), 'r-')
%     end
    for ii = 1:length(Correct_start)
        idx = Correct_start(ii):Correct_end(ii);
        plot3(x(idx) ,y(idx),ts(idx), 'g.')
    end
end
%%
lose1_win2  = isCorrect(1:end-1)==0 & isCorrect(2:end)==1;
lose1_lose2 = isCorrect(1:end-1)==0 & isCorrect(2:end)==0;
win1_win2   = isCorrect(1:end-1)==1 & isCorrect(2:end)==1;
win1_lose2  = isCorrect(1:end-1)==1 & isCorrect(2:end)==0;
choice_perf = zeros(2,2);
choice_perf(1,1) = sum(win1_win2)/length(isCorrect);
choice_perf(1,2) = sum(win1_lose2)/length(isCorrect);
choice_perf(2,2) = sum(lose1_lose2)/length(isCorrect);
choice_perf(2,1) = sum(lose1_win2)/length(isCorrect);

lr_swap = isLeft(1:end-1)==1 & isLeft(2:end)==0;
rl_swap = isLeft(2:end)==1   & isLeft(1:end-1)==0;
dir_swap = [0; lr_swap | rl_swap];
dir_swap_idx = find(dir_swap==1);

lose1_win2  = isCorrect(dir_swap_idx-1)==0 & isCorrect(dir_swap_idx)==1;
lose1_lose2 = isCorrect(dir_swap_idx-1)==0 & isCorrect(dir_swap_idx)==0;
win1_win2   = isCorrect(dir_swap_idx-1)==1 & isCorrect(dir_swap_idx)==1;
win1_lose2  = isCorrect(dir_swap_idx-1)==1 & isCorrect(dir_swap_idx)==0;
swap_perf = zeros(2,2);
swap_perf(1,1) = sum(win1_win2)/length(dir_swap_idx);
swap_perf(1,2) = sum(win1_lose2)/length(dir_swap_idx);
swap_perf(2,2) = sum(lose1_lose2)/length(dir_swap_idx);
swap_perf(2,1) = sum(lose1_win2)/length(dir_swap_idx);

% % bootstrapping by scrambling seems wrong direction
% swap_perf_rand = zeros(2,2, 100);
% for ii = 1:1000
% %     [~, randord] = sort(rand(length(isLeft),1));
% %     isLeft_rand    = isLeft(randord);
% %     [~, randord] = sort(rand(length(isLeft),1));
% %     isCorrect_rand = isCorrect(randord);
%     isLeft_rand    = circshift(isLeft, randi([-length(isLeft)+1 length(isLeft)-1]));
%     isCorrect_rand = circshift(isCorrect, randi([-length(isCorrect)+1 length(isCorrect)-1]));
%     lr_swap = isLeft_rand(1:end-1)==1 & isLeft_rand(2:end)==-1;
%     rl_swap = isLeft_rand(2:end)==1   & isLeft_rand(1:end-1)==-1;
%     dir_swap = [0; lr_swap | rl_swap];
%     dir_swap_idx_rand = find(dir_swap==1);
%
%     lose1_win2  = isCorrect_rand(dir_swap_idx_rand-1)==-1 & isCorrect_rand(dir_swap_idx_rand)==1;
%     lose1_lose2 = isCorrect_rand(dir_swap_idx_rand-1)==-1 & isCorrect_rand(dir_swap_idx_rand)==-1;
%     win1_win2   = isCorrect_rand(dir_swap_idx_rand-1)==1 & isCorrect_rand(dir_swap_idx_rand)==1;
%     win1_lose2  = isCorrect_rand(dir_swap_idx_rand-1)==1 & isCorrect_rand(dir_swap_idx_rand)==-1;
%     swap_perf_rand(1,1, ii) = sum(win1_win2)/length(dir_swap_idx_rand);
%     swap_perf_rand(1,2, ii) = sum(win1_lose2)/length(dir_swap_idx_rand);
%     swap_perf_rand(2,2, ii) = sum(lose1_lose2)/length(dir_swap_idx_rand);
%     swap_perf_rand(2,1, ii) = sum(lose1_win2)/length(dir_swap_idx_rand);
% end
%%%%%% assuming lose-shift just yields swap_perf_fake=[0 0; 1 0] of course
% fake_dir_swap = dir_swap*0;
% fake_dir_swap(1+find(isCorrect(1:end-1)==-1 & isCorrect(2:end)==1)) = 1; % assume lose win strat
% dir_swap_idx = find(fake_dir_swap==1);
%
% lose1_win2  = isCorrect(dir_swap_idx-1)==-1 & isCorrect(dir_swap_idx)==1;
% lose1_lose2 = isCorrect(dir_swap_idx-1)==-1 & isCorrect(dir_swap_idx)==-1;
% win1_win2   = isCorrect(dir_swap_idx-1)==1 & isCorrect(dir_swap_idx)==1;
% win1_lose2  = isCorrect(dir_swap_idx-1)==1 & isCorrect(dir_swap_idx)==-1;
% swap_perf_fake = zeros(2,2);
% swap_perf_fake(1,1) = sum(win1_win2)/length(dir_swap_idx);
% swap_perf_fake(1,2) = sum(win1_lose2)/length(dir_swap_idx);
% swap_perf_fake(2,2) = sum(lose1_lose2)/length(dir_swap_idx);
% (2,1) = sum(lose1_win2)/length(dir_swap_idx);

%%
Correct_start = sort(Correct_start, 'ascend');
Correct_end = sort(Correct_end, 'ascend');
Error_start = sort(Error_start, 'ascend');
Error_end = sort(Error_end, 'ascend');

dt = 5*60; % minutes
t2 = round(ts(end)/60)*60;
t_inds = 0:dt:t2;
prob_L = zeros(length(t_inds)-1, 1);
prob_C = zeros(length(t_inds)-1, 1);

n_thresh = 4; % if less than this in period, will use the last n
for ind= 2:length(t_inds)
    ls = ts(L_start);
    rs = ts(R_start);
    cs = ts(Correct_start);
    es = ts(Error_start);
    
    sub_left = sum(ls <= t_inds(ind) & ls > t_inds(ind-1));
    sub_right = sum(rs <= t_inds(ind) & rs > t_inds(ind-1));
    sub_correct = sum(cs <= t_inds(ind) & cs > t_inds(ind-1));
    sub_incorrect = sum(es <= t_inds(ind) & es > t_inds(ind-1));
    n_trials = (sub_left + sub_right);
    if n_trials >= n_thresh
        %         sub_left      = ls <= t_inds(ind);
        %         sub_left      = sum(sub_left(1:n_trials));
        %         sub_right     = rs <= t_inds(ind);
        %         sub_right     = sum(sub_right(1:n_trials));
        %
        %         sub_correct   = cs <= t_inds(ind);
        %         sub_correct   = sum(sub_correct(1:n_trials));
        %         sub_incorrect = es <= t_inds(ind);
        %         sub_incorrect = sum(sub_incorrect(end-n_trials:end));
        prob_C(ind-1) = sub_correct / (sub_correct + sub_incorrect);
        prob_L(ind-1) = sub_left / (sub_left + sub_right);
    end
end
if debug
    figure; hold on;
    plot(t_inds(2:end)/60, prob_L, 'b.-')
    plot(t_inds(2:end)/60, prob_C, 'go--')
    plot(t_inds(2:end)/60, prob_C*0 + .5, 'k:')
    ylim([ -.1 1.1])
    xlim([ t_inds(1), t_inds(end)+dt]/60)
end
% time_to_retrieve = NaN(length(room.entrance_start_idx),1);

ZoneStruct = [];
ZoneStruct.arena_zones.angle_cut = angle_cut;
ZoneStruct.arena_zones.rho_cut   = rho_cut;
ZoneStruct.arena_zones.center = center_zone;
ZoneStruct.arena_zones.top    = top_choice_zone;
ZoneStruct.arena_zones.bot    = bot_choice_zone;
ZoneStruct.arena_zones.left   = left_zone;
ZoneStruct.arena_zones.right  = right_zone;
ZoneStruct.arena_zones_seq    = zone_sequence;
ZoneStruct.arena_zones_labels = {'center_1', 'either_choice_2', 'right_3', 'left_4'};
ZoneStruct.arena_zones.left_start  = L_start;
ZoneStruct.arena_zones.left_end    = L_end;
ZoneStruct.arena_zones.right_start = R_start;
ZoneStruct.arena_zones.right_end   = R_end;
ZoneStruct.metrics.Leftarm_bias    = idxL/(idxL+idxR);
ZoneStruct.metrics.sequence_prob_Left   = prob_L;
ZoneStruct.metrics.sequence_prob_Correct   = prob_C;
ZoneStruct.metrics.sequence_prob_times   = t_inds(2:end);
ZoneStruct.metrics.choice_performance      = choice_perf;
ZoneStruct.metrics.swap_performance        = swap_perf;
ZoneStruct.metrics.choice_swap_key    = {'win-win', 'win-lose'; 'lose-win', 'lose-lose'};
ZoneStruct.metrics.average_streak     = mean(streaks(streaks>0));
ZoneStruct.metrics.average_streak_err = mean(streaks_err(streaks_err>0));
ZoneStruct.choices = choices;

if isfield(room.DAT_fileinfo, 'ReinforcedSector') && ~isfield(arena.DAT_fileinfo, 'ReinforcedSector')
    ZoneStruct.metrics.reinforced_zone_preference   = room.reinforced_pref;
    ZoneStruct.metrics.mean_retrieval_time          = nanmean(room.time_to_retrieve);
    ZoneStruct.metrics.mean_retrieval_time_fake     = nanmean(room.time_to_retrieve_fake);
elseif ~isfield(room.DAT_fileinfo, 'ReinforcedSector') && isfield(arena.DAT_fileinfo, 'ReinforcedSector')
    ZoneStruct.metrics.reinforced_zone_preference   = arena.reinforced_pref;
    ZoneStruct.metrics.mean_retrieval_time          = nanmean(arena.time_to_retrieve);
    ZoneStruct.metrics.mean_retrieval_time_fake     = nanmean(arena.time_to_retrieve_fake);
else % no reinforcement, habituation
    ZoneStruct.metrics.reinforced_zone_preference   = NaN;
    ZoneStruct.metrics.mean_retrieval_time          = NaN;
    ZoneStruct.metrics.mean_retrieval_time_fake     = NaN;
end

%% %%%%%%%%%%%%%%
if true
ang_bins = linspace(-pi, pi, 13);
rho_thresh = 0;
% inds = abs(zone_sequence)==3;
inds = any(abs(zone_sequence)==[2,3], 2);
x = room.x; y = room.y;
x(inds==false) = NaN;
y(inds==false) = NaN;
[r_pxy, r_pt, rt, room_resultant] = path_visit_occupancy(x,  y,  params.pos_bins, ang_bins, rho_thresh, L_start, L_end, R_start, R_end);
% r_pt = histcounts(rt, ang_bins, 'Normalization', 'probability');

x = arena.x; y = arena.y;
x(inds==false) = NaN;
y(inds==false) = NaN;
[a_pxy, a_pt, at, arena_resultant] = path_visit_occupancy(x, y, params.pos_bins, ang_bins, rho_thresh, L_start, L_end, R_start, R_end);
% a_pt = histcounts(at, ang_bins, 'Normalization', 'probability');
ang_bins_center  = ang_bins(1:end-1) + median(abs(diff(ang_bins)));
xs = rad2deg(ang_bins_center);

selectivity_room = mean(r_pt(r_pt>=quantile(r_pt, .75))) - mean(r_pt(r_pt<=quantile(r_pt, .25)));
selectivity_arena = mean(r_pt(r_pt>=quantile(r_pt, .75))) - mean(r_pt(r_pt<=quantile(r_pt, .25)));

ZoneStruct.metrics.sampling_room = r_pt;
ZoneStruct.metrics.sampling_room_xy = r_pxy;
ZoneStruct.metrics.room_resultant_prob = [resultantAngle(xs(r_pt>=quantile(r_pt, .75))), selectivity_room]; % room_resultant
ZoneStruct.metrics.arena_resultant_prob =[resultantAngle(xs(a_pt>=quantile(a_pt, .75))), selectivity_arena]; % arena_resultant;
ZoneStruct.metrics.sampling_arena = a_pt;
ZoneStruct.metrics.sampling_arena_xy = a_pxy;
ZoneStruct.metrics.sampling_room_arena_bins = ang_bins;

mindist = 0;
valid = ZoneStruct.choices.choiceDist>=mindist;%
sessmin = room.timestamps(end)/60000;
ZoneStruct.metrics.lapsPerMin = length(ZoneStruct.choices.isCorrect(valid))/sessmin;
ZoneStruct.metrics.p_left     = mean(ZoneStruct.choices.isLeft(valid)==1);
ZoneStruct.metrics.p_correct  = mean(ZoneStruct.choices.isCorrect(valid)==1);

% figure; 
% set(gcf, 'Position', [680   802   641   176])
% im = [a_pxy(:,1:4)*NaN, r_pxy*-1, a_pxy(:,1:4)*NaN, a_pxy, a_pxy(:,1:4)*NaN];
% im = [im(1:4, :)*NaN; im; im(1:4, :)*NaN];
% subplot(1,2,1);
% h = imagesc(im, [-.8 .8]); colorbar
% axis image
% colormap(lbmap(100,'RedBlue'));
% set(h, 'AlphaData', ~isnan(im))
% subplot(1,2,2); hold on; 
% % a_pt = [a_pt(end) a_pt];
% % r_pt = [r_pt(end) r_pt];
% % xs = [xs(end) xs];
% plot(xs, r_pt, 'r'); 
% plot(xs, a_pt, 'b')
% set(gca, 'XTick', [-180, -90 0 90 180])
% ylim([0 1])
% drawnow

end
%%



%% %%%%%%%%%%%%%%

if debug
    str1 = {'arena', 'room'};
    str2 = {'Correct', 'Error', 'L', 'R'};
    clrs = {'g' 'm' 'b' 'r'};
    for ii = 1:2
        temp = eval('%s', str1{ii});
        x=temp.x; y=temp.y; ts=temp.timestamps/1000;
        figure;
        for jj = 1:4
            subplot(2,2,jj); hold on
            plot3(x, y, ts, 'Color', [.8, .8, .8], 'LineWidth', .25)
            a = eval(sprintf('%s_start', str2{jj}));
            b = eval(sprintf('%s_end',   str2{jj}));
            for i =1:length(a)
                inds = a(i):b(i);
                plot3(x(inds), y(inds), ts(inds), 'Color', clrs{jj}, 'LineWidth', 2)
            end
            title(sprintf('N %s = %d', str2{jj}, length(a)));
            axis([-40 40 -40 40]); set(gca, 'View', [16, 58])
        end
    end
end

end


function [starts, ends, idx] = seq_finder(sequence, zs, zs_ind1, zs_ind2)
% L_seq = [1,2,3,2,1;  1,-2,3,-2,1; ]';
idx = 0;
starts = zeros(length(zs),1);
ends   = zeros(length(zs),1);

[seq_len, num_seqs] = size(sequence);
for ind = seq_len:length(zs)-1
    seq_check = zs(ind-(seq_len-1):ind);
    for sq = 1:num_seqs
        if all(seq_check==sequence(:,sq)) % L_seq
            %         seq_check
            idx = idx+1;
            starts(idx) = zs_ind1(ind-(seq_len-1));
            ends(idx)   = zs_ind2(ind);
        end
    end
end
starts = starts(1:idx);
ends = ends(1:idx);
end

function [Correct_start, Correct_end, Error_start, Error_end, frame_type, is_path_reward] = ...
    correct_choice_eval(L_start, L_end, R_start, R_end, room_rewards, arena_rewards)
nL = length(L_start);
nR = length(R_start);
n = nL+nR;
Correct_start = zeros(n,1);
Correct_end = zeros(n,1);
Error_start = zeros(n,1);
Error_end = zeros(n,1);
frame_type = zeros(n,1);
is_path_reward = false(length(room_rewards),1);
c_ind=0;
e_ind=0;

for ind = 1:nL
    a = L_start(ind);
    b = L_end(ind);
    rew_r = room_rewards>=a & room_rewards<=b;
    rew_a = arena_rewards>=a & arena_rewards<=b;
    if any(rew_r) || any(rew_a)
        c_ind = c_ind+1;
        if any(rew_r)
            frame_type(c_ind) = 1;
            is_path_reward(rew_r) = true;
        else
            frame_type(c_ind) = 2;
        end
        Correct_start(c_ind) = a;
        Correct_end(c_ind) = b;
    else
        e_ind = e_ind+1;
        Error_start(e_ind) = a;
        Error_end(e_ind) = b;
    end
    
end
for ind = 1:nR
    a = R_start(ind);
    b = R_end(ind);
    rew_r = room_rewards>=a & room_rewards<=b;
    rew_a = arena_rewards>=a & arena_rewards<=b;
    if any(rew_r) || any(rew_a)
        c_ind = c_ind+1;
        if any(rew_r)
            frame_type(c_ind) = 1;
            is_path_reward(rew_r) = true;
        else
            frame_type(c_ind) = 2;
        end
        Correct_start(c_ind) = a;
        Correct_end(c_ind) = b;
    else
        e_ind = e_ind+1;
        Error_start(e_ind) = a;
        Error_end(e_ind) = b;
    end
end
%%%%%%%%%%%
Correct_start = Correct_start(1:c_ind);
Correct_end = Correct_end(1:c_ind);
Error_start = Error_start(1:e_ind);
Error_end  = Error_end(1:e_ind);
frame_type = frame_type(1:c_ind);
%%%%%%%%%%%
end

function [pxy, pt, t, mean_ang] = path_visit_occupancy(x, y, xy_bins, ang_bins, rho_thresh, L_start, L_end, R_start, R_end)
[t,r] = cart2pol(x,y);
x(r<=rho_thresh) = NaN;
y(r<=rho_thresh) = NaN;
t(r<=rho_thresh) = NaN;
mean_ang = resultantAngle(t);
pxy = zeros(length(xy_bins)-1, length(xy_bins)-1);
pt  = zeros(1, length(ang_bins)-1);
for ii = 1:length(L_start)
    ispath = L_start(ii):L_end(ii);
    h = histcounts2(y(ispath), x(ispath), xy_bins, xy_bins);
    pxy = pxy + (h>0);
    h = histcounts(t(ispath), ang_bins);
    pt = pt + (h>0);
end
for ii = 1:length(R_start)
    ispath = R_start(ii):R_end(ii);
    h = histcounts2(y(ispath), x(ispath), xy_bins, xy_bins);
    pxy = pxy + (h>0);
    h = histcounts(t(ispath), ang_bins);
    pt = pt + (h>0);
end
pt = pt/(length(L_start) + length(R_start));
pxy = pxy/(length(L_start) + length(R_start));
pxy(pxy(:)==0) = NaN;
end

function [rew_diff, rew_diff_idx] = Tjunction_difficulty(x, y, ts, rr, zone_sequence, R_start, R_end, plotting)

% x = room.x;  y = room.y;  ts = room.timestamps./1000;
[t,r] = cart2pol(x,y);
d_t = abs(diff(unwrap(t))); d_t = [d_t(1); d_t];
d_t = movmedian(d_t, 11);
dd_t = (diff(d_t)); dd_t = [dd_t(1); dd_t];
% d_r = abs(diff(r)); d_r = [d_r(1); d_r];
% d_r = movmedian(d_r, 11);
% dd_r = (diff(d_r)); dd_r = [dd_r(1); dd_r];

if plotting
    % figure;
subplot(121); hold on; 
plot3(x, y, ts, 'k-')
end
% choice_idx = 0;
nP = length(R_start);
rew_diff = zeros(nP,1);
rew_diff_idx = zeros(nP,1);
for choice_idx = 1:nP
    % path index
    idx = R_start(choice_idx):R_end(choice_idx);
    % zone sequence, regardless of arm choice
    zs_sub    = abs(zone_sequence(idx));
    % goes from center (1) to choice (2)
    one_two   = find(zs_sub(1:end-1)==1 & zs_sub(2:end)==2, 1, 'first');
    % goes from choice (2) to arm (3)
    two_three = find(zs_sub(one_two:end-1)==2 & zs_sub(one_two+1:end)==3) + one_two;
    % full behavior indices
    tjunc_idx = R_start(choice_idx)+one_two:R_start(choice_idx)+two_three;
    if plotting; plot3(x(tjunc_idx), y(tjunc_idx), ts(tjunc_idx), 'r.-'); end
    % angular distance of animal with goal zone
    adiff = angular_distance(t(tjunc_idx), rr(tjunc_idx));
    % angular acceleration
    d_t1 = dd_t(tjunc_idx); % d_r1 = dd_r(tjunc_idx);
    % max consecutive positive samples of angular acceleration
    ispositive = double( d_t1(1:end-1)>0 & d_t1(2:end)>0);
    tjunc_midpoint = round(length(tjunc_idx)/2);
    %     ispos(1:tjunc_midpoint) = 0;
    % find number of consecutive positive samples
    for iii = 1:length(ispositive)
        if ispositive(iii)>0
            next0 = find(ispositive(iii:end)==0, 1);
            if isempty(next0); next0 = length(ispositive) - iii; end
            ispositive(iii) = next0;
        end
    end
    if max(ispositive)>=5
        idx = find(max(ispositive) == ispositive);
    else
        idx = tjunc_midpoint;
    end
    rew_diff(choice_idx) = nanmedian( adiff(idx) );
    rew_diff_idx(choice_idx) = round(mean(idx)) + tjunc_idx(1);
end
end