function [ZoneStruct] = DUAL8MAZE_performance_eval(room, arena, params, debug)
% debug = true;


%% define trajectories along the arena
% x=arena.rawx; y=arena.rawy; ts=arena.timestamps/1000;
[room.pol_theta, room.pol_rho] = cart2pol(room.x, room.y);
[arena.pol_theta, arena.pol_rho] = cart2pol(arena.x, arena.y);

x=arena.x; y=arena.y; ts=arena.timestamps/1000;
% x=room.x; y=room.y; ts=room.timestamps/1000;
t = arena.pol_theta -pi/2;
t(t<-pi)=t(t<-pi)+2*pi;
t = (t)*180/pi;
rho = arena.pol_rho;

angle_cut = 15; % cutoff angle to define choice and arm sections
rho_cut   = 20; % cutoff distance between choice and center
center_zone = rho <= rho_cut;
top_choice_zone = (rho >rho_cut) & (abs(t) <= angle_cut);
bot_choice_zone = (rho >rho_cut) & (abs(t) >= 180-angle_cut);
either_choice_zone = top_choice_zone | bot_choice_zone;
left_zone  = (rho > rho_cut) & (t < -angle_cut) & (t > -180+angle_cut);
right_zone = (rho > rho_cut) & (t >  angle_cut) & (t <  180-angle_cut);

% find entries into the room reinforcement zone
z1 = mod(room.DAT_fileinfo.RoomTrackReinforcedSector(1)-room.DAT_fileinfo.RoomTrackReinforcedSector(2)/2, 360);
z2 = mod(room.DAT_fileinfo.RoomTrackReinforcedSector(1)+room.DAT_fileinfo.RoomTrackReinforcedSector(2)/2, 360);
d = mod(rad2deg(room.pol_theta+pi)+180, 360);
in_zone_theta =  ( d >= z1 )& ( d <= z2 );
r1 = room.DAT_fileinfo.RoomTrackReinforcedSector(3);
% r2 = room.DAT_fileinfo.RoomTrackReinforcedSector(4);
r = room.pol_rho*room.DAT_fileinfo.TrackerResolution_PixPerCM;
in_zone_rho =  ( r >= r1 );% & ( r <=  r2);

in_real_zone = in_zone_theta & in_zone_rho;

% find entries into the a fake zone on the opposite area
z1 = mod(room.DAT_fileinfo.RoomTrackReinforcedSector(1)-room.DAT_fileinfo.RoomTrackReinforcedSector(2)/2, 360);
z2 = mod(room.DAT_fileinfo.RoomTrackReinforcedSector(1)+room.DAT_fileinfo.RoomTrackReinforcedSector(2)/2, 360);
faked = mod(rad2deg(room.pol_theta+pi+pi)+180, 360);
in_zone_theta =  ( faked >= z1 )& ( faked <= z2 );

in_fake_zone = in_zone_theta & in_zone_rho;


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
end

zone_sequence = arena.x*0;
zone_change = arena.x*0;
zone_sequence(center_zone) = 1;
zone_sequence(either_choice_zone) = 2;
zone_sequence(right_zone) = 3;
zone_sequence(left_zone) = 4;
zone_change(2:end) = zone_sequence(2:end) ~= zone_sequence(1:end-1);

if debug
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
end

if debug
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
%%

zs = zone_sequence(zone_change==true);
zs_ind = find(zone_change==true);
L_seq = [1,2,3,2,1]';
R_seq = [1,2,4,2,1]';
idxR = 0;
idxL = 0;
L_start = zeros(500,1);
L_end = zeros(500,1);
R_start = zeros(500,1);
R_end = zeros(500,1);

for ind = 5:length(zs)-1
    seq_check = zs(ind-4:ind);
    if all(seq_check==L_seq) % L_seq
        %         seq_check
        idxL = idxL+1;
        L_start(idxL) = zs_ind(ind-3);
        L_end(idxL)   = zs_ind(ind+1);
    elseif all(seq_check==R_seq) % R_seq
        %         seq_check
        idxR = idxR+1;
        R_start(idxR) = zs_ind(ind-3);
        R_end(idxR)   = zs_ind(ind+1);
    end
end
L_start = L_start(1:idxL);
L_end = L_end(1:idxL);
R_start = R_start(1:idxR);
R_end = R_end(1:idxR);

Correct_start = zeros(500,1);
Correct_end = zeros(500,1);
Error_start = zeros(500,1);
Error_end = zeros(500,1);

c_ind=0;
e_ind=0;
% figure;
% x=room.x; y=room.y; ts=room.timestamps/1000;
% hold on;
% subplot(121); plot3(x,y, ts)
% subplot(122); plot3(x,y, ts)
ts=room.timestamps/1000;
for ind = 1:idxL
    a = L_start(ind);
    b = L_end(ind);
    rew = room.entrance_start_idx>=a & room.entrance_start_idx<=b;
    if any(rew)
        c_ind = c_ind+1;
        Correct_start(c_ind) = a;
        Correct_end(c_ind) = b;
    else
        e_ind = e_ind+1;
        Error_start(e_ind) = a;
        Error_end(e_ind) = b;
    end
end
for ind = 1:idxR
    a = R_start(ind);
    b = R_end(ind);
    rew = room.entrance_start_idx>=a & room.entrance_start_idx<=b;
    if any(rew)
        c_ind = c_ind+1;
        Correct_start(c_ind) = a;
        Correct_end(c_ind) = b;
    else
        e_ind = e_ind+1;
        Error_start(e_ind) = a;
        Error_end(e_ind) = b;
    end
end

Correct_start = Correct_start(1:c_ind);
Correct_end = Correct_end(1:c_ind);
Error_start = Error_start(1:e_ind);
Error_end = Error_end(1:e_ind);


%%
isCorrect = zeros(length(room.x), 1);
isLeft = zeros(length(room.x), 1);
isCorrect(Correct_start)    =  1;
isCorrect(Error_start)      = -1;
isLeft(L_start)             =  1;
isLeft(R_start)             = -1;
isIDX                       = find(isCorrect~=0);
isCorrect                   = isCorrect(isIDX);
isLeft                      = isLeft(isIDX);
isTime                      = arena.timestamps(isIDX)./1000;

streaks     = isCorrect*0;
streaks_err = isCorrect*0;
for ii = 1:length(isCorrect)
    idx = find(isCorrect(1:ii)==-1, 1, 'last');
    if ~isempty(idx)
        streaks(ii) = sum(isCorrect(idx:ii)==1);
    end
    idx = find(isCorrect(1:ii)==1, 1, 'last');
    if ~isempty(idx)
        streaks_err(ii) = sum(isCorrect(idx:ii)==-1);
    end
end



choices = [];
choices.L_start_end        = [L_start, L_end];
choices.R_start_end        = [R_start, R_end];
choices.Correct_start_end  = [Correct_start, Correct_end];
choices.Error_start_end    = [Error_start, Error_end];
choices.isLeft    = isLeft;
choices.isCorrect = isCorrect;
choices.isIDX     = isIDX ;
choices.isTime    = isTime;
choices.streaks   = streaks;
choices.streaks_err = streaks_err;



if debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure; 
    subplot(121)
    hold on;
    % correct_runav = cumsum(isCorrect)/sum(abs(isCorrect));
    % left_runav = cumsum(isLeft)/sum(abs(isLeft));
    correct_runav = cumsum(isCorrect)/length(abs(isCorrect));
    left_runav = cumsum(isLeft)/length(abs(isLeft));
    plot(isTime/60, correct_runav, 'g-', 'Color', [.4 1 .4]/2);
    plot(isTime/60, left_runav, 'b--')
    plot([0 ts(end)]/60, [0 0], 'k:')
    plot([0 ts(end)]/60, [0 1], 'k:')
    plot([0 ts(end)]/60, [0 -1], 'k:')
    ylim([-1 1])
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
    for ii = 1:length(Correct_start)
        idx = Correct_start(ii):Correct_end(ii);
        plot3(x(idx) ,y(idx),ts(idx), 'g.')
    end
end

lose1_win2  = isCorrect(1:end-1)==-1 & isCorrect(2:end)==1;
lose1_lose2 = isCorrect(1:end-1)==-1 & isCorrect(2:end)==-1;
win1_win2   = isCorrect(1:end-1)==1 & isCorrect(2:end)==1;
win1_lose2  = isCorrect(1:end-1)==1 & isCorrect(2:end)==-1;
choice_perf = zeros(2,2);
choice_perf(1,1) = sum(win1_win2)/length(isCorrect);
choice_perf(1,2) = sum(win1_lose2)/length(isCorrect);
choice_perf(2,2) = sum(lose1_lose2)/length(isCorrect);
choice_perf(2,1) = sum(lose1_win2)/length(isCorrect);

lr_swap = isLeft(1:end-1)==1 & isLeft(2:end)==-1;
rl_swap = isLeft(2:end)==1   & isLeft(1:end-1)==-1;
dir_swap = [0; lr_swap | rl_swap];
dir_swap_idx = find(dir_swap==1);

lose1_win2  = isCorrect(dir_swap_idx-1)==-1 & isCorrect(dir_swap_idx)==1;
lose1_lose2 = isCorrect(dir_swap_idx-1)==-1 & isCorrect(dir_swap_idx)==-1;
win1_win2   = isCorrect(dir_swap_idx-1)==1 & isCorrect(dir_swap_idx)==1;
win1_lose2  = isCorrect(dir_swap_idx-1)==1 & isCorrect(dir_swap_idx)==-1;
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
time_to_retrieve = NaN(length(entrance_real_idx),1);
retrieve_idx = NaN(length(entrance_real_idx),1);
time_to_retrieve_fake = NaN(length(entrance_fake_idx),1);
fake_retrieve_idx = NaN(length(entrance_fake_idx),1);
for ind = 1:length(entrance_real_idx)%length(room.entrance_start_idx)
    %     e = room.entrance_start_idx(ind);
    e = entrance_real_idx(ind);
    r = find(center_zone(e:end), 1, 'first')+e;
    if ~isempty(r)
        time_to_retrieve(ind) = ts(r) - ts(e);
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
ZoneStruct.arena_zones.left_start   = L_start;
ZoneStruct.arena_zones.left_end   = L_end;
ZoneStruct.arena_zones.right_start   = R_start;
ZoneStruct.arena_zones.right_end   = R_end;
ZoneStruct.metrics.Leftarm_bias   = idxL/(idxL+idxR);
ZoneStruct.metrics.sequence_prob_Left   = prob_L;
ZoneStruct.metrics.sequence_prob_Correct   = prob_C;
ZoneStruct.metrics.sequence_prob_times   = t_inds(2:end);
ZoneStruct.metrics.choice_performance      = choice_perf;
ZoneStruct.metrics.swap_performance        = swap_perf;
ZoneStruct.metrics.choice_swap_key  = {'win-win', 'win-lose'; 'lose-win', 'lose-lose'};
ZoneStruct.metrics.average_streak = mean(streaks(streaks>0));
ZoneStruct.metrics.average_streak_err = mean(streaks_err(streaks_err>0));
ZoneStruct.choices = choices;

if isfield(room.DAT_fileinfo, 'RoomTrackReinforcedSector')
    ZoneStruct.room_reward.reinforcedenter_idx   = entrance_real_idx;
    ZoneStruct.room_reward.opposite_enter_idx    = entrance_fake_idx;
    ZoneStruct.room_reward.reinforced_zone_pref  = pref;
    ZoneStruct.room_reward.time_to_retrieve      = time_to_retrieve;
    ZoneStruct.room_reward.time_to_retrieve_fake = time_to_retrieve_fake;
    
    ZoneStruct.metrics.reinforced_zone_preference   = pref;
    ZoneStruct.metrics.mean_retrieval_time          = nanmean(time_to_retrieve);
    ZoneStruct.metrics.mean_retrieval_time_fake     = nanmean(time_to_retrieve_fake);
end
if isfield(arena.DAT_fileinfo, 'ArenaTrackReinforcedSector')
    error('Arena not yet implemented')
    % ZoneStruct.room_reward.reinforcedenter_idx   = entrance_real_idx;
    % ZoneStruct.room_reward.opposite_enter_idx    = entrance_fake_idx;
    % ZoneStruct.room_reward.reinforced_zone_pref  = pref;
    % ZoneStruct.room_reward.time_to_retrieve      = time_to_retrieve;
    % ZoneStruct.room_reward.time_to_retrieve_fake = time_to_retrieve_fake;
end

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
