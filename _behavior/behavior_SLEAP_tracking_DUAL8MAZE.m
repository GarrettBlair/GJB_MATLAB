function [room, arena, params, sleap] = behavior_SLEAP_tracking_DUAL8MAZE(room_tracking_fname, arena_tracking_fname, sleap_fname, sleap_timestamp_fname, params)
% global params
debug = false;

[room , room_params]  = read_APA_csv(room_tracking_fname, []);
[arena, arena_params] = read_APA_csv(arena_tracking_fname, []);
%
arena_radius_cm = 100*arena.DAT_fileinfo.ArenaDiameter_m/2;

if isfield(params, 'camera_positions')
    pos_center_X = params.camera_positions.pos_center_X;
    pos_center_Y = params.camera_positions.pos_center_Y;
    arena_size_px = params.camera_positions.arena_size_px;
else
    [pos_center_X, pos_center_Y, arena_size_px] = get_pos_centers_arena_room(arena, room, arena_radius_cm);
    params.camera_positions.pos_center_X = pos_center_X;
    params.camera_positions.pos_center_Y = pos_center_Y;
    params.camera_positions.arena_size_px = arena_size_px;
end
% video_center = [480,640]/2;
% video_pixpercm = 370/80;
%
ax=arena.rawx; ay=arena.rawy; ats=arena.timestamps;
rx=room.rawx;  ry=room.rawy;  rts=room.timestamps;

%     ax = ax(ismember(ats,rts)); ay = ay(ismember(ats,rts));
%     rx = rx(ismember(rts,ats)); ry = ry(ismember(rts,ats));
goods = ax>0 & ay>0 & rx>0 & ry>0;
ax = arena_radius_cm*(ax-pos_center_X)/arena_size_px; % -room.DAT_fileinfo.ArenaCenterXY(1);
ay = arena_radius_cm*(ay-pos_center_Y)/arena_size_px; % -room.DAT_fileinfo.ArenaCenterXY(2);
rx = arena_radius_cm*(rx-pos_center_X)/arena_size_px; % -room.DAT_fileinfo.ArenaCenterXY(1);
ry = arena_radius_cm*(ry-pos_center_Y)/arena_size_px; % -room.DAT_fileinfo.ArenaCenterXY(2);

angdiff1 = compute_arena_angle(ax, ay, rx, ry, rts, goods, 15, false);
%


if ~isempty(sleap_fname) % isfile(sleap_fname)== true && 

% %     figure; hold on
% %     plot3(ax, ay, ats/1000, '.', 'Color', [.4, .4, 1], 'LineWidth', .25)
% %     plot3(rx, ry, rts/1000, '.', 'Color', [1, .4, .4], 'LineWidth', .25)
%     ad = NaN(length(ax),1);
%     % d(goods) = atan2(ay(goods)-ry(goods), ax(goods)-rx(goods));
%     ad(goods) = atan2(ay(goods), ax(goods));
%     rd = NaN(length(ax),1);
%     % d(goods) = atan2(ay(goods)-ry(goods), ax(goods)-rx(goods));
%     rd(goods) = atan2(ry(goods), rx(goods));
%     
%     angdiff1 = ad-rd; % mod(ad-rd, pi);
%     angdiff1 = unwrap(angdiff1, [pi]);
%     angdiffspeed = abs(diff(angdiff1));
%     counter=0;
%     while any(angdiffspeed>.1) && counter<1000 % remove large jumps
%         counter = counter+1;
%         angdiffspeed = abs(diff(angdiff1));
%         angdiff1(angdiffspeed>.1)=NaN;
%         nanind = find(isnan(angdiff1));
%         angdiff1(nanind) = interp1(find(~isnan(angdiff1)), angdiff1(~isnan(angdiff1)), nanind, 'linear');
%     end
%     % s=15;
%     % angdiff1(s/2:end-s/2) = conv(angdiff1, ones(s,1)/s,'valid');
% %     angdiff1 = mod(angdiff1+pi, -2*pi)+pi;
% %     figure; plot(ats/1000, angdiff1)
%     %
    varnames = {'node_names' 'edge_names' 'edge_inds' 'tracks' 'track_occupancy'...
        'point_scores' 'instance_scores' 'labels_path' 'video_path' 'video_ind'};
    sleap = [];
    for ii = 1:length(varnames)
        structstr = sprintf('sleap.%s = h5read(sleap_fname, ''/%s'');', varnames{ii}, varnames{ii});
        eval(structstr)
    end
    sleap_ts = csvread(sleap_timestamp_fname);
    nsamples_sleap = size(sleap.tracks,1);
    sleap_ts = sleap_ts(1:nsamples_sleap);
    angdiff_sleap = interp1(rts, angdiff1, sleap_ts, 'linear');
    ledx = 580*ones(nsamples_sleap,1); ledy =240*ones(nsamples_sleap,1);
    
    valid = 1:nsamples_sleap; % length(room.timestamps);
% %     if length(room.timestamps) ~= length(sleap_ts)
% %         if abs(length(room.timestamps) - length(sleap_ts)) <= 1
% %             if length(room.timestamps) < length(sleap_ts)
% %                 valid = 1:length(room.timestamps);
% %             else
% %                 valid = 1:length(sleap_ts);
% % %                 error('file size difference, unhandled exception, must be equal')
% %             end
% %         else
% %             error('file size difference, unhandled exception, must be equal')
% %         end
% %     end
%     fn = sleap.video_path{1};
%     fn_ts = [fn(1:strfind(fn, '.avi')-1) '_timeStamps.csv'];
%     sleap_ts = dlmread(fn_ts);
    if isfield(params.camera_positions, 'sleap_pos_center_X')
        sleap_pos_center_X  = params.camera_positions.sleap_pos_center_X;
        sleap_pos_center_Y  = params.camera_positions.sleap_pos_center_Y;
        sleap_arena_size_px = params.camera_positions.sleap_arena_size_px;
    else
        [sleap_pos_center_X, sleap_pos_center_Y, sleap_arena_size_px] = get_pos_centers_SLEAP(sleap, arena_radius_cm);
        params.camera_positions.sleap_pos_center_X  = sleap_pos_center_X;
        params.camera_positions.sleap_pos_center_Y  = sleap_pos_center_Y;
        params.camera_positions.sleap_arena_size_px = sleap_arena_size_px;
    end

    sleap = correct_sleap_jumps(sleap);
    sleap.timeStamps = sleap_ts;
    body_idx = find(contains(sleap.node_names, 'body_center'));
    head_idx = find(contains(sleap.node_names, 'head_base'));
    tail_idx = find(contains(sleap.node_names, 'tail_base'));
%     mxy = squeeze(mean(sleap.tracks,2));
%     mxy = squeeze(sleap.tracks(:,[head_idx, body_idx, tail_idx],:));
    mxy = squeeze(sleap.tracks(:,[body_idx],:));
    mxy_score = sleap.instance_scores/length(sleap.node_names); % squeeze(min(sleap.point_scores(:,[head_idx, body_idx, tail_idx],:), [], 2));
%     mxy(mxy_score<=.5) = NaN;
    room.rawx = mxy(valid,1);
    room.rawy = mxy(valid,2);
    room.sleap = sleap;
    
    [at, arho] = cart2pol(mxy(valid,1)-sleap_pos_center_X, mxy(valid,2)-sleap_pos_center_Y); 
    
    at = at+angdiff_sleap(valid);
    [ax2, ay2] = pol2cart(at, arho);
    arena.rawx = ax2+sleap_pos_center_X;
    arena.rawy = ay2+sleap_pos_center_Y;
    
    [at, arho] = cart2pol(ledx(valid)-sleap_pos_center_X, ledy(valid)-sleap_pos_center_Y); 
    at = at+angdiff_sleap(valid);
    [ax2, ay2] = pol2cart(at, arho);
    ledx_room = ax2+sleap_pos_center_X;
    ledy_room = ay2+sleap_pos_center_Y;
    
    room.dat_timestamps  = room.timestamps; %    room.timestamps(valid);
    arena.dat_timestamps  = arena.timestamps; %    room.timestamps(valid);
    room.timestamps  = sleap.timeStamps; %    room.timestamps(valid);
    arena.timestamps = sleap.timeStamps; % arena.timestamps(valid);
%     arena.angular_offset = NaN(length(arena.timestamps), 1);
    arena.angular_offset = angdiff_sleap(valid);
    centerX = params.camera_positions.sleap_pos_center_X; 
    centerY = params.camera_positions.sleap_pos_center_Y;
    arena_size_px = params.camera_positions.sleap_arena_size_px;
else
    warning('No SLEAP, Using DAT files: \n\t%s...\n\t%s...', room_tracking_fname, arena_tracking_fname)
    centerX = 0; centerY = 0;
    arena_radius_cm = 100*arena.DAT_fileinfo.ArenaDiameter_m/2;
    
    centerX = params.camera_positions.pos_center_X;
    centerY = params.camera_positions.pos_center_Y;
    arena_size_px = params.camera_positions.arena_size_px;
    params.camera_positions.sleap_pos_center_X  = 0;
    params.camera_positions.sleap_pos_center_Y  = 0;
    params.camera_positions.sleap_arena_size_px = 0;
    room.dat_timestamps  = room.timestamps; %    room.timestamps(valid);
    arena.dat_timestamps  = arena.timestamps; %    room.timestamps(valid);
    
    arena.angular_offset = angdiff1;
%     [room , room_params] = read_APA_csv(room_tracking_fname, []);
%     [arena, arena_params] = read_APA_csv(arena_tracking_fname, []);
end


if contains(arena_tracking_fname, '_HC')
    warning('Homecage Arena file, no rotation, setting room and arena structures ==')
    arena = room;
    
% elseif contains(arena_tracking_fname, 'DAT Files\NAPA_24510_IL1_')
%     % this arena file was over
%     warning('weird file fix pos')
end


if false % debug
figure; 
subplot(121)
hold on
x=arena.rawx; y=arena.rawy; ts=arena.timestamps;
plot3(x, y, ts/1000, '-', 'Color', [.4, .4, 1], 'LineWidth', .25)
try 
plot3(ledx, ledy, ts/1000, '-.', 'Color', [.4, .4, .6], 'LineWidth', .25)
end
% plot3(ledx_room, ledy_room, ts/1000, '-.', 'Color', [.6, .4, .4], 'LineWidth', .25)
subplot(122)
hold on
x=room.rawx; y=room.rawy; ts=room.timestamps;
plot3(x, y, ts/1000, '-', 'Color', [1, .4, .4], 'LineWidth', .25)

% plot3(ledx, ledy, ts/1000, '-.', 'Color', [.4, .4, .6], 'LineWidth', .25)
try
    plot3(ledx_room, ledy_room, ts/1000, '-.', 'Color', [.6, .4, .4], 'LineWidth', .25)
end% axis([-100 600 -100 600 0 100])
end



%%
room.params = params;
room.params.headerSize = room_params.headerSize;
arena.params = params;
arena.params.headerSize = arena_params.headerSize;

ksize = round(params.behav_fps*params.behav_smoothing_interval);
kern = ones(ksize, 1); kern = kern./sum(kern(:));

% remove any internal reflections in tracking
testx = arena.rawx; testx(testx==0)=NaN;
testy = arena.rawy; testy(testy==0)=NaN;
% center it
testx = params.arena_radius*(testx-centerX)./arena_size_px; % params.pixpercm;
testy = params.arena_radius*(testy-centerY)./arena_size_px; % params.pixpercm;

bads = define_valid_arenapos(testx, testy);
if sum(bads)/length(bads) >= .1*length(bads)
    error('BAD TRACKING')
elseif false % debug == true
    figure; 
    subplot(121); hold on
    plot(arena.rawx, arena.rawy, 'r.')
    subplot(122); hold on
    plot(room.rawx, room.rawy, 'r.')
    arena.rawx(bads) = NaN;
    arena.rawy(bads) = NaN;
    room.rawx(bads)  = NaN;
    room.rawy(bads)  = NaN;
    subplot(121); hold on
    plot(arena.rawx, arena.rawy, 'k.')
    subplot(122); hold on
    plot(room.rawx, room.rawy, 'k.')
    
    figure(8); clf; 
    subplot(111); hold on
    t = round(arena.timestamps/1000);
    plot3(arena.rawx, arena.rawy, t, 'b.-')
    plot3(room.rawx, room.rawy, t, 'r.-')
    subplot(122); hold on
    plot(room.rawx, room.rawy, 'r.')
    arena.rawx(bads) = NaN;
    arena.rawy(bads) = NaN;
    room.rawx(bads)  = NaN;
    room.rawy(bads)  = NaN;
    subplot(121); hold on
    plot(arena.rawx, arena.rawy, 'k.')
    subplot(122); hold on
    plot(room.rawx, room.rawy, 'k.')
end
%%
strnames = {'room', 'arena'};
for strLoop = 1:length(strnames)
    %%
    % input structure
    eval(sprintf('behav = %s;', strnames{strLoop}));
    nanind = behav.rawx==0 & behav.rawy == 0;
    if ~any(nanind==0) && strcmp(strnames{strLoop}, 'arena')==1
        behav.x = interp1(room.timestamps, room.x, behav.timestamps, 'linear');
        behav.y = interp1(room.timestamps, room.y, behav.timestamps, 'linear');
        x = behav.x;
        y = behav.y;
        behav.wasnan = nanind;

    else
        %%
        t = behav.timestamps;
        if length(unique(t)) ~= length(t)
            ind = find(diff(t)==0)+1;
            t(ind) = t(ind)+1;
        end
        ts = t./1000; % convert from ms to seconds
        x = behav.rawx;
        y = behav.rawy;
        if isfield(behav, 'Angle')
            a = behav.Angle;
            behav.raw_headangle = a;
            headnan = a==-1;
        end
        x(nanind) = NaN;
        y(nanind) = NaN;
        
        smooth_jumps=false;
        if smooth_jumps % smooth_jumps  % not quite yet 
            figure; plot(x,y)
            jump_thresh   = 10;
            return_thresh = 30;
            twin_seconds = 15;
            [x,y] = smooth_jumps_apa(x, y, ts, params.pixpercm, jump_thresh, return_thresh, twin_seconds);
            hold on; plot(x,y);
        end
        
        %
        if params.nan_interp && any(~nanind)
            nanind = (isnan(x) & isnan(y));
            xn = interp1(ts(~nanind), x(~nanind), ts(nanind), 'linear');
            yn = interp1(ts(~nanind), y(~nanind), ts(nanind), 'linear');
            x(nanind) = xn; y(nanind) = yn;
            if isfield(behav, 'Angle') && any(~headnan)
                a(headnan) = interp1(ts(~headnan), a(~headnan), ts(headnan), 'linear');
            end
            if any(isnan([x(1) x(end) y(1) y(end)]))
                nanind = isnan(x.*y);
                xn = interp1(ts(~nanind), x(~nanind), ts(nanind), 'nearest', 'extrap');
                yn = interp1(ts(~nanind), y(~nanind), ts(nanind), 'nearest', 'extrap');
                x(nanind) = xn; y(nanind) = yn;
            end
        end
        behav.wasnan = nanind;
        
%         x = (x-params.arena_center(1))./params.pixpercm;
%         y = (y-params.arena_center(2))./params.pixpercm;
%         x = arena_radius_cm*(x-sleap_pos_center_X)./sleap_arena_size_px;
%         y = arena_radius_cm*(y-sleap_pos_center_Y)./sleap_arena_size_px;
        x = arena_radius_cm*(x-centerX)./arena_size_px;
        y = arena_radius_cm*(y-centerY)./arena_size_px;

        xn = conv(x, kern, 'same');
        yn = conv(y, kern, 'same');
%         tails = [1:floor(ksize/2)+1, length(x)-floor(ksize/2)-1:length(x)];
        xn(1:floor(ksize/2)+1) = xn(floor(ksize/2)+2);
        xn(length(x)-floor(ksize/2):length(x)) = xn(length(x)-floor(ksize/2)-1);
        yn(1:floor(ksize/2)+1) = yn(floor(ksize/2)+2);
        yn(length(x)-floor(ksize/2):length(x)) = yn(length(x)-floor(ksize/2)-1);
%         xn(tails) = x(tails);
%         yn(tails) = y(tails);
        x = xn;
        y = yn*-1; % flip y axis to align with cpu view
        
        behav.x = x; behav.y = y;
        if isfield(behav, 'Angle')
            a = deg2rad(a);
            an = (unwrap(a));
            an = conv(an, kern, 'same');
            an(1:floor(ksize/2)+1) = an(floor(ksize/2)+2);
            an(length(x)-floor(ksize/2):length(x)) = an(length(x)-floor(ksize/2)-1);
            an = wrapTo2Pi(an);
%             an(an>360) = 360; an(an<0) = 0;
            a = an;
            behav.Angle = rad2deg(a);
        end
    end
    [behav.pol_theta, behav.pol_rho] = cart2pol(x,y);
    
    behav_dt = [median(diff(behav.timestamps)); diff([behav.timestamps])]/1000;
    behav.dt = behav_dt;
    if params.correct_dt
        % sometimes the camera will disconnect and reconnect, with large jumps
        % in the timestamp file
        bad_dt_thresh = mean(behav_dt)+10*std(behav_dt);
        bad_vals = behav_dt >= bad_dt_thresh;
        if ( sum(bad_vals)/length(behav_dt) )>.01
           warning('Many bad values found in timestamp dt, should check!') 
        end
        behav_dt(bad_vals) = median(behav_dt);
    end
    behav.dt_corrected = behav_dt;
    behav.speed = sqrt(diff([x(1); x]).^2 + diff([y(1); y]).^2)./behav_dt;
    behav.speed(1) = behav.speed(2);
    
    behav.entrance_start_idx = find(behav.State(1:end-1)==1 & behav.State(2:end)==2) + 1; % 1 == Entrance latency
    lostinds = find(behav.State(1:end-1)==5 & behav.State(2:end)==2) + 1; % 5 == lost tracking
    for jjj = 1:length(lostinds)
        % this ensures that inteval latency lost tracking is not counted as
        % a new entrance
        prevStates = behav.State( lostinds(jjj)-100:lostinds(jjj)-1 );
        prevStates = prevStates(prevStates~=5); % look at the previous states before loss
        if ~isempty(prevStates) && prevStates(end) == 1
            behav.entrance_start_idx(end+1) = lostinds(jjj);
        end
    end
    behav.entrance_start_idx = sort(behav.entrance_start_idx, 'ascend');
    behav.shock_start_idx = find(diff([behav.State(1); behav.State] == 2)==1); % 2 == Shock on
%     behav.shock_start_idx = find(behav.CurrentLevel(1:end-1)==0 & behav.CurrentLevel(2:end)>0 ) + 1; % 2 == Shock on
    
    behav.num_entrances = length(behav.entrance_start_idx);
    behav.num_shocks = length(behav.shock_start_idx);
    if behav.num_entrances>0
        behav.first_entrance = behav.timestamps(behav.entrance_start_idx(1))./1000;
    else
        behav.first_entrance = NaN;
    end
    
    % output structure
    eval(sprintf('%s = behav;', strnames{strLoop}));
end

if exist('ledx_room', 'var')
room.ledx = arena_radius_cm*(ledx_room-centerX)./arena_size_px;
room.ledy = arena_radius_cm*(ledy_room-centerY)./arena_size_px;
arena.ledx = arena_radius_cm*(ledx-centerX)./arena_size_px;
arena.ledy = arena_radius_cm*(ledy-centerY)./arena_size_px;
end
% debug = false;
if debug
figure(); 
subplot(121)
hold on
x=arena.x; y=arena.y; ts=arena.timestamps;
plot3(x, y, ts/1000, '-', 'Color', [.4, .4, 1], 'LineWidth', .25)
plot3(arena.ledx, arena.ledy, ts/1000, '-.', 'Color', [.4, .4, .6], 'LineWidth', .25)
subplot(122)
hold on
x=room.x; y=room.y; ts=room.timestamps;
plot3(x, y, ts/1000, '-', 'Color', [1, .4, .4], 'LineWidth', .25)
plot3(room.ledx, room.ledy, ts/1000, '-.', 'Color', [.6, .4, .4], 'LineWidth', .25)
% axis([-100 600 -100 600 0 100])
end

if exist('sleap', 'var')
    sleap_room_pos  = sleap.tracks;
    sleap_arena_pos = sleap.tracks*0;
    ts = sleap.timeStamps;
    for ii = 1:length(sleap.node_names)
        x = sleap_room_pos(:,ii,1);
        y = sleap_room_pos(:,ii,2);
        nanind = sleap.instance_scores/length(sleap.node_names) <.5;
        x(nanind) = NaN;
        y(nanind) = NaN;
        if params.nan_interp && any(~nanind)
            nanind = (isnan(x) & isnan(y));
            xn = interp1(ts(~nanind), x(~nanind), ts(nanind), 'linear');
            yn = interp1(ts(~nanind), y(~nanind), ts(nanind), 'linear');
            x(nanind) = xn; y(nanind) = yn;
            if any(isnan([x(1) x(end) y(1) y(end)]))
                nanind = isnan(x.*y);
                xn = interp1(ts(~nanind), x(~nanind), ts(nanind), 'nearest', 'extrap');
                yn = interp1(ts(~nanind), y(~nanind), ts(nanind), 'nearest', 'extrap');
                x(nanind) = xn; y(nanind) = yn;
            end
        end
        
        x = arena_radius_cm*(x-centerX)./arena_size_px;
        y = arena_radius_cm*(y-centerY)./arena_size_px;
        
        xn = conv(x, kern, 'same');
        yn = conv(y, kern, 'same');
        %         tails = [1:floor(ksize/2)+1, length(x)-floor(ksize/2)-1:length(x)];
        xn(1:floor(ksize/2)+1) = xn(floor(ksize/2)+2);
        xn(length(x)-floor(ksize/2):length(x)) = xn(length(x)-floor(ksize/2)-1);
        yn(1:floor(ksize/2)+1) = yn(floor(ksize/2)+2);
        yn(length(x)-floor(ksize/2):length(x)) = yn(length(x)-floor(ksize/2)-1);
        %         xn(tails) = x(tails);
        %         yn(tails) = y(tails);
        x = xn;
        y = yn;
        sleap_room_pos(:,ii,1) = x;
        sleap_room_pos(:,ii,2) = y*-1;
        
        [at, arho] = cart2pol(x, y);
        
        at = at+arena.angular_offset;
        [sleap_arena_pos(:,ii,1), sleap_arena_pos(:,ii,2)] = pol2cart(at, arho);
    end
    room.sleap_pos = sleap_room_pos;
    arena.sleap_pos = sleap_arena_pos;
    if false % debug
        figure
        subplot(121)
        plot3(arena.x, arena.y, arena.timestamps/1000)
        for ii = 1:length(sleap.node_names)
            hold on;
            plot3(arena.sleap_pos(:,ii,1), arena.sleap_pos(:,ii,2), arena.timestamps/1000)
        end
        legend(['pos'; sleap.node_names], 'Interpreter', 'none')
        
        subplot(122)
        plot3(room.x, room.y, room.timestamps/1000)
        for ii = 1:length(sleap.node_names)
            hold on;
            plot3(room.sleap_pos(:,ii,1), room.sleap_pos(:,ii,2), room.timestamps/1000)
        end
        legend(['pos'; sleap.node_names], 'Interpreter', 'none', 'Location', 'northeast')
    end
else
    sleap = [];
end
end

function bads = define_valid_arenapos(x,y)
% lx =[...
%    -5.4378;...
%    -3.5945;...
%    -2.8571;...
%    -2.8571;...
%    -3.5945;...
%    -5.4378;...
%   -24.0553;...
%   -32.3502;...
%   -32.3502;...
%   -16.1290];
% ly = [...
%   -30.0875;...
%   -21.9242;...
%    -4.4315;...
%     7.6968;...
%    21.9242;...
%    30.0875;...
%    21.4577;...
%     6.2974;...
%    -9.7959;...
%   -27.9883];
% rx = [...
%     6.1751;...
%     4.3318;...
%     4.8848;...
%     5.4378;...
%     6.9124;...
%    25.1613;...
%    32.5346;...
%    29.2166;...
%    18.7097;...
%     7.8341];
% ry = [...
%   -29.6210;...
%   -22.1574;...
%     7.6968;...
%    23.0904;...
%    30.5539;...
%    19.3586;...
%     2.7988;...
%   -10.0292;...
%   -24.9563;...
%   -29.8542];
rx = [7.49, 8.18, 17.17, 28.46, 33.76, 35.83, 30.76, 25.46, 16.71, 9.79, 6.11, 6.34, 7.26, 7.26, 7.03, 7.26];
ry = [-27.11, -34.11, -32.07, -23.32, -10.20, 2.62, 14.58, 23.32, 29.15, 31.20, 22.74, -4.08, -11.95, -18.95, -24.78, -26.24];
lx = rx*-1;
ly = ry;
bads_left  = inpolygon(x,y, lx, ly);
bads_right = inpolygon(x,y, rx, ry);
bads = bads_left | bads_right;

end

function [pos_center_X, pos_center_Y, arena_size_px] = get_pos_centers_arena_room(arena, room, arena_scale_cm)

ax=arena.rawx; ay=arena.rawy; ats=arena.timestamps;
rx=room.rawx;  ry=room.rawy;  rts=room.timestamps;
goods = ax>0 & ay>0 & rx>0 & ry>0;
qx1 = min([ax(goods); rx(goods)]);
qx2 = median([ax(goods); rx(goods)]);
qx3 = max([ax(goods); rx(goods)]);
qx4 = (qx3-qx1)/2 + qx1;

qy1 = min([ay(goods); ry(goods)]);
qy2 = median([ay(goods); ry(goods)]);
qy3 = max([ay(goods); ry(goods)]);
qy4 = (qy3-qy1)/2 + qy1;
q=-1;
figure(1927); clf; hold on;
set(gcf, 'Position', [192         489        1269         469])
title('click center and border, space to stop')

% a_scale = 100*room.DAT_fileinfo.ArenaDiameter_m/2;

while q~=32
    clf
    subplot(121);
    hold on
    scatter(room.rawx(goods), room.rawy(goods), 10,'r.')
    scatter(arena.rawx(goods), arena.rawy(goods), 4,'b.')
    plot([qx1 qx2 qx3], [qy1 qy2 qy3], '-om')
    plot([qx1 qx4 qx3], [qy1 qy4 qy3], '-oc')
    
    plot([qx1 qx2 qx3], [qy3 qy2 qy1], '-om')
    plot([qx1 qx4 qx3], [qy3 qy4 qy1], '-oc')
    rectangle('Position', [qx1 qy1 qx3-qx1 qy3-qy1], 'EdgeColor', 'c')
    axis([0 250 0 250])
    [px, py, ~] = ginput(1);
    scatter(px, py, 200,'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .1)
    scatter(px, py, 10,'k.')
    title('Click center of arena and edge!')
    [px2, py2, ~] = ginput(1);
    scatter(px2, py2, 200,'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .1)
    scatter(px2, py2, 10,'k.')
    
    plot([px px2], [py py2], '-k')
    a_size = sqrt((px-px2)^2 + (py-py2)^2);
    text(mean([px px2]), mean([py py2])-10, sprintf('%d pixels', round(a_size)))
    
    subplot(122);
    hold on
    tempax = arena_scale_cm*(ax-px)/a_size;
    tempay = arena_scale_cm*(ay-py)/a_size;
    temprx = arena_scale_cm*(rx-px)/a_size;
    tempry = arena_scale_cm*(ry-py)/a_size;
    scatter(temprx(goods), tempry(goods), 10,'r.')
    scatter(tempax(goods), tempay(goods), 4,'b.')
    axis([-50 50 -50 50])
    
    title('Space to confirm!')
    [~, ~, q] = ginput(1);

end

pos_center_X = round(px);
pos_center_Y = round(py);
arena_size_px = a_size;
end

function [pos_center_X, pos_center_Y, arena_size_px] = get_pos_centers_SLEAP(sleap, arena_scale_cm)
%%
mxy = squeeze(mean(sleap.tracks,2));

rx = mxy(:,1);
ry = mxy(:,2);
goods = ~isnan(rx) & ~isnan(ry);
qx1 = min([rx(goods)]);
qx2 = median([rx(goods)]);
qx3 = max([rx(goods)]);
qx4 = (qx3-qx1)/2 + qx1;

qy1 = min([ry(goods)]);
qy2 = median([ry(goods)]);
qy3 = max([ry(goods)]);
qy4 = (qy3-qy1)/2 + qy1;
q=-1;
figure(1928); clf; hold on;
set(gcf, 'Position', [192         489        1269         469])
title('click center and border, space to stop')

% a_scale = 100*room.DAT_fileinfo.ArenaDiameter_m/2;

while q~=32
    %%
    clf
    subplot(121);
    hold on
    scatter(rx(goods), ry(goods), 10,'r.')
%     scatter(arena.rawx(goods), arena.rawy(goods), 4,'b.')
    plot([qx1 qx2 qx3], [qy1 qy2 qy3], '-om')
    plot([qx1 qx4 qx3], [qy1 qy4 qy3], '-oc')
    
    plot([qx1 qx2 qx3], [qy3 qy2 qy1], '-om')
    plot([qx1 qx4 qx3], [qy3 qy4 qy1], '-oc')
    rectangle('Position', [qx1 qy1 qx3-qx1 qy3-qy1], 'EdgeColor', 'c')
    [px, py, ~] = ginput(1);
    scatter(px, py, 200,'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .1)
    scatter(px, py, 10,'k.')
    
    [px2, py2, ~] = ginput(1);
    scatter(px2, py2, 200,'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .1)
    scatter(px2, py2, 10,'k.')
    
    plot([px px2], [py py2], '-k')
    a_size = sqrt((px-px2)^2 + (py-py2)^2);
    text(mean([px px2]), mean([py py2])-10, sprintf('%d pixels', round(a_size)))
    
    subplot(122);
    hold on
%     tempax = a_scale*(ax-px)/a_size;
%     tempay = a_scale*(ay-py)/a_size;
    temprx = arena_scale_cm*(rx-px)/a_size;
    tempry = arena_scale_cm*(ry-py)/a_size;
    scatter(temprx(goods), tempry(goods), 10,'r.')
    plot([0,0], [-50 50], 'k')
    plot([-50 50], [0,0], 'k')
%     scatter(tempax(goods), tempay(goods), 4,'b.')
    axis equal square
    
    
    [~, ~, q] = ginput(1);

end

pos_center_X = round(px);
pos_center_Y = round(py);
arena_size_px = a_size;
end

function angdiff = compute_arena_angle(ax, ay, rx, ry, ts, goods, rho_thresh, debug)
%%

angdiff = ax*NaN;
[at,ar] = cart2pol(ax,ay);
[rt,rr] = cart2pol(rx,ry);
thresh = rho_thresh; % 15; %quantile([rr; ar], cuts(iii));

% goods = ax>0 & ay>0 & rx>0 & ry>0;
valid = goods & (ar>thresh & rr>thresh);
angdiff(valid) = unwrap(at(valid) - rt(valid));
angdiff(isnan(angdiff)) = interp1(ts(~isnan(angdiff)), angdiff(~isnan(angdiff)), ts(isnan(angdiff)), 'linear', NaN);
angdiff(isnan(angdiff)) = interp1(ts(~isnan(angdiff)), angdiff(~isnan(angdiff)), ts(isnan(angdiff)), 'nearest', 'extrap');

angdiffspeed = abs(diff(angdiff));
if debug==true
figure(); hold on
plot(ts/1000,  mod(angdiff, pi), '-', 'color', 'k');
end
% remove jumps
counter=0;
while any(angdiffspeed>.1) && counter<1000 % remove large jumps
    counter = counter+1;
    angdiffspeed = abs(diff(angdiff));
    angdiff(angdiffspeed>.1)=NaN;
    nanind = find(isnan(angdiff));
    angdiff(nanind) = interp1(find(~isnan(angdiff)), angdiff(~isnan(angdiff)), nanind, 'linear');
end
% smooth speedy jumps
s=150;
kern = gausswin(s); kern = kern/sum(kern(:));
angdiff(s/2:end-s/2) = conv(angdiff, kern, 'valid');

if debug==true
plot(ts/1000,  mod(angdiff, pi), '-', 'color', 'b');
end
angdiff = mod(angdiff, 2*pi);
end
