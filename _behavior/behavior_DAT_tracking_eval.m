function [room, arena, params] = behavior_DAT_tracking_eval(room_tracking_fname, arena_tracking_fname, params)
% global params
[room , room_params] = read_APA_csv(room_tracking_fname, []);
[arena, arena_params] = read_APA_csv(arena_tracking_fname, []);
if contains(arena_tracking_fname, '_HC')
    warning('Homecage Arena file, no rotation, setting room and arena structures ==')
    arena = room;
end

room.params = params;
room.params.headerSize = room_params.headerSize;
arena.params = params;
arena.params.headerSize = arena_params.headerSize;

ksize = round(params.behav_fps*params.behav_smoothing_interval);
kern = ones(ksize, 1); kern = kern./sum(kern(:));

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
        %
        x = behav.rawx;
        y = behav.rawy;
        if isfield(behav, 'Angle')
            a = behav.Angle;
            behav.raw_headangle = a;
            headnan = a==-1;
        end
        x(nanind) = NaN;
        y(nanind) = NaN;
        t = behav.timestamps;
        if length(unique(t)) ~= length(t)
            ind = find(diff(t)==0)+1;
            t(ind) = t(ind)+1;
        end
        ts = t./1000; % convert from ms to seconds
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
        
        x = (x-params.arena_center(1))./params.pixpercm;
        y = (y-params.arena_center(2))./params.pixpercm;
    % % %     x = movmedian(x, round(params.behav_fps/2));
    % % %     y = movmedian(y, round(params.behav_fps/2));
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
%%

%%
end
