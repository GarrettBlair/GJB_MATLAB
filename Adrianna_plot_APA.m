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
experiment_folder = 'D:\ALL_DAT\';
animals = {'Hipp18240' 'Acc19947' 'Acc20832' 'HPCACC24500' 'HPCACC24501' 'HPCACC24502' 'HPCACC24503'};
animaltype = {'HPC' 'ACC' 'ACC' 'HPCACC' 'HPCACC' 'HPCACC' 'HPCACC'};
maxtr = 50;%length(trials);
na = length(animals);
name_mat       = cell(maxtr, na);
trial_mat      = cell(maxtr, na);
trial_ind      = cell(maxtr, na);
trial_type      = cell(maxtr, na);
relative_daytime        = NaN(maxtr, na);
sessName       = cell(maxtr, na);
numShocks      = NaN(maxtr, na);
numEntr        = NaN(maxtr, na);
maxtimeAvoid   = NaN(maxtr, na);
fistEntr       = NaN(maxtr, na);
occ_diff_room  = NaN(maxtr, na);
occ_diff_arena = NaN(maxtr, na);
sessMinutes    = NaN(maxtr, na);
time_origin = datetime(2020, 1, 1, 0, 0, 0); % 
behav_early = [];
behav_late = [];

for anum = 1:na
    animal_name = animals{anum};
%     proc_dir            = sprintf('%s%s\\processed_files\', experiment_folder, animal_name);
    cd(experiment_folder)
    an_file = dir([animal_name '*.dat']);
    
    name = [];
    trial = [];
    room_name = [];
    arena_name = [];
    fdate = [];
    ind = 0;
    for i = 1:length(an_file)
            unds = strfind(an_file(i).name, '_');
        frametemp = an_file(i).name(unds(2)+1:end-4);
        if strcmp(frametemp,'Room')
            ind = ind+1;
            name{ind,1} = animal_name;
            trial{ind,1} = an_file(i).name(unds(1)+1:unds(2)-1);
            room_name{ind,1} = an_file(i).name;%(unds(2)+1:end-4);
            arenatemp = [an_file(i).name(1:unds(2)) 'Arena.dat'];
            fdate{ind,1} = an_file(i).date;
            if isfile(arenatemp)
            arena_name{ind,1} = arenatemp;
            else
                warning('no arena match')
            end
        end
    end
    %%
    nt = ind;
    tr_count = 0;
    for trloop = 1:nt
        %%
        is_TR  = contains(trial{trloop}, 'TR') && ~contains(trial{trloop}, 'WTR');
        is_HAB = contains(trial{trloop}, 'HAB');
        if is_TR % (is_TR || is_HAB)
            tr_count = tr_count+1;
            [room, arena, params] = sub_behavior_DAT_tracking_eval(room_name{trloop}, arena_name{trloop}, params);
            sessMinutes( tr_count, anum ) = (room.timestamps(end) - room.timestamps(1))/(60000);
            relative_daytime( tr_count, anum ) = days(room.DAT_fileinfo.datetime - time_origin);
            time_ints = abs(diff([0; room.timestamps(room.entrance_start_idx); room.timestamps(end)]))./1000;
            maxtimeAvoid( tr_count, anum ) = max(time_ints);
            sessName{ tr_count, anum } = trial{trloop};
            fistEntr( tr_count, anum ) = room.first_entrance;
            numShocks( tr_count, anum ) = room.num_shocks;
            numEntr( tr_count, anum ) = room.num_entrances;
            hr = histcounts(room.pol_theta,  [-pi:pi/16:pi], 'Normalization', 'probability');
            ha = histcounts(arena.pol_theta, [-pi:pi/16:pi], 'Normalization', 'probability');
            hr = sort(hr);
            ha = sort(ha);
            null_h = ones(size(hr))./length(hr);
            occ_diff_room( tr_count, anum ) = sum(abs(hr-null_h));
            occ_diff_arena( tr_count, anum ) = sum(abs(ha-null_h));
            
            if strcmp(animal_name, 'HPCACC24500') && tr_count==1
                behav_early.animal = 'HPCACC24500';
                behav_early.trial = trial{trloop};
                behav_early.room = room;
                behav_early.arena = arena;
                behav_early.params = params;
            elseif strcmp(animal_name, 'HPCACC24500') && tr_count==7
                behav_late.animal = 'HPCACC24500';
                behav_late.trial = trial{trloop};
                behav_late.room = room;
                behav_late.arena = arena;
                behav_late.params = params;
            end
        end
    end
    relative_daytime(:,anum) = relative_daytime(:,anum) - nanmin(relative_daytime(:,anum));
    [val, ord] = sort(relative_daytime(:,anum), 'ascend');
    ord = ord(1:tr_count);
    sessMinutes(1:tr_count, anum )      = sessMinutes( ord, anum );
    sessName(1:tr_count, anum )         = sessName( ord, anum );
    relative_daytime(1:tr_count, anum ) = relative_daytime( ord, anum );
    maxtimeAvoid(1:tr_count, anum )     = maxtimeAvoid( ord, anum );
    fistEntr(1:tr_count, anum )         = fistEntr( ord, anum );
    numShocks(1:tr_count, anum )        = numShocks( ord, anum );
    numEntr(1:tr_count, anum )          = numEntr( ord, anum );
    occ_diff_room(1:tr_count, anum )    = occ_diff_room( ord, anum );
    occ_diff_arena(1:tr_count, anum )   = occ_diff_arena( ord, anum );
end
%%
first_tr = [1 1 1 1 1 1 1];
first_good_sess = [16 6 5 5 7 16 12];
vars2save = {'animals', 'animaltype', 'sessMinutes', 'sessName', 'relative_daytime',...
    'behav_early', 'behav_late', 'first_good_sess', 'fistEntr', 'numShocks',...
    'maxtimeAvoid', 'numEntr', 'occ_diff_room', 'occ_diff_arena'};
fname = 'D:\Adrianne\behav_data.mat';
% save(fname, vars2save{:})
vars2eval = vars2save(find(strcmp(vars2save,'fistEntr')):end);
ns = length(vars2eval);
figure(10); clf
dt = [];
anum=7
for j = 1:anum
    dt(j,1) = eval(sprintf('relative_daytime(first_good_sess(j),j);'));
end
for i = 1:ns
    a=NaN(anum,1); b=a;
    for j = 1:anum
        a(j,1) = eval(sprintf('%s(first_tr(j),j);', vars2eval{i}));
        b(j,1) = eval(sprintf('%s(first_good_sess(j),j)', vars2eval{i}));
    end
    subplot(2,3,i); hold on
    bar(mean([a b], 1))
    plot([a b]')
    title(vars2eval{i}, 'Interpreter', 'none')
end

% figure(1023); clf
% set(gcf, 'Position', [128 646 1045 297])
% subplot(1,3,1);
% plot(numEntr./sessMinutes)
% set(gca, 'XTick', [1:nt], 'XTickLabel', trials, 'XTickLabelRotation', -90)
% subplot(1,3,2);
% plot(occ_diff_room)
% set(gca, 'XTick', [1:nt], 'XTickLabel', trials, 'XTickLabelRotation', -90)
% subplot(1,3,3);
% plot(occ_diff_arena)
% legend(animals)   
% set(gca, 'XTick', [1:nt], 'XTickLabel', trials, 'XTickLabelRotation', -90)
% save('D:\Sample Data\HPCACC_behav.mat', 'animals', 'trials', 'numShocks', 'numEntr', 'params')

function [room, arena, params] = sub_behavior_DAT_tracking_eval(room_tracking_fname, arena_tracking_fname, params)
% global params
[room , room_params] = read_APA_csv(room_tracking_fname, []);
[arena, arena_params] = read_APA_csv(arena_tracking_fname, []);
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
    
    behav.entrance_start_idx = find(behav.State(1:end-1)==1 & behav.State(2:end)==2); % 1 == Entrance latency
    behav.shock_start_idx = find(diff([behav.State(1); behav.State] == 2)==1); % 2 == Shock on
    
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

end