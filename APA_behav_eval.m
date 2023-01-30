function [entrances, numShocks, path_length] = APA_behav_eval(experiment_folder, animals, params)
% DAT_Dir = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_water/DAT_files/';
% animals = {'Hipp17861', 'Hipp18240'};
% maxsess = 8; % NaN = habituation sess HAB
% wtr_sess = 8; % Which session number is Water manip (WTR#)
% Plot_habituation = false;
if isempty(animals)
    animals = {'Hipp17861', 'Hipp18240'};
end
nsub = length(animals);
% if Plot_habituation == true
%     sessmatch = [-1, 1:2:maxsess-1; 0, 2:2:maxsess]'; % index of session on the same day
%     allsess = [NaN, 0:maxsess]; % 'NaN' to inluclude habituation session
%     hab_diff = 2;
% else
%     sessmatch = [0, 1:2:maxsess-1; 0, 2:2:maxsess]'; % index of session on the same day
% %     sessmatch = [0 0; 1 2; 3 4; 5 6]; % index of session on the same day
%     allsess = [0:maxsess];
%     hab_diff = 1;
% end

AnimalDir = cell(nsub,1);
dir_list_fname = '_directory_list.csv';
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

params.plotting                 = false;
%%
% sumstats = table('Name', 'FirstEnt', 'NumEnt', 'NumShk');
entrances = NaN(nsub, 1);
numShocks = NaN(nsub, 1);
% first_entrance = NaN(nsub, length(allsess));
path_length = NaN(nsub, 1);
max_sess = 0;
isHAB = cell(nsub, 1);
isWTR = cell(nsub, 1);
isDRK = cell(nsub, 1);
sessLabel_indiv = cell(nsub, 1);
for i = 1:nsub
    %%
    dir_file            = sprintf('%s%s/%s%s', experiment_folder, animals{i}, animals{i}, dir_list_fname);
    try
    AnimalDir{i} = setup_imaging_Sessionfiles(animals{i}, dir_file, experiment_folder);%DAT_Dir, processedDir, contourDir);
        skip_date_create = false;
    catch 
        [AnimalDir{i}] = getDATfiles([experiment_folder 'DAT_files/'], animals{i});
        % use a manual template from the yoked animal if possible
        load(['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\' animals{i} '_matchingsess.mat'])
            skip_date_create = true;
    end
    if ~skip_date_create
        matching_day = zeros(AnimalDir{i}.numSess, 1);
        matching_time = zeros(AnimalDir{i}.numSess, 1);
        matching_sess = zeros(AnimalDir{i}.numSess, 1);
        matching_trial = cell(AnimalDir{i}.numSess, 1);
        isHAB{i} = false(AnimalDir{i}.numSess, 1);
        isWTR{i} = false(AnimalDir{i}.numSess, 1);
        isDRK{i} = false(AnimalDir{i}.numSess, 1);
    dayInd = 0;
    timeInd = 1;
    end
    for j = 1:AnimalDir{i}.numSess
        if ~skip_date_create
            [msInfo] = read_JSON_ms_metaData([AnimalDir{i}.SessionDir{j} 'metaData.json']);
            dateVal = msInfo.recordingStartTime.year*10000 + msInfo.recordingStartTime.month*100 + msInfo.recordingStartTime.day;
            timeVal = msInfo.recordingStartTime.hour*100 + msInfo.recordingStartTime.second;
            if any(matching_day == dateVal)
                timeInd = timeInd+1;
                matching_time(dayInd, timeInd) = timeVal;
            else
                dayInd = dayInd+1;
                timeInd = 1;
            end
            matching_day(dayInd,1) = dateVal;
            matching_time(dayInd, timeInd) = timeVal;
            matching_sess(dayInd, timeInd) = j;
            matching_trial{dayInd, timeInd} = AnimalDir{i}.SessType{j};
            if contains(AnimalDir{i}.SessType{j}, 'HAB')
                isHAB{i}(dayInd, timeInd) = true;
            end
            if contains(AnimalDir{i}.SessType{j}, 'WTR')
                isWTR{i}(dayInd, timeInd) = true;
            end
            if contains(AnimalDir{i}.SessType{j}, 'DRK')
                isDRK{i}(dayInd, timeInd) = true;
            end
        end
        room_dat_fname = AnimalDir{i}.tracking_room{j};
        arena_dat_fname = AnimalDir{i}.tracking_arena{j};
        [room, arena, ~] = behavior_DAT_tracking_eval(room_dat_fname, arena_dat_fname, params);
        pathlength = sum(arena.speed.*arena.dt)./100; % path length in meters
        path_length(i,j) = pathlength; % path length in meters
        entrances(i,j) = room.num_entrances + arena.num_entrances;
        numShocks(i,j) = room.num_shocks + arena.num_shocks;
        sessLabel_indiv{i,j} = AnimalDir{i}.SessType{j};
    end
    if ~skip_date_create
        matching_day   = matching_day(1:dayInd, :);
        matching_time  = matching_time(1:dayInd, :);
        matching_sess  = matching_sess(1:dayInd, :);
        matching_trial = matching_trial(1:dayInd, :);
        isHAB{i} = isHAB{i}(1:dayInd, :);
        isWTR{i} = isWTR{i}(1:dayInd, :);
        isDRK{i} = isDRK{i}(1:dayInd, :);
        sessDays{i}   = matching_day;
        sessTimes{i}  = matching_time;
        sessInd{i}    = matching_sess;
        sessTrials{i} = matching_trial;
    end
    max_sess = max(max_sess, max(matching_sess(:)));
    clearvars room
end

%%
entrances_pm = entrances./path_length;
shocks_pm = numShocks./path_length;
rng(142)
wtr_width = .4;
y_max_scale = 1.5;
allsess_plotting = 1:max_sess;
% if ~Plot_habituation
%     for i = 1:nsub
%        ind = isHAB{i}==1;
%        sessInd{i}(ind) = NaN;
%     end
% end
wtr_sess = [];
for i = 1:nsub
    ind = isWTR{i}==1;
    wtr_sess = cat(2, wtr_sess, sessInd{i}(ind));
end
wtr_sess = unique(wtr_sess);
% end
drk_sess = [];
for i = 1:nsub
    ind = isDRK{i}==1;
    drk_sess = cat(2, drk_sess, sessInd{i}(ind));
end
drk_sess = unique(drk_sess);

nsess = length(allsess_plotting);
acolor = scramble_mat(jet(nsub*3)/1.2, 1);
sessLabel = sessLabel_indiv(1,:);
% sessLabel = cell(nsess,1);
% for i = 1:nsess
%     isHAB{i}
%     if allsess_plotting(i) == 1
%         sessLabel{i} = 'HAB';
%     else
%         sessLabel{i} = sprintf('TR%d', allsess_plotting(i)-2);
%    end
% end
figure(44); clf; 
% vars = {'path_length', 'entrances_pm', 'shocks_pm'};
vars = {'entrances', 'numShocks', 'path_length'};
for vLoop = 1:length(vars)
    eval(sprintf('var%d = %s;', vLoop, vars{vLoop}));
end
subplot(131);  hold on
plot_rects(allsess_plotting, wtr_sess, max(var1(:))*y_max_scale, wtr_width, [.85 .85 1], 'Water')
plot_rects(allsess_plotting, drk_sess, max(var1(:))*y_max_scale, wtr_width, [.25 .25 .25], 'Dark')
% plot_wtr(allsess_plotting, wtr_sess, max(var1(:))*y_max_scale, wtr_width)
norm_axis(allsess_plotting, sessLabel)
ylim([0 max(var1(:))*y_max_scale])
title('Total Entrances')
ind = [-5, -5];
plot(ind, [Inf, Inf], 'k')
plot(ind, [Inf, Inf], 'k:')
lgd_line = legend({'Same Day', 'Separate'});


%
subplot(132);  hold on
title('Total Shocks')
plot_rects(allsess_plotting, wtr_sess, max(var2(:))*y_max_scale, wtr_width, [.85 .85 1], 'Water')
plot_rects(allsess_plotting, drk_sess, max(var2(:))*y_max_scale, wtr_width, [.25 .25 .25], 'Dark')
% plot_wtr(allsess_plotting, wtr_sess, max(var2(:))*y_max_scale, wtr_width)
norm_axis(allsess_plotting, sessLabel)
ylim([0 max(var2(:))*y_max_scale])

subplot(133);  hold on
title('Path Length (m)')
for i = 1:nsub
        ind = [-5, -5];
        plot(ind, [Inf, Inf], '-', 'Color', acolor(i,:))
end
lgd = legend(animals);
plot_rects(allsess_plotting, wtr_sess, max(var3(:))*y_max_scale, wtr_width, [.85 .85 1], 'Water')
plot_rects(allsess_plotting, drk_sess, max(var3(:))*y_max_scale, wtr_width, [.25 .25 .25], 'Dark')
% plot_wtr(allsess_plotting, wtr_sess, max(var3(:))*y_max_scale, wtr_width)
norm_axis(allsess_plotting, sessLabel)
ylim([0 y_max_scale*max(var3(:))])
%
for i = 1:nsub
    for j = 1:size(sessDays{i})
        subplot(131);
        ind = sessInd{i}(j,:); ind = ind(ind>0);
        plot(ind, var1(i, ind), '-', 'Color', acolor(i,:))
        subplot(132);
        plot(ind, var2(i, ind), '-', 'Color', acolor(i,:))
        subplot(133);
        plot(ind, var3(i, ind), '-', 'Color', acolor(i,:))
    end
    subplot(131);
    plot(allsess_plotting, var1(i,:), ':', 'Color', acolor(i,:))
    scatter(allsess_plotting, var1(i,:), 10, 'o', 'MarkerFaceColor', acolor(i,:), 'MarkerEdgeColor', 'k');
    subplot(132);
    plot(allsess_plotting, var2(i,:), ':', 'Color', acolor(i,:))
    scatter(allsess_plotting, var2(i,:), 10, 'o', 'MarkerFaceColor', acolor(i,:), 'MarkerEdgeColor', 'k');
    subplot(133);
    plot(allsess_plotting, var3(i,:), ':', 'Color', acolor(i,:))
    scatter(allsess_plotting, var3(i,:), 10, 'o', 'MarkerFaceColor', acolor(i,:), 'MarkerEdgeColor', 'k');
end
lgd_line.String = lgd_line.String(1:2);
lgd.String = lgd.String(1:nsub);
end
%%
function plot_rects(allsess, wtr_sess, height, spc, color_vec, val_name)
w = find(ismember(allsess, wtr_sess));
if ~isempty(w)
rectangle('Position', [.5, height*.85, 3, height*.1], 'EdgeColor', 'none', 'FaceColor', color_vec)
text(1, height*.9, val_name, 'Color', color_vec/3);
for i = 1:length(w)
    rectangle('Position', [w(i)-spc, 0, 2*spc, height], 'EdgeColor', 'none', 'FaceColor', color_vec)
end
end
end

function norm_axis(allsess, sessLabel)
nsess = length(allsess);
set(gca, 'XTickLabel', sessLabel, 'XTick', allsess, 'XTickLabelRotation', 90)
axis square
xlim([min(allsess)-1 max(allsess)+1])
end
function [file_info] = getDATfiles(datdir, animal)
% datdir  ='C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\DAT_files\';
% animal = 'Hipp17861';
ddir = dir([datdir animal '_*']);
% file_part = cell(size(ddir,1)/2, 1);
%%

% SessionDir = cell(size(ddir,1)/2, 1);
file_info.SessType = cell(size(ddir,1)/2, 1);
file_info.SessInd = NaN(size(ddir,1)/2, 1);
file_info.SessNum = NaN(size(ddir,1)/2, 1);
file_info.tracking_arena = cell(size(ddir,1)/2, 1);
file_info.tracking_room = cell(size(ddir,1)/2, 1);
% AnimalDir{i}.numSess
% AnimalDir{i}.SessType{j}
% AnimalDir{i}.tracking_room{j}
% AnimalDir{i}.tracking_arena{j}
for i = 1:size(ddir,1)
    fn = [ddir(i).folder '\' ddir(i).name];
    if contains(fn, '_Room.dat')% || contains(fn, 'rena.dat')
        x = strfind(fn, animal) + length(animal)+1;
        temp = fn(x:end);
        x2 = strfind(temp, '_')-1;
        sess = temp(1:x2);
        if contains(sess,'HAB')
            ind = 1;
            num = NaN;
        elseif contains(sess,'DRK')
            x1 = strfind(temp, 'DRK')+3;
            num = str2double(temp(x1:x2));
            ind = num+2;
        elseif contains(sess,'TR')
            x1 = strfind(temp, 'TR')+2;
            num = str2double(temp(x1:x2));
            ind = num+2;
        else
            error('unkown session type!')
        end
        file_info.tracking_arena{ind} = [ddir(i).folder '\' animal '_' sess '_Arena.dat'];
        file_info.tracking_room{ind} = [ddir(i).folder '\' animal '_' sess '_Room.dat'];
        file_info.SessType{ind} = sess;
        file_info.SessInd(ind) = ind;
        file_info.SessNum(ind) = num;
    end
end
file_info.numSess = max(file_info.SessInd);

end
