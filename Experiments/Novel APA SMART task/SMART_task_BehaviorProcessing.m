%% ideas to implement
% look at all turns-1:turn, correct vs incorrect
% look at all arena choice changes (turn-1 was L, turn0 is R), correct vs incorrect
% deliberation time at choice point as a function of goal angular distance
% angular distance of choice with percent correct
% half vs ful moon on correct and incorrect trials - see HPCACC34990_RAR21
%
% ideas for later
% need a different model for dark sessions
% increase entry time required for reward
%
% great sessions
% 24524_IL12
%
% problems-
% weird spiral positions in arena position due to inference of position
% from room and arena offset

%%
clear
close all

batch = 3
% cd('E:\RecordingData\GarrettBlair\DREADDs')
switch batch %if true
    case 1
        ddir = "C:\Users\gjb326\Desktop\TRACKER DOCS\\DAT Files\";
        outdir = "C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\";
        anames = {'NAPA_24510', 'NAPA_24511', 'NAPA_24512', 'NAPA_24513', 'NAPA_24522', 'NAPA_24523', 'NAPA_24524', 'NAPA_24525'};
        sleap_dir = 'C:\Users\gjb326\Desktop\SLEAP\predicted outputs 662\';
        sleap_modelname = 'SmartyRats 662 predict';
        sleap_vid_timestamp_dir = 'C:\Users\gjb326\Desktop\SLEAP\SmartyRats 662 input videos\';
        sesstypes = {'IL' 'ILCON'};
        max_tr = 30;
        day0 = datetime(2024,12,19);
        cam_info_filename = 'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\misc\cam_info.mat';
    case 2
        ddir = "C:\Users\gjb326\Desktop\TRACKER DOCS\\DAT Files\";
        outdir = "C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\";
        anames = {'HPCACC34990'};
        sleap_dir = 'C:\Users\gjb326\Desktop\SLEAP\SMART task imaging\predicted_outputs\';
        sleap_modelname = 'HPCACC34990_844-labels__';
        sleap_vid_timestamp_dir = 'C:\Users\gjb326\Desktop\SLEAP\SMART task imaging\';
        sesstypes = {'HAB', 'IL', 'RAR', 'CON'};
        max_tr = 30;
        day0 = datetime(2025,4,4);
        cam_info_filename = 'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\misc\cam_info.mat';
    case 3 % lesion animals
        ddir = "E:\SuaKim\SMART task data\DAT_Files\";
        outdir = "E:\SuaKim\SMART task data\behavior output\";
        anames = {  'SUA34991', 'SUA34992', 'SUA34993', 'SUA34994', 'SUA34995', 'SUA34996',...
            'SUA35008', 'SUA35009', 'SUA35010', 'SUA35011', 'SUA35012', 'SUA35013'};
        % note: 35013 passed away after surgery, no data after session 16
        sleap_dir = 'C:\Users\gjb326\Desktop\SLEAP\SMART task imaging\predicted_outputs\';
        sleap_modelname = '';
        sleap_vid_timestamp_dir = '';
        sesstypes = {'HAB', 'IL', 'RET', 'CON'};
        max_tr = 26;
        day0 = datetime(2025,6,25);
        cam_info_filename = 'E:\SuaKim\SMART task data\behavior output\misc\cam_info.mat';
end
% dreadd_region = {'NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN'};
% anames = {'HPCACC34990'};
dreadd_region = {'NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN'};
% sleap_dir = 'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\sleap_data\';
% ideas
% performance during CCW vs CW2 in saline/c21
% behavior stereotypy diff in exploration under drugs

numan = length(anames);

MAX_TIME_VAL = NaN;
% inj = ["NONE" "SAL" "C21" "CNO" "DCZ"];

params = [];
APA_rat_imaging_params_current;

params.strict_errors_only  = false; % evaluate full circle errors or not

sessdate    = string(zeros(numan, max_tr +1));
filenames   = string(zeros(numan, max_tr +1));
sesstime    = string(zeros(numan, max_tr +1));
anname      = string(zeros(numan, max_tr +1));
region      = string(zeros(numan, max_tr +1));
sessname    = string(zeros(numan, max_tr +1));
sessname_full    = string(zeros(numan, max_tr +1));
sesstype    = string(zeros(numan, max_tr +1));

% sessnum = ones(4,1)*[0:max_tr];

trainNum        = NaN(numan, max_tr +1);
trainNum_sessType= NaN(numan, max_tr +1);
dayNum          = NaN(numan, max_tr +1);
entrpermin      = NaN(numan, max_tr +1);
first_entr      = NaN(numan, max_tr +1);
meantime2entry  = NaN(numan, max_tr +1);
distRan         = NaN(numan, max_tr +1);
dirpref         = NaN(numan, max_tr +1);
sessmins        = NaN(numan, max_tr +1);

reinforced_zone_preference  = NaN(numan, max_tr +1);
right_preference_room       = NaN(numan, max_tr +1);
right_preference_arena      = NaN(numan, max_tr +1);
mean_retrieval_time         = NaN(numan, max_tr +1);
mean_retrieval_time_fake    = NaN(numan, max_tr +1);
swap_performance_all        = NaN(numan, max_tr +1, 2, 2);
swap_performance_win_win    = NaN(numan, max_tr +1);
swap_performance_lose_win   = NaN(numan, max_tr +1);
average_streak              = NaN(numan, max_tr +1);
lapsPerMin                  = NaN(numan, max_tr +1);
probability_correct         = NaN(numan, max_tr +1);
probability_left            = NaN(numan, max_tr +1);
arena_displacement_max      = NaN(numan, max_tr +1);
zone_sampling_angdiff       = NaN(numan, max_tr +1);
zone_sampling_selectivity   = NaN(numan, max_tr +1);


distbins = linspace(0, pi, 19);
correct_by_dist      = NaN(numan, max_tr +1, length(distbins)-1);

cm = lbmap(36,'RedBlue'); 
rats_cmap = [cm(1:2:11,:);cm(end:-2:end-11,:)]; 
c = [.7 .7 .7]; e = [1 .4 0];
rats_cmap = cat(1, c, c, c, c, c, e, c, e, e, e, c, e);

is_expt = [0 0 0 0 0 1 0 1 1 1 0 1]==1;
is_fem  = [1 1 1 1 1 1 0 0 0 0 0 0]==1;
mindist_scale = 1; % mult the size of the reinforcement zone
clrs = [.7 .7 .7; 1, .3, .3];
output_file_name = sprintf('%s%s', outdir, '\SUA_behavior_metrics.mat');
vars2save = whos();
% save('test.mat', vars2save.name)
% SLEAP_resave_tracker_avis;
%
varnames_cam = {'unique_sname', 'pos_center_X', 'pos_center_Y', 'arena_size_px', 'sleap_pos_center_X', 'sleap_pos_center_Y', 'sleap_arena_size_px'};
% cam_info = table(string(unique_sname), ...
%     params_sub.camera_positions.pos_center_X,...
%     params_sub.camera_positions.pos_center_Y,...
%     params_sub.camera_positions.arena_size_px,...
%     params_sub.camera_positions.sleap_pos_center_X,...
%     params_sub.camera_positions.sleap_pos_center_Y,...
%     params_sub.camera_positions.sleap_arena_size_px,...    
%     'VariableNames', varnames(:));
% save(cam_info_filename, 'cam_info');
if isfile(cam_info_filename)
    load(cam_info_filename, 'cam_info');
end
%
%%
% close all
% all_h = NaN(10000, 2);
% all_h2 = NaN(10000, 2);
resave_processing = false;
params.strict_errors_only  = false; % evaluate full circle errors or not
params.during_rotation_only = true;
noSLEAP_ok = true;

%
% close all
% figure(201); clf
figure(202); clf
for a = 1:numan % [6,8,9,10] % [1:5,7,11] % 1:numan
    %%
    fprintf('\n\t____  %s  ____\n', anames{a});
    for s = 0:27
        %%
%         figure(s+100)
        params_sub = params;
        firststr = sprintf('%s%s', ddir, anames{a});
        sname = [];
        for jj = 1:length(sesstypes)
            r_fn = sprintf('%s%s_%s%d_Room.dat', ddir, anames{a}, sesstypes{jj}, s);
            a_fn = sprintf('%s%s_%s%d_Arena.dat', ddir, anames{a}, sesstypes{jj}, s);
            if isfile(r_fn) && isfile(a_fn)
                sname = [sesstypes{jj} num2str(s)];
                stype = [sesstypes{jj}];
                unique_sname = sprintf('%s_%s%d', anames{a}, sesstypes{jj}, s);
                break
            end
        end
        if ~isempty(sname)
            fname = sprintf('%s@SMARTies_%s_%s.mat', outdir, anames{a}, sname);
            sleap_fname = sprintf('%s%s*_%s_%s_TrackerVideo.analysis.h5', sleap_dir, sleap_modelname, anames{a}, sname);
            temp = dir(sleap_fname);
            if ~isempty(temp)
                sleap_fname = sprintf('%s\\%s', temp.folder, temp.name);
            else 
                sleap_fname = [];
            end
            sleap_timestamp_fname = sprintf('%s%s_%s_TrackerVideo_timeStamps.csv', sleap_vid_timestamp_dir, anames{a}, sname);
            sessname(a, s+1)    = sname;
            sesstype(a, s+1)    = stype;
            
            matchIdx = find(strcmp(cam_info.unique_sname, unique_sname));
            if ~isempty(matchIdx)
                save_cam_info=false;
                for ii = 2:length(varnames_cam)
                    eval(sprintf('params_sub.camera_positions.%s = cam_info.%s(matchIdx);', varnames_cam{ii}, varnames_cam{ii}))
                end
            else
                save_cam_info=true;
                disp('no match')
            end

            disp(unique_sname)
            [room, arena, params_sub, sleap] = behavior_SLEAP_tracking_DUAL8MAZE(r_fn, a_fn, sleap_fname, sleap_timestamp_fname, params_sub);
            if ~isempty(sleap) &&(length(room.dat_timestamps) ~= length(sleap.timeStamps))
                figure
                warning(unique_sname)
            end
            if ~isempty(sleap_fname)
                room.pos_file_name  = sleap_fname; 
                arena.pos_file_name = sleap_fname;
            elseif noSLEAP_ok == false
                warning('No sleap file, using .DAT')
%                 fprintf('\t\t %s \n', unique_sname)
%                 [room, arena, params_sub] = behavior_DAT_tracking_eval_DUAL8MAZE(r_fn, a_fn, params_sub);
                room.pos_file_name  = r_fn; 
                arena.pos_file_name = a_fn;
            end
            
            if save_cam_info==true
                n_idx = size(cam_info,1)+1;
                warning('off', 'MATLAB:table:RowsAddedExistingVars')
                cam_info.unique_sname(n_idx) = unique_sname;
                warning('on', 'MATLAB:table:RowsAddedExistingVars')
                for ii = 2:length(varnames_cam)
                    eval(sprintf('cam_info.%s(n_idx) = params_sub.camera_positions.%s;', varnames_cam{ii}, varnames_cam{ii}))
                end   
%                 for ii = 2:length(varnames_cam)
%                     fprintf('cam_info.%s(n_idx) = params_sub.camera_positions.%s;\n', varnames_cam{ii}, varnames_cam{ii})
%                 end 
                save(cam_info_filename, 'cam_info');
            end
            [zonestruct] = DUAL8MAZE_performance_eval(room, arena, params_sub, false);
            set(gcf, 'Name', sprintf('%d - %s', s, unique_sname))
%             figure(10000+a); hold on
            [right_bias, out_pref] = smart_task_right_bias(room, arena, 35, pi/2, false);
            binedge = zonestruct.metrics.sampling_room_arena_bins;%    linspace(-pi, pi, 37);
            bincenters = binedge(2:end) - mean(abs(diff(binedge)));
            reinf_area = abs(rad2deg(bincenters))<=30; 
            opp_area = abs(rad2deg(bincenters))>=150; 
            if false % strcmp(sname, 'RET16') % strcmp(stype, 'RET')
                figure(200);
                if a==1; clf; end
                binedge = zonestruct.metrics.sampling_room_arena_bins;%    linspace(-pi, pi, 37);
                bincenters = binedge(2:end) - mean(abs(diff(binedge)));
                bincenters = [bincenters, bincenters(1)];
                y = zonestruct.metrics.sampling_room; y = [y, y(1)];
%                 subplot(121); 
                polarplot(bincenters, y, 'Color', rats_cmap(a,:))
                hold on; drawnow
%                 polarplot(bincenters, zonestruct.metrics.sampling_arena(2,:), 'Color', rats_cmap(a,:))
%                 plot(zonestruct.metrics.sampling_room(1,:),  'Color', rats_cmap(a,:))
%                 subplot(122); hold on;
%                 plot(zonestruct.metrics.sampling_arena(2,:), 'Color', rats_cmap(a,:))
%                 ylim([-.1 1.1])
            end
            opp_area = abs(rad2deg(bincenters))>=150; 
            if true % strcmp(sname, 'RET16') % strcmp(stype, 'RET')
%                 figure(201); hold on
%                 plot(s + (a-6)/24, zonestruct.metrics.room_resultant_prob(1), 'o', 'Color', clrs(is_expt(a)+1,:))
%                 axis([13 20 -190 190])
                
                figure(202); hold on
                xx = s + (a-6)/100; 
                xx = [xx, xx + zonestruct.metrics.room_resultant_prob(2)];
                yy = zonestruct.metrics.room_resultant_prob(1); 
                yy = abs([yy, yy]);% + zonestruct.metrics.room_resultant_prob(2)];
                plot(xx, yy, 'Color', clrs(is_expt(a)+1,:), 'LineWidth', 2);
                scatter(xx(1), yy(1), 'Marker', 'o', 'MarkerEdgeColor', clrs(is_expt(a)+1,:), 'MarkerFaceColor', 'none');
                scatter(xx(2), yy(2), 'Marker', 'd', 'MarkerEdgeColor', clrs(is_expt(a)+1,:), 'MarkerFaceColor', 'none');
                axis([-1 25 -10 190])
            end
%             set(gcf, 'Name', unique_sname, 'Position', [399   656   841   322])
%             set(gca, 'View', [-15   54])
            t = room.DAT_fileinfo.datetime;
            sessday             = sprintf('%d_%02.0f_%02.0f', t.Year, t.Month, t.Day);
            runtime             = sprintf('%d_%02.0f_%02.0f', t.Hour, t.Minute, t.Second);
            sessdate(a, s+1)    = sessday;
            filenames(a, s+1)   = fname;
            anname(a, s+1)      = anames(a);
            sessname_full(a, s+1)      = unique_sname;
%             region(a, s+1)      = dreadd_region(a);
            sesstime(a, s+1)    = runtime;
            dayNum(a, s+1)      = floor(days( room.DAT_fileinfo.datetime - day0 ));
            trainNum(a, s+1)    = s;
            trainNum_sessType(a, s+1)    = sum(contains( sesstype(a, 1:s), stype));
            distRan(a, s+1)        = sum(arena.speed.*arena.dt)/100; % total meters
            [theta, rho]        = cart2pol(arena.x, arena.y);
            u                   = unwrap(theta);
            dirpref(a, s+1)     = sum(diff(u))/sum(abs(diff(u)));
            sessmins(a, s+1)    = room.timestamps(end)/60000;
            if room.num_entrances > 0
                et = room.timestamps(room.entrance_start_idx);
                m = median(abs(diff([0; et; room.timestamps(end)])))/1000;
                entrpermin(a, s+1) = room.num_entrances/sessmins(a, s+1);
                first_entr(a, s+1) = room.first_entrance;
                meantime2entry(a, s+1) = m;
            elseif arena.num_entrances > 0
                %%%%%
                et = arena.timestamps(arena.entrance_start_idx);
                m = median(abs(diff([0; et; arena.timestamps(end)])))/1000;
                entrpermin(a, s+1) = arena.num_entrances/sessmins(a, s+1);
                first_entr(a, s+1) = arena.first_entrance;
                meantime2entry(a, s+1) = m;
            end
            zone_sampling_angdiff(a, s+1)       = zonestruct.metrics.room_resultant_prob(1);
            zone_sampling_selectivity(a, s+1)   = zonestruct.metrics.room_resultant_prob(2);
            reinforced_zone_preference(a, s+1)  = zonestruct.metrics.reinforced_zone_preference;
            right_preference_room(a, s+1)       = mean(zonestruct.metrics.sampling_room(reinf_area));% - mean(zonestruct.metrics.sampling_room_arena(1,opp_area));
            right_preference_arena(a, s+1)      = mean(zonestruct.metrics.sampling_arena(reinf_area));% - mean(zonestruct.metrics.sampling_room_arena(2,opp_area));
            mean_retrieval_time(a, s+1)         = zonestruct.metrics.mean_retrieval_time;
            mean_retrieval_time_fake(a, s+1)    = zonestruct.metrics.mean_retrieval_time_fake;
            swap_performance_win_win(a, s+1)    = zonestruct.metrics.swap_performance(1,1);
            swap_performance_lose_win(a, s+1)   = zonestruct.metrics.swap_performance(2,1);
            average_streak(a, s+1)              = zonestruct.metrics.average_streak;
            swap_performance_all(a, s+1,:,:)    = zonestruct.metrics.swap_performance;
            arena_displacement_max(a, s+1)      = max(abs(arena.angular_offset));
            hc = binned_statistic1d(zonestruct.choices.choiceDist, zonestruct.choices.isCorrect, distbins, 'mean');
            correct_by_dist(a, s+1, :)          = hc; % histcounts(zonestruct.choices.choiceDist(zonestruct.choices.isCorrect), distbins);
%             probability_correct(a, s+1)         = mean(zonestruct.choices.isCorrect==1);
%             probability_left(a, s+1)            = mean(zonestruct.choices.isLeft==1);
            
            if strcmp(room.DAT_fileinfo.Frame, 'RoomFrame')
                mindist = room.DAT_fileinfo.ReinforcedSector(2); % reinforced sector size, degrees
                mindist = deg2rad(mindist);
            else
                mindist = 0;
            end
%             valid = zonestruct.choices.choiceDist>=0; 
            valid = zonestruct.choices.choiceDist>=mindist*mindist_scale;% & zonestruct.choices.choiceDist<=(pi-mindist);
            %             valid = zonestruct.choices.choiceDist>pi/6; %  & zonestruct.choices.choiceDist>=5*pi/6;
            lapsPerMin(a, s+1)                  = length(zonestruct.choices.isCorrect(valid))/sessmins(a, s+1);
            probability_correct(a, s+1)         = mean(zonestruct.choices.isCorrect(valid)==1);
            probability_left(a, s+1)            = mean(zonestruct.choices.isLeft(valid)==1);
%             if probability_correct(a, s+1)>=.65 && entrpermin(a, s+1)>1
%                 idx1 = find(isnan(all_h(:,1)),1);
%                 idx2 = idx1+length(zonestruct.choices.choiceDist)-1;
%                 all_h(idx1:idx2, :) = [zonestruct.choices.choiceDist, zonestruct.choices.isCorrect];
%             elseif probability_correct(a, s+1)<.65 && entrpermin(a, s+1)>1
%                 idx1 = find(isnan(all_h2(:,1)),1);
%                 idx2 = idx1+length(zonestruct.choices.choiceDist)-1;
%                 all_h2(idx1:idx2, :) = [zonestruct.choices.choiceDist, zonestruct.choices.isCorrect];
%             end
%             h = binned_statistic1d(zonestruct.choices.choiceDist, zonestruct.choices.isCorrect, distbins, 'mean');
            if ~isfile(fname) || resave_processing==true
                vars2save = {'room', 'arena', 'sleap', 'zonestruct', 'params_sub', 'sname', 'stype', 'unique_sname', 'sessday'};
                fprintf('   SAVING %s \n\n', fname)
                save(fname, vars2save{:})
            else
%                 fprintf('    skipped save %s \n', fname)
            end
        else
            %             disp('Uknown sesstype')
            sessname(a, s+1)    = '?';
        end
    end
end
% vars2save = {'trainNum', 'anname', 'sessname_full', 'probability_correct', 'probability_left',...
%     'average_streak', 'swap_performance_lose_win', 'swap_performance_win_win', 'swap_performance_all', 'zone_sampling_angdiff', 'zone_sampling_selectivity', '',...
%     'distRan', 'sessmins', 'entrpermin', 'right_preference_room', 'right_preference_arena', 'correct_by_dist', 'lapsPerMin', '', }; 
% save(sprintf('%s%s', outdir, '\SUA_behavior_metrics.mat'), vars2save{:})
% save(output_file_name, vars2save.name)

%%
SMART_rats_sua_plotting;


%%
figure;
if false
    %%
    temp = load('C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_HPCACC34990_IL8.mat')
    cam_pos = cam_info(find(cam_info.unique_sname == unique_sname), :);
%     [v, vt] = plot_sleap_tracking(temp.sleap, [], 15, 0, temp.sleap.video_path{1}, cam_pos, temp.room);
%     acc = load("D:\GarrettBlair\APA\HPCACC34990\simple_files\2025_04_10_H11_21_44_IL7_@placecells_ACC_miniscope2.mat");
%     hpc = load("D:\GarrettBlair\APA\HPCACC34990\simple_files\2025_04_10_H11_21_44_IL7_@placecells_HPC_miniscope1.mat");
%     acc = load("D:\GarrettBlair\APA\HPCACC34990\processed_files\2025_04_10_H11_21_44_IL7_@placecells_ACC_miniscope2.mat");
%     hpc = load("D:\GarrettBlair\APA\HPCACC34990\processed_files\2025_04_10_H11_21_44_IL7_@placecells_HPC_miniscope1.mat");
%     acc = load("D:\GarrettBlair\APA\HPCACC34990\processed_files\2025_04_07_H11_53_39_HAB2_@placecells_ACC_miniscope2.mat");
%     hpc = load("D:\GarrettBlair\APA\HPCACC34990\processed_files\2025_04_07_H11_53_39_HAB2_@placecells_HPC_miniscope1.mat");
%     acc = load("D:\GarrettBlair\APA\HPCACC34990\processed_files\2025_04_10_H16_32_48_IL8_@placecells_ACC_miniscope2.mat");
%     hpc = load("D:\GarrettBlair\APA\HPCACC34990\processed_files\2025_04_10_H16_32_48_IL8_@placecells_HPC_miniscope1.mat");
    acc = load("D:\GarrettBlair\APA\HPCACC34990\processed_files\2025_04_09_H12_32_43_IL5_@placecells_ACC_miniscope2.mat");
    hpc = load("D:\GarrettBlair\APA\HPCACC34990\processed_files\2025_04_09_H12_32_43_IL5_@placecells_HPC_miniscope1.mat");
    
    %%
    params.skip_ensemble = true;
    [acc.ms.room.momentary_pos_info,  acc_rtemp]  = Fenton_ipos(acc.ms, params.ipos_int_time, 'room', params); disp('acc room')
    [acc.ms.arena.momentary_pos_info,  acc_atemp] = Fenton_ipos(acc.ms, params.ipos_int_time, 'arena', params); disp('acc arena')
    
    [hpc.ms.room.momentary_pos_info, hpc_rtemp]   = Fenton_ipos(hpc.ms, params.ipos_int_time, 'room', params); disp('hpc room')
    [hpc.ms.arena.momentary_pos_info, hpc_atemp]  = Fenton_ipos(hpc.ms, params.ipos_int_time, 'arena', params); disp('hpc arena')
    %%
%     correct_ms = interp1(temp.room.timestamps/1000, single(temp.zonestruct.choices.paths.correct), hpc_rtemp.t, 'nearest')==1;
    correct_ms = interp1(temp.room.timestamps/1000, single(temp.zonestruct.choices.paths.correct), hpc_rtemp.t, 'nearest')==0;
    acc_ipos = abs(acc.ms.room.momentary_pos_info) - abs(acc.ms.arena.momentary_pos_info);
    acc_mean = nanmean(acc_ipos,1); 
%     x = acc_rtemp.x; y = acc_rtemp.y;
    bins = [-50:10:50]; % params.pos_bins
    x = acc_atemp.x; y = acc_atemp.y;
    [ha] = binned_statistic2d(x(correct_ms), y(correct_ms), acc_mean(correct_ms), bins, bins, 'median');
    hpc_ipos = abs(hpc.ms.room.momentary_pos_info) - abs(hpc.ms.arena.momentary_pos_info);    
    hpc_mean = nanmean(hpc_ipos,1); 
%     x = hpc_rtemp.x; y = hpc_rtemp.y;
    x = hpc_atemp.x; y = hpc_atemp.y;
    [hh] = binned_statistic2d(x(correct_ms), y(correct_ms), hpc_mean(correct_ms), bins, bins, 'median');
    
    clrmap = lbmap(100, 'redblue');
    clrmap = [clrmap; [.6 .7 .6]];
    clrmap = clrmap(end:-1:1, :);
    h = [ha hh]; alphamask = ~isnan(h);
    figure; 
    im = imagesc(bins, bins, h, [-.05 .05]);
%     set(im, 'AlphaData', alphamask)
    colormap(clrmap); colorbar
    axis image
    title('median ipos, HPC (left) and ACC (right)');
%     [reduction,umap,clusterIdentifiers,extras]=run_umap(acc_ipos');
    %%
    [accspks, ~] = average_spks_time(acc.spks, 250, acc.time_ms, false, 'mean');
    [hpcspks, ~] = average_spks_time(hpc.spks, 250, hpc.time_ms, false, 'mean');
    [hpccraw, ~] = average_spks_time(hpc.craw, 250, hpc.time_ms, false, 'mean');
    [acccraw, ~] = average_spks_time(acc.craw, 250, acc.time_ms, false, 'mean');
    
    tvec = [0:size(hpcspks,2)-1]*250;
    time2show = 0:250:60*1000;
    time2show = 0:250:vt(end);
    accspks = normalize_rows(2*(accspks>0)+acccraw);
    hpcspks = normalize_rows(2*(hpcspks>0)+hpccraw);
    
    accspks = ones(size(accspks,1),1)*nansum(accspks,1);
    hpcspks = ones(size(hpcspks,1),1)*nansum(hpcspks,1);
    accspks = normalize_rows(accspks);
    hpcspks = normalize_rows(hpcspks);
    
    h_im = cat(3, hpcspks*0, hpcspks, hpcspks);
    a_im = cat(3, accspks, accspks*0, accspks);
%     spks = cat(1, h_im, a_im);
    spks = cat(1, h_im(1:4:end,:,:), a_im(1:4:end,:,:));
    
    correctpath = temp.zonestruct.choices.paths.correct;
    correctpath = conv(correctpath, ones(90,1), 'same')>0;
    time2show = temp.room.timestamps(correctpath);
    time2show = time2show(1:4:end);
    time2show = time2show(time2show<vt(end));
%     time2show = time2show(1:1000);
    %
    twin = [10, 30];
    figure(1); set(gcf, 'Position', [23   489   772   416], 'Color', 'k')
%     time2show = 0:250:160*1000;
    for ii = 1+twin(1):length(time2show)-twin(2)
        %%
        clf
        bidx = find( min(abs(time2show(ii)-vt)) == abs(time2show(ii)-vt), 1 );
        spkidx = find( min(abs(time2show(ii)-tvec)) == abs(time2show(ii)-tvec), 1 );
        
        idx = (spkidx-twin(1)):(spkidx+twin(2));
        idx = idx(idx>1 & idx<=size(spks,2));
        subplot_tight(1,3,[1,2],[0,0])
        image(v{bidx})
        text(370, 25, sprintf('%4.0f sec', vt(bidx)/1000), 'Color', 'w')

        axis image off
        subplot_tight(1,3,3,[0,0])
        image(spks(:,idx,:))
        hold on;
        plot([twin(1), twin(1)]+1, [0 size(spks,1)], 'r--', 'LineWidth', 2)
        text([twin(2)+twin(1)]*.8, 20, sprintf('%4.0f sec', tvec(spkidx)/1000), 'Color', 'w')
        axis image off
        drawnow
        
        
    end
    
end

%%
minc = .5;
otherc = .2;
dt = 300;
for a = 1:numan-1
    figure(1020+a); clf
    set(gcf, 'Position', [100+a*100   442   850   536], 'Name', anames{a})
    hold on
    used = false(length(sessname(a,:)), 1);
    for s = 0:max_tr
        fname = filenames(a,s+1);
        if isfile(fname)
            sess = sesstype(a,s+1);
            sessb = trainNum_sessType(a,s+1)+1;
            
            maxs = max(trainNum_sessType(a, contains(sesstype(a,:), sess)))+1;
            d=load(fname);
            all_previous = true;
            [outvars] = SMART_useful_figs(d.zonestruct, 0, dt, all_previous, false);
            y1 = outvars.correct_runav_time;
            y2 = outvars.left_runav_time;
            x  = outvars.correct_runav_time_xv;
            if sess == 'IL' % yellow
                clr = [linspace(minc, 1, maxs)', ones(maxs,1)*otherc, ones(maxs,1)*otherc];
            elseif sess == 'CON' % cyan
                clr = [ones(maxs,1)*otherc, linspace(minc, 1, maxs)', ones(maxs,1)*otherc];
            elseif sess == 'RAR' % magenta
                clr = [ones(maxs,1)*otherc, ones(maxs,1)*otherc, linspace(minc, 1, maxs)'];
            else
                clr = [];
            end
            if ~isempty(clr)
                used(s+1) = true;
                subplot(1,2,1); hold on
                plot(x,y1, 'Color', clr(sessb, :), 'LineWidth', 2)
                scatter(2000+sessb*10, probability_correct(a, s+1), 'o', 'MarkerFaceAlpha', .7, 'MarkerFaceColor',...
                    clr(sessb, :), 'MarkerEdgeColor', 'k')
                
                subplot(1,2,2); hold on
                plot(x,y2, 'Color', clr(sessb, :), 'LineWidth', 2)
%                 scatter(2000+sessb*50, probability_left(a, s+1), 'o', 'MarkerFaceAlpha', .7, 'MarkerFaceColor',...
%                     clr(sessb, :), 'MarkerEdgeColor', 'k')
                
%                 probability_correct(a, s+1)
%                 probability_left(a, s+1)
            end
        end
    end
end
subplot(1,2,1)
axis([-100 3000 0 1])
title('Percent Correct')
subplot(1,2,2)
title('Percent Left')
axis([-100 3000 0 1])
legend(sessname(1,used))
subplot(1,2,2)

%%
vars = {...
    'Date',...
    'Name',...
    'SessName',...
    'EntrPerMin',...
    'FirstEntr',...
    'MedianEntryTime',...
    'DistanceRan',...
    'DirPreference',...
    'DayNum',...
    'TRNum',...
    'Type'...
    'reinforced_zone_preference',...
    'mean_retrieval_time',...
    'mean_retrieval_time_fake',...
    'arena_displacement_max',...
    'probability_correct',...
    'probability_left',...
    'swap_performance_win_win',...
    'swap_performance_lose_win',...
    'swap_performance_lose_lose',...
    'swap_performance_win_lose',...
    'average_streak'};


B = table(...
    sessdate(:),...
    anname(:),...
    sessname(:),...
    entrpermin(:),...
    first_entr(:),...
    meantime2entry(:),...
    dist(:),...
    dirpref(:),...
    dayNum(:),...
    trainNum(:),...
    sesstype(:),...
    reinforced_zone_preference(:),...
    mean_retrieval_time(:),...
    mean_retrieval_time_fake(:),...
    arena_displacement_max(:),...
    probability_correct(:),...
    probability_left(:),...
    swap_performance_win_win(:),...
    swap_performance_lose_win(:),...
    swap_performance_lose_lose(:),...
    swap_performance_win_lose(:),...
    average_streak(:),...
    'VariableNames', vars);
if sum(contains(B.SessName, '?'))>0
    fprintf('\nRemoving %d rows, could not identify session.\n', sum(contains(B.SessName, '?')))
    B = B(contains(B.SessName, '?')==false, :);
end
% writetable(B, sprintf('%sHPCACC_imaging_behav_summary.csv', outdir))
writetable(B, sprintf('%sbehav_summary.csv', outdir))

%%
vars  = {   'Date',    'Name',   'DREADDs',    'SessName',     'EntrPerMin',  'FirstEntr', 'MaxTimeAvoided', 'DistanceRan', 'DirPreference', 'DayNum', 'TRNum', 'Type'};
xval  = 'TRNum';
yvals = {'MedianEntryTime', 'DistanceRan', 'EntrPerMin', 'average_streak', 'reinforced_zone_preference'};
ylims = {[0, 300], [0 600], [-.1 6.5], [0 10], [0, 1]};
% yvals = {'mean_retrieval_time', 'mean_retrieval_time_fake', 'swap_performance_win_win', 'swap_performance_lose_win'};
% ylims = {[0, 60], [0 60], [-.1 1.1], [-.1 1.1]};
jitter = rand(numan,1)/numan;
cmap = plasma(numan*2);
% cmap = [cmap(1:numan/2,:); cmap(size(cmap,1)-(numan/2)+1:end,:)];
% cmap = cmap(end:-1:1,:);

sessns2show = 1:4;
for i_y = 1:length(yvals)
    figure(i_y); clf; hold on
    for i_sub = 1:numan
        animal_idx = contains(B.Name,anames{i_sub});
        y = eval(sprintf( 'B.%s(animal_idx)', yvals{i_y}));
        if isempty(sessns2show)
            x = eval(sprintf( 'B.%s(animal_idx)', xval));
            x = x+jitter(i_sub);
        else
            x = sessns2show;
            y = y(sessns2show);
        end
        plot(x,y, '-', 'Color', cmap(i_sub,:))
        scatter(x,y, 100, 'o', 'MarkerFaceColor', cmap(i_sub,:),...
            'MarkerFaceAlpha', .7, 'MarkerEdgeColor', 'k')
    end
    if isempty(sessns2show)
        xlim([-2, max_tr+2])
    else
        xlim([-2, max(sessns2show)+2])
    end
    ylim(ylims{i_y})
    title(sprintf('%s', yvals{i_y}), 'Interpreter', 'none')
end



%% Rotation manipulation
% clear 
close all
ddir = 'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\';
if false % NAPA animals behavior only pilot
    expt_str = '@activeplacepref_NAPA_';
    names = {'24510' '24511' '24512' '24513' '24522' '24523' '24524' '24525'};
    sesstype = 'IL';
    sexmark = {'o', 'o', 'o', 'o', '^', '^', '^', '^'};
    nnum = [5 7 10 7 9 8 10 7];
    cmap_sub = cmap([1,1,1,1,1,1,1,1]==1,:);
end
expt_str = '@activeplacepref_HPCACC';
names = {'34990'};
sesstype = 'IL';
sexmark = {'o'};
nnum = [4];
cmap_sub = cmap([1]==1,:);


n = length(names);
%
prev_sess = 4; post_sess = 0;
rew_rate = NaN(n, prev_sess+post_sess+1);
corr_prob = NaN(n, prev_sess++post_sess+1);
for i = 1:n
    sessns = nnum(i)-prev_sess:nnum(i)+post_sess;
    for pn = 1:length(sessns)
        fn = sprintf('%s%s%s_%s%d.mat', ddir, expt_str, names{i}, sesstype, sessns(pn));
        % figure(i);  clf; hold on;
        if isfile(fn)
%             disp(fn)
            temp = load(fn);
            nr = temp.room.num_entrances;
            sesstime = temp.room.timestamps(end) / 60000;
            rew_rate(i, pn) = nr/sesstime;
            nc = temp.zonestruct.metrics.reinforced_zone_preference;
            corr_prob(i, pn) = nc;
        else
            fprintf('\n no file %s', fn)
        end
    end
end
rew_rate;

%%
prob_L = NaN(n, 60);
prob_C = NaN(n, 60);
t_inds = NaN(n, 60);
for i = 1:n
    fn = sprintf('%s%s%s_%s%d.mat', ddir, expt_str, names{i}, sesstype, nnum(i));
    % figure(i);  clf; hold on;
    d=load(fn);
%     nv = length(d.zonestruct.metrics.sequence_prob_Left);
%     prob_L(i, 1:nv) = d.zonestruct.metrics.sequence_prob_Left;
%     prob_C(i, 1:nv) = d.zonestruct.metrics.sequence_prob_Correct;
%     t_inds(i, 1:nv) = d.zonestruct.metrics.sequence_prob_times;
    all_previous = false;
    [outvars] = SMART_useful_figs(d.zonestruct, 0, 300, all_previous, false);
    nv = length(outvars.correct_runav_time);
    prob_C(i, 1:nv) = outvars.correct_runav_time;
    prob_L(i, 1:nv) = outvars.left_runav_time;
    t_inds(i, 1:nv) = outvars.correct_runav_time_xv; 
    

end
prob_C = prob_C(:,~all(isnan(t_inds), 1));
prob_L = prob_L(:,~all(isnan(t_inds), 1));
t_inds = t_inds(:,~all(isnan(t_inds), 1));
%
figure(100); clf;
set(gcf, 'Color', 'w', 'Position', [300   364   223   553], 'Name', 'SMART acquisition')

subplot_tight(2,1,1, [.1 .1]); hold on
title('Task Aquisition - Reward rate')
xs = -prev_sess:0;
clr_rwr = [1 .6 .6; 1 .5 .5];
sz = 50;
% plot([ -10, 10], [1 1], 'k:')

s1= scatter(xs*1000, rew_rate(1, :), sz, 'Marker', sexmark{1}, 'MarkerFaceColor',...
    clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9);
s2= scatter(xs*1000, rew_rate(end, :), sz, 'Marker', sexmark{end}, 'MarkerFaceColor',...
    clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9);

% legend({s1, s2}, {'Females', 'Males'}, 'Location', 'northwest')
meanL = nanmean(rew_rate,1);
stdL = nanstd(rew_rate, [], 1);
shadedErrorBar(xs, meanL, stdL, 'LineProps', {'Color', clr_rwr(2,:)/10, 'LineWidth', 3})

for i = 1:n
    
    mv = 1:length(xs);
    plot(xs, rew_rate(i, mv), '-',  'Color', clr_rwr(2,:), 'LineWidth', 2)
    scatter(xs, rew_rate(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     plot(xs, prob_C(i, mv)*0 + .5, 'k:')
end
ylim([ -.2 4.7])
xlim([ -4.5, 0.5])
% axis square

set(gca, 'XTick', [xs], 'YTick', [0:1:6], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
ylabel('Rewards per minute', 'FontWeight', 'Bold')
xlabel('Session before rotation', 'FontWeight', 'Bold')
legend([s1, s2], {'Females', 'Males'}, 'Location', 'northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot_tight(2,1,2, [.1 .1]); hold on
title('Task Aquisition - Reinforced arm preference')
% plot([ -10, 10], [1 1], 'k:')

s1= scatter(xs*1000, corr_prob(1, :), sz, 'Marker', sexmark{1}, 'MarkerFaceColor',...
    clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9);
s2= scatter(xs*1000, corr_prob(end, :), sz, 'Marker', sexmark{end}, 'MarkerFaceColor',...
    clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9);

% legend({s1, s2}, {'Females', 'Males'}, 'Location', 'northwest')
meanL = nanmean(corr_prob,1);
stdL = nanstd(corr_prob, [], 1);
shadedErrorBar(xs, meanL, stdL, 'LineProps', {'Color', clr_rwr(2,:)/10, 'LineWidth', 3})

for i = 1:n
    
    mv = 1:length(xs);
    plot(xs, corr_prob(i, mv), '-',  'Color', clr_rwr(2,:), 'LineWidth', 2)
    scatter(xs, corr_prob(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     plot(xs, prob_C(i, mv)*0 + .5, 'k:')
end
plot(xs, corr_prob(1, mv)*0 + .5, 'k:')

ylim([ .4 .9])
xlim([ -4.5, 0.5])
% axis square

set(gca, 'XTick', [xs], 'YTick', [0:.1:1], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
ylabel('Reinforced prob', 'FontWeight', 'Bold')
xlabel('Session before rotation', 'FontWeight', 'Bold')
legend([s1, s2], {'Females', 'Males'}, 'Location', 'northwest')

%
figure(101); clf;
set(gcf, 'Color', 'w', 'Position', [300   554   577   363], 'Name', 'SMART acquisition')

subplot_tight(1,2,1, [.1 .1]); hold on
title('Task Aquisition')
xs = -prev_sess:0;
clr_rwr = [1 .6 .6; 1 .5 .5];
sz = 50;
% plot([ -10, 10], [1 1], 'k:')

s1= scatter(xs*1000, rew_rate(1, :), sz, 'Marker', sexmark{1}, 'MarkerFaceColor',...
    clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9);
s2= scatter(xs*1000, rew_rate(end, :), sz, 'Marker', sexmark{end}, 'MarkerFaceColor',...
    clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9);

% legend({s1, s2}, {'Females', 'Males'}, 'Location', 'northwest')
meanL = nanmean(rew_rate,1);
stdL = nanstd(rew_rate, [], 1);
shadedErrorBar(xs, meanL, stdL, 'LineProps', {'Color', clr_rwr(2,:)/10, 'LineWidth', 3})

for i = 1:n
    
    mv = 1:length(xs);
    plot(xs, rew_rate(i, mv), '-',  'Color', clr_rwr(2,:), 'LineWidth', 2)
    scatter(xs, rew_rate(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     plot(xs, prob_C(i, mv)*0 + .5, 'k:')
end
ylim([ -.2 4.7])
xlim([ -5.5, 1.5])
% axis square

text(1, 4.5, sprintf('n=%d', size(rew_rate,1)), 'FontWeight', 'Bold')
set(gca, 'XTick', [xs], 'YTick', [0:1:6], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
ylabel('Rewards per minute', 'FontWeight', 'Bold')
xlabel('Session before rotation', 'FontWeight', 'Bold')
legend([s1, s2], {'Females', 'Males'}, 'Location', 'northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot_tight(1,2,2, [.1 .1]); hold on
title('Behavior shift during rotation')
clr_corr = [.4 1 .8; .3 .8 .6];
clr_left = [.6 .5 1; .4 .3 .8];
sz = 50;
lw = 1;

xmax = 30*60 ;

mv = all(t_inds <= xmax, 1);
meanL = mean(prob_L(:, mv),1);
% stdL = std(prob_L(:, mv), [], 1);%/sqrt(n);
quant = .25;
stdL = abs(meanL-[quantile(prob_L(:, mv), 1-quant, 1); quantile(prob_L(:, mv), quant, 1)]);%/sqrt(n);
meanC = mean(prob_C(:, mv),1);
% stdC = std(prob_C(:, mv), [], 1);%/sqrt(n);
stdC = abs(meanC-[quantile(prob_C(:, mv), 1-quant, 1); quantile(prob_C(:, mv), quant, 1)]);%/sqrt(n);
x = t_inds(1, mv)/60;
offset = 0;
m1 = shadedErrorBar(x-offset, meanL, stdL, 'LineProps', {'Color', clr_left(2,:)/1.5, 'LineWidth', lw*3});
m2 = shadedErrorBar(x+offset, meanC, stdC, 'LineProps', {'Color', clr_corr(2,:)/1.5, 'LineWidth', lw*3});
for i = 1:n
    mv = t_inds(i,:) <= xmax;
    xs = t_inds(i, mv)/60;
    if true
        plot(xs-offset, prob_L(i, mv), '-',  'Color', clr_left(2,:), 'LineWidth', lw)
        plot(xs+offset, prob_C(i, mv), '-', 'Color', clr_corr(2,:), 'LineWidth', lw)
       
    else % fancier
    plot(xs-offset, prob_L(i, mv), '-',  'Color', clr_left(2,:), 'LineWidth', lw)
    plot(xs+offset, prob_C(i, mv), '--', 'Color', clr_corr(2,:), 'LineWidth', lw)
    scatter(xs-offset, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_left(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
    scatter(xs+offset, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_corr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
    end
%     scatter(xs, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
%     scatter(xs, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
end
i=n;
mv = t_inds(i,:) <= xmax;
xs = t_inds(i, mv)/60;
plot(xs, prob_C(i, mv)*0 + .5, 'k:')
% ylim([ -.1 1.1])
% ylim([ .35 .93])
ylim([ -.1 1.2])
xxmax = (xmax + t_inds(i,2)-t_inds(i,1))/60;
xlim([ 0, xxmax])
% xlim([ 0, 30])

rectangle('Position', [10, .45, 15, .01], 'EdgeColor', 'none', 'FaceColor', 'r')
text(13.3, .48,'rotation on', 'Color', 'r')
% axis square
set(gca, 'XTick', [5:5:xmax/60], 'YTick', [.4:.1:1], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
% set(gca, 'XTick', [5:5:25], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off')
ylabel('Arm Choice prob.', 'FontWeight', 'Bold')
xlabel('Session Epoch (minutes)', 'FontWeight', 'Bold')
legend({'Choose Left', 'Choose Correct'}, 'Location', 'northwest', 'LineWidth', .25)
text(xxmax, 0, sprintf('n=%d', size(rew_rate,1)), 'FontWeight', 'Bold')

prob_L(:, mv)
%% CONFLICT SESS
% clear 
close all
ddir = 'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\';
expt_str = '@activeplacepref_NAPA_';
sesstype = 'ILCON';
sexmark = {'o', 'o', 'o', '^', '^', '^', '^'};
names = {'24510' '24511' '24513' '24522' '24523' '24524' '24525'};
nnum = [11 13 11 13 13 14 11];
numan=8;
cmap = plasma(numan*2);
cmap = [cmap(1:numan/2,:); cmap(size(cmap,1)-(numan/2)+1:end,:)];
cmap = cmap(end:-1:1,:);
cmap_sub = cmap([1,1,0,1,1,1,1,1]==1,:);

% dd = {'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL4.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24511_IL5.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24513_IL6.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24525_IL6.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24523_IL7.mat'};
% dd = {'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL5.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24511_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24513_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24525_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24523_IL8.mat'};

n = length(names);
%
prev_sess = 0;
rew_rate = NaN(n, prev_sess+2);
for i = 1:n
    sessns = nnum(i)-prev_sess:nnum(i);
    for pn = 1:length(sessns)
        fn = sprintf('%s%s%s_%s%d.mat', ddir, expt_str, names{i}, sesstype, sessns(pn));
        % figure(i);  clf; hold on;
        if isfile(fn)
%             disp(fn)
            temp = load(fn);
            nr = temp.room.num_entrances;
            sesstime = temp.room.timestamps(end) / 60000;
            rew_rate(i, pn) = nr/sesstime;
        else
            fprintf('\n no file %s', fn)
        end
    end
end
rew_rate;
%
%%
prob_L = NaN(n, 60);
prob_C = NaN(n, 60);
t_inds = NaN(n, 60);
for i = 1:n
    fn = sprintf('%s%s%s_%s%d.mat', ddir, expt_str, names{i}, sesstype, nnum(i));
    % figure(i);  clf; hold on;
    d=load(fn);
%     nv = length(d.zonestruct.metrics.sequence_prob_Left);
%     prob_L(i, 1:nv) = d.zonestruct.metrics.sequence_prob_Left;
%     prob_C(i, 1:nv) = d.zonestruct.metrics.sequence_prob_Correct;
%     t_inds(i, 1:nv) = d.zonestruct.metrics.sequence_prob_times;
    
    d=load(fn);
    all_previous = false;
    ns=5;
    [outvars] = SMART_useful_figs(d.zonestruct, ns, 300, all_previous, false);
    nv = length(outvars.correct_runav_time);
    prob_C(i, 1:nv) = outvars.correct_runav_time;
    prob_L(i, 1:nv) = outvars.left_runav_time;
    t_inds(i, 1:nv) = outvars.correct_runav_time_xv/60;
    
%     nv = ns; % length(outvars.correct_runav_prop_xv);
%     prob_C(i, 1:nv) = outvars.correct_runav_prop;
%     prob_L(i, 1:nv) = outvars.left_runav_prop;
%     t_inds(i, 1:nv) = 1/ns:1/ns:1; % outvars.correct_runav_prop_xv;
end
prob_C = prob_C(:,~all(isnan(t_inds), 1));
prob_L = prob_L(:,~all(isnan(t_inds), 1));
t_inds = t_inds(:,~all(isnan(t_inds), 1));

%
figure(102); clf;
set(gcf, 'Color', 'w', 'Position', [300 364 340   298], 'Name', 'SMART acquisition')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot_tight(1,2,1, [.1 .1]); hold on
title('Behavior shift during conflict, \color[rgb]{.4 .3 .8}LEFT')
clr_corr = [.4 1 .8; .3 .8 .6];
clr_left = [.6 .5 1; .4 .3 .8];
sz = 50;
lw = 1;

xmax = 20;
% xmax = 1;

% mv = all(t_inds <= xmax, 1);
mv = any(t_inds <= xmax, 1);
meanL = nanmean(prob_L(:, mv),1);
% stdL = std(prob_L(:, mv), [], 1);%/sqrt(n);
quant = .25;
stdL = abs(meanL-[quantile(prob_L(:, mv), 1-quant, 1); quantile(prob_L(:, mv), quant, 1)]);%/sqrt(n);
meanC = nanmean(prob_C(:, mv),1);
% stdC = std(prob_C(:, mv), [], 1);%/sqrt(n);
stdC = abs(meanC-[quantile(prob_C(:, mv), 1-quant, 1); quantile(prob_C(:, mv), quant, 1)]);%/sqrt(n);
x = t_inds(1, mv);%/60;
offset = 0;
m1 = shadedErrorBar(x-offset, meanL, stdL, 'LineProps', {'Color', clr_left(2,:)/1.5, 'LineWidth', lw*3});
% m2 = shadedErrorBar(x+offset, meanC, stdC, 'LineProps', {'Color', clr_corr(2,:)/1.5, 'LineWidth', lw*3});
for i = 1:n
    mv = t_inds(i,:) <= xmax;
    xs = t_inds(i, mv);%/60;
    plot(xs-offset, prob_L(i, mv), '-',  'Color', clr_left(2,:), 'LineWidth', lw)
%     plot(xs+offset, prob_C(i, mv), '--', 'Color', clr_corr(2,:), 'LineWidth', lw)
%     scatter(xs-offset, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
%         clr_left(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     scatter(xs+offset, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
%         clr_corr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     scatter(xs, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
%     scatter(xs, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
end
plot(xs, prob_C(i, mv)*0 + .5, 'k:')
ylim([ .11 1.1])
xlim([ 0, xmax + t_inds(i,2)-t_inds(i,1)])%/60)
% xlim([ 0, 30])
% axis square
% set(gca, 'XTick', [0:1/ns:1], 'YTick', [0:.1:1], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
set(gca, 'XTick', [5:5:xmax], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
% set(gca, 'XTick', [5:5:25], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off')
ylabel('Arm Choice prob.', 'FontWeight', 'Bold')
xlabel('Session Epoch (minutes)', 'FontWeight', 'Bold')
% legend({'Choose Left', 'Choose Correct'}, 'Location', 'northwest', 'LineWidth', .25)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot_tight(1,2,2, [.1 .1]); hold on
title('Behavior shift during conflict, \color[rgb]{.3 .8 .6}CORRECT')

% m1 = shadedErrorBar(x-offset, meanL, stdL, 'LineProps', {'Color', clr_left(2,:)/1.5, 'LineWidth', lw*3});
m2 = shadedErrorBar(x+offset, meanC, stdC, 'LineProps', {'Color', clr_corr(2,:)/1.5, 'LineWidth', lw*3});
for i = 1:n
    mv = t_inds(i,:) <= xmax;
    xs = t_inds(i, mv);%/60;
%     plot(xs-offset, prob_L(i, mv), '-',  'Color', clr_left(2,:), 'LineWidth', lw)
    plot(xs+offset, prob_C(i, mv), '-', 'Color', clr_corr(2,:), 'LineWidth', lw)
%     scatter(xs-offset, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
%         clr_left(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     scatter(xs+offset, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
%         clr_corr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     scatter(xs, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
%     scatter(xs, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
end
plot(xs, prob_C(i, mv)*0 + .5, 'k:')
ylim([ .11 1.1])
xlim([ 0, xmax + t_inds(i,2)-t_inds(i,1)])%/60)
% xlim([ 0, 30])
% axis square
% set(gca, 'XTick', [0:1/ns:1], 'YTick', [0:.1:1], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
set(gca, 'XTick', [5:5:xmax], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
% set(gca, 'XTick', [5:5:25], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off')
ylabel('Arm Choice prob.', 'FontWeight', 'Bold')
xlabel('Session Epoch (minutes)', 'FontWeight', 'Bold')
% legend({'Choose Left', 'Choose Correct'}, 'Location', 'northwest', 'LineWidth', .25)

[cdiff, pdiff] = ttest2(prob_C(:,end), prob_C(:,1));
text(t_inds(1,find(mv,1, 'last')), .25, sprintf('e-s p=%1.3f', pdiff))

%% REVERSE DIR SESS
% clear 
close all
ddir = 'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\';
expt_str = '@activeplacepref_NAPA_';
sesstype = 'IL';
sexmark = {'o', 'o', 'o', 'o', '^', '^', '^'};
names = {'24510' '24511' '24512' '24513' '24523' '24524' '24525'};
nnum = [9 11 12 9 10 12 9];

cmap = plasma(numan*2);
cmap = [cmap(1:numan/2,:); cmap(size(cmap,1)-(numan/2)+1:end,:)];
cmap = cmap(end:-1:1,:);
cmap_sub = cmap([1,1,1,1,0,1,1,1]==1,:);

% dd = {'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL4.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24511_IL5.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24513_IL6.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24525_IL6.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24523_IL7.mat'};
% dd = {'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL5.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24511_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24513_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24525_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24523_IL8.mat'};

n = length(names);
%
prev_sess = 0;
rew_rate = NaN(n, prev_sess+2);
for i = 1:n
    sessns = nnum(i)-prev_sess:nnum(i);
    for pn = 1:length(sessns)
        fn = sprintf('%s%s%s_%s%d.mat', ddir, expt_str, names{i}, sesstype, sessns(pn));
        % figure(i);  clf; hold on;
        if isfile(fn)
%             disp(fn)
            temp = load(fn);
            nr = temp.room.num_entrances;
            sesstime = temp.room.timestamps(end) / 60000;
            rew_rate(i, pn) = nr/sesstime;
        else
            fprintf('\n no file %s', fn)
        end
    end
end
rew_rate;
%

prob_L = NaN(n, 10);
prob_C = NaN(n, 10);
t_inds = NaN(n, 10);
for i = 1:n
    fn = sprintf('%s%s%s_%s%d.mat', ddir, expt_str, names{i}, sesstype, nnum(i));
    % figure(i);  clf; hold on;
    if isfile(fn)
        d=load(fn);
        [temp] = SMART_useful_figs(d.zonestruct, 0, 15*60, false, true);
        
        nv = length(temp.correct_runav_time);
        prob_L(i, 1:nv) = temp.left_runav_time;
        prob_C(i, 1:nv) = temp.correct_runav_time;
        t_inds(i, 1:nv) = temp.left_runav_time_xv - (15*60)/2;
    else
        warning('no file %s', fn)
    end
end
%
figure(103); clf;
set(gcf, 'Color', 'w', 'Position', [300 364 1093 553], 'Name', 'SMART acquisition')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot_tight(1,2,1, [.1 .1]); hold on
title('Rotation reversed at 15 min, \color[rgb]{.4 .3 .8}LEFT')
clr_corr = [.4 1 .8; .3 .8 .6];
clr_left = [.6 .5 1; .4 .3 .8];
sz = 50;
lw = 1;

xmax = 30*60 ;

% mv = all(t_inds <= xmax, 1);
mv = any(t_inds <= xmax, 1);
meanL = nanmean(prob_L(:, mv),1);
% stdL = std(prob_L(:, mv), [], 1);%/sqrt(n);
quant = .25;
stdL = abs(meanL-[quantile(prob_L(:, mv), 1-quant, 1); quantile(prob_L(:, mv), quant, 1)]);%/sqrt(n);
meanC = nanmean(prob_C(:, mv),1);
% stdC = std(prob_C(:, mv), [], 1);%/sqrt(n);
stdC = abs(meanC-[quantile(prob_C(:, mv), 1-quant, 1); quantile(prob_C(:, mv), quant, 1)]);%/sqrt(n);
x = nanmedian(t_inds(:, mv), 1)/60;
offset = 0;
m1 = shadedErrorBar(x-offset, meanL, stdL, 'LineProps', {'Color', clr_left(2,:)/1.5, 'LineWidth', lw*3});
% m2 = shadedErrorBar(x+offset, meanC, stdC, 'LineProps', {'Color', clr_corr(2,:)/1.5, 'LineWidth', lw*3});
for i = 1:n
    mv = t_inds(i,:) <= xmax;
    xs = t_inds(i, mv)/60;
    plot(xs-offset, prob_L(i, mv), '-',  'Color', clr_left(2,:), 'LineWidth', lw)
%     plot(xs+offset, prob_C(i, mv), '--', 'Color', clr_corr(2,:), 'LineWidth', lw)
    scatter(xs-offset, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_left(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     scatter(xs+offset, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
%         clr_corr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     scatter(xs, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
%     scatter(xs, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
end
plot(xs, prob_C(i, mv)*0 + .5, 'k:')
plot([15 15], [-1 2], 'k--')
ylim([ -.1 1.1])
xlim([ 0, xmax + t_inds(i,2)-t_inds(i,1)]/60)
% xlim([ 0, 30])
axis square
set(gca, 'XTick', [5:5:xmax/60], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
% set(gca, 'XTick', [5:5:25], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off')
ylabel('Arm Choice prob.', 'FontWeight', 'Bold')
xlabel('Session Epoch (minutes)', 'FontWeight', 'Bold')
% legend({'Choose Left', 'Choose Correct'}, 'Location', 'northwest', 'LineWidth', .25)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot_tight(1,2,2, [.1 .1]); hold on
title('Rotation reversed at 15 min, \color[rgb]{.3 .8 .6}CORRECT')

% m1 = shadedErrorBar(x-offset, meanL, stdL, 'LineProps', {'Color', clr_left(2,:)/1.5, 'LineWidth', lw*3});
m2 = shadedErrorBar(x+offset, meanC, stdC, 'LineProps', {'Color', clr_corr(2,:)/1.5, 'LineWidth', lw*3});
for i = 1:n
    mv = t_inds(i,:) <= xmax;
    xs = t_inds(i, mv)/60;
%     plot(xs-offset, prob_L(i, mv), '-',  'Color', clr_left(2,:), 'LineWidth', lw)
    plot(xs+offset, prob_C(i, mv), '-', 'Color', clr_corr(2,:), 'LineWidth', lw)
%     scatter(xs-offset, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
%         clr_left(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
    scatter(xs+offset, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_corr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     scatter(xs, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
%     scatter(xs, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
end
plot(xs, prob_C(i, mv)*0 + .5, 'k:')
plot([15 15], [-1 2], 'k--')
ylim([ -.1 1.1])
xlim([ 0, xmax + t_inds(i,2)-t_inds(i,1)]/60)
% xlim([ 0, 30])
axis square
set(gca, 'XTick', [5:5:xmax/60], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
% set(gca, 'XTick', [5:5:25], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off')
ylabel('Arm Choice prob.', 'FontWeight', 'Bold')
xlabel('Session Epoch (minutes)', 'FontWeight', 'Bold')
% legend({'Choose Left', 'Choose Correct'}, 'Location', 'northwest', 'LineWidth', .25)
saveas(gcf, 'C:\Users\gjb326\Desktop\TRACKER DOCS\figures\reverse_dir_session.svg')

%% LIGHT DARK LIGHT SESS
% clear 
close all
ddir = 'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\';
expt_str = '@activeplacepref_NAPA_';
sesstype = {'IL', 'ILCON'};
sexmark = {'o', 'o'};
names = {'24510' '24511'};
% nnum = [9 19];
nnum = [9 10];

cmap_sub = cmap([1,1,0,0,0,0,0,0]==1,:);

% dd = {'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL4.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24511_IL5.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24513_IL6.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24525_IL6.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24523_IL7.mat'};
% dd = {'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL5.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24511_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24513_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24525_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24523_IL8.mat'};

n = length(names);
%
prev_sess = 2;
rew_rate = NaN(n, prev_sess+1);
corr_prob = NaN(n, prev_sess+1);
for i = 1:n
    sessns = nnum(i)-prev_sess:nnum(i);
    for pn = 1:length(sessns)
        for stype = 1:length(sesstype)
        fn = sprintf('%s%s%s_%s%d.mat', ddir, expt_str, names{i}, sesstype{stype}, sessns(pn));
        % figure(i);  clf; hold on;
        if isfile(fn)
%             disp(fn)
            temp = load(fn);
            
            nr = temp.room.num_entrances;
            nc = temp.zonestruct.metrics.reinforced_zone_preference;
            sesstime = temp.room.timestamps(end) / 60000;
            rew_rate(i, pn) = nr/sesstime;
            corr_prob(i, pn) = nc;
        else
            fprintf('\n no file %s', fn)
        end
        end
    end
end
%

prob_L = NaN(n, 10);
prob_C = NaN(n, 10);
t_inds = NaN(n, 10);
for i = 1:n
    for stype = 1:length(sesstype)
        fn = sprintf('%s%s%s_%s%d.mat', ddir, expt_str, names{i}, sesstype{stype}, nnum(i));
        % figure(i);  clf; hold on;
        if isfile(fn)
            d=load(fn);
            nv = length(d.zonestruct.metrics.sequence_prob_Left);
            prob_L(i, 1:nv) = d.zonestruct.metrics.sequence_prob_Left;
            prob_C(i, 1:nv) = d.zonestruct.metrics.sequence_prob_Correct;
            t_inds(i, 1:nv) = d.zonestruct.metrics.sequence_prob_times;
        else
            warning('no file %s', fn)
        end
    end
end
%%
figure(104); clf;
set(gcf, 'Color', 'w', 'Position', [300 364 1093 553], 'Name', 'SMART acquisition')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot_tight(1,2,1, [.1 .1]); hold on
title('Task Aquisition')
xs = [-1 0 1];
clr_rwr = [1 .6 .6; 1 .5 .5];
sz = 50;
% plot([ -10, 10], [1 1], 'k:')

s1= scatter(xs*1000, rew_rate(1, :), sz, 'Marker', sexmark{1}, 'MarkerFaceColor',...
    clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9);
s2= scatter(xs*1000, rew_rate(end, :), sz, 'Marker', sexmark{end}, 'MarkerFaceColor',...
    clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9);

% legend({s1, s2}, {'Females', 'Males'}, 'Location', 'northwest')
meanL = nanmean(rew_rate,1);
stdL = nanstd(rew_rate, [], 1);
shadedErrorBar(xs, meanL, stdL, 'LineProps', {'Color', clr_rwr(2,:)/10, 'LineWidth', 3})

for i = 1:n
    
    mv = 1:length(xs);
    plot(xs, rew_rate(i, mv), '-',  'Color', clr_rwr(2,:), 'LineWidth', 2)
    scatter(xs, rew_rate(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     plot(xs, prob_C(i, mv)*0 + .5, 'k:')
end
ylim([ -.2 7.7])
xlim([ -1.5, 1.5])
axis square

set(gca, 'XTick', [xs], 'XTickLabel', {'Light', 'Dark', 'Light'}, 'YTick', [0:1:6], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
ylabel('Rewards per minute', 'FontWeight', 'Bold')
xlabel('Session before rotation', 'FontWeight', 'Bold')
legend([s1, s2], {'Females', 'Males'}, 'Location', 'northwest')


subplot_tight(1,2,2, [.1 .1]); hold on
title('Task Aquisition')
% plot([ -10, 10], [1 1], 'k:')

s1= scatter(xs*1000, corr_prob(1, :), sz, 'Marker', sexmark{1}, 'MarkerFaceColor',...
    clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9);
s2= scatter(xs*1000, corr_prob(end, :), sz, 'Marker', sexmark{end}, 'MarkerFaceColor',...
    clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9);

% legend({s1, s2}, {'Females', 'Males'}, 'Location', 'northwest')
meanL = nanmean(corr_prob,1);
stdL = nanstd(corr_prob, [], 1);
shadedErrorBar(xs, meanL, stdL, 'LineProps', {'Color', clr_rwr(2,:)/10, 'LineWidth', 3})

for i = 1:n
    
    mv = 1:length(xs);
    plot(xs, corr_prob(i, mv), '-',  'Color', clr_rwr(2,:), 'LineWidth', 2)
    scatter(xs, corr_prob(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_rwr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     plot(xs, prob_C(i, mv)*0 + .5, 'k:')
end
plot(xs, corr_prob(1, mv)*0 + .5, 'k:')
ylim([ -.1 1.1])
xlim([ -1.5, 1.5])
axis square

set(gca, 'XTick', [xs], 'XTickLabel', {'Light', 'Dark', 'Light'}, 'YTick', [0:1:6], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
ylabel('Rewards per minute', 'FontWeight', 'Bold')
xlabel('Session before rotation', 'FontWeight', 'Bold')
legend([s1, s2], {'Females', 'Males'}, 'Location', 'northwest')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot_tight(1,2,1, [.1 .1]); hold on
title('Rotation reversed at 15 min, \color[rgb]{.4 .3 .8}LEFT')
clr_corr = [.4 1 .8; .3 .8 .6];
clr_left = [.6 .5 1; .4 .3 .8];
sz = 50;
lw = 1;

xmax = 30*60 ;

% mv = all(t_inds <= xmax, 1);
mv = any(t_inds <= xmax, 1);
meanL = nanmean(prob_L(:, mv),1);
% stdL = std(prob_L(:, mv), [], 1);%/sqrt(n);
quant = .25;
stdL = abs(meanL-[quantile(prob_L(:, mv), 1-quant, 1); quantile(prob_L(:, mv), quant, 1)]);%/sqrt(n);
meanC = nanmean(prob_C(:, mv),1);
% stdC = std(prob_C(:, mv), [], 1);%/sqrt(n);
stdC = abs(meanC-[quantile(prob_C(:, mv), 1-quant, 1); quantile(prob_C(:, mv), quant, 1)]);%/sqrt(n);
x = nanmedian(t_inds(:, mv), 1)/60;
offset = 0;
m1 = shadedErrorBar(x-offset, meanL, stdL, 'LineProps', {'Color', clr_left(2,:)/1.5, 'LineWidth', lw*3});
% m2 = shadedErrorBar(x+offset, meanC, stdC, 'LineProps', {'Color', clr_corr(2,:)/1.5, 'LineWidth', lw*3});
for i = 1:n
    mv = t_inds(i,:) <= xmax;
    xs = t_inds(i, mv)/60;
    plot(xs-offset, prob_L(i, mv), '-',  'Color', clr_left(2,:), 'LineWidth', lw)
%     plot(xs+offset, prob_C(i, mv), '--', 'Color', clr_corr(2,:), 'LineWidth', lw)
    scatter(xs-offset, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_left(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     scatter(xs+offset, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
%         clr_corr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     scatter(xs, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
%     scatter(xs, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
end
plot(xs, prob_C(i, mv)*0 + .5, 'k:')
plot([15.5 15.5], [-1 2], 'k--')
ylim([ -.1 1.1])
xlim([ 0, xmax + t_inds(i,2)-t_inds(i,1)]/60)
% xlim([ 0, 30])
axis square
set(gca, 'XTick', [5:5:xmax/60], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
% set(gca, 'XTick', [5:5:25], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off')
ylabel('Arm Choice prob.', 'FontWeight', 'Bold')
xlabel('Session Epoch (minutes)', 'FontWeight', 'Bold')
% legend({'Choose Left', 'Choose Correct'}, 'Location', 'northwest', 'LineWidth', .25)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot_tight(1,2,2, [.1 .1]); hold on
title('Rotation reversed at 15 min, \color[rgb]{.3 .8 .6}CORRECT')

% m1 = shadedErrorBar(x-offset, meanL, stdL, 'LineProps', {'Color', clr_left(2,:)/1.5, 'LineWidth', lw*3});
m2 = shadedErrorBar(x+offset, meanC, stdC, 'LineProps', {'Color', clr_corr(2,:)/1.5, 'LineWidth', lw*3});
for i = 1:n
    mv = t_inds(i,:) <= xmax;
    xs = t_inds(i, mv)/60;
%     plot(xs-offset, prob_L(i, mv), '-',  'Color', clr_left(2,:), 'LineWidth', lw)
    plot(xs+offset, prob_C(i, mv), '-', 'Color', clr_corr(2,:), 'LineWidth', lw)
%     scatter(xs-offset, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
%         clr_left(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
    scatter(xs+offset, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor',...
        clr_corr(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .9)
%     scatter(xs, prob_L(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
%     scatter(xs, prob_C(i, mv), sz, 'Marker', sexmark{i}, 'MarkerFaceColor', cmap_sub(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
end
plot(xs, prob_C(i, mv)*0 + .5, 'k:')
plot([15.5 15.5], [-1 2], 'k--')
ylim([ -.1 1.1])
xlim([ 0, xmax + t_inds(i,2)-t_inds(i,1)]/60)
% xlim([ 0, 30])
axis square
set(gca, 'XTick', [5:5:xmax/60], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off', 'FontWeight', 'Bold')
% set(gca, 'XTick', [5:5:25], 'YTick', [0:.25:1], 'FontName', 'Arial', 'Box', 'off')
ylabel('Arm Choice prob.', 'FontWeight', 'Bold')
xlabel('Session Epoch (minutes)', 'FontWeight', 'Bold')
% legend({'Choose Left', 'Choose Correct'}, 'Location', 'northwest', 'LineWidth', .25)







%% ROOM ARENA ROOM manips
clear
close all
cd('C:\Users\gjb326\Desktop\TRACKER DOCS\figures')

ddir = 'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\';
expt_str = '@activeplacepref_NAPA_';
sesstype = {'ILCON'};
sexmark = {'o', 'o', 'o', 'o', '^', '^', '^'};
names = {'24513' '24525'};
% nnum = [9 19];
anum =  [11 11 13 13 22 25 25];
aidx =  [2  2  4  4  5   8  8];
nnum =  [20 23 17 21 16 18 21];
a_marker = {'o', 'o', 'd', 'd', 's', '^', '^'};
a_marker = {'o', 'o', 'o', 'o', 'o', 'o', 'o'};

numan=8;
cmap = plasma(numan*2);
cmap = [cmap(1:numan/2,:); cmap(size(cmap,1)-(numan/2)+1:end,:)];
cmap = cmap(end:-1:1,:);
cmap_sub = cmap;

% dd = {'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL4.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24511_IL5.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24513_IL6.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24525_IL6.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24523_IL7.mat'};
% dd = {'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL5.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24511_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24513_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24525_IL7.mat';...
%     'C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24523_IL8.mat'};

n = length(anum);
%
rew_rate = NaN(n, 3);
corr_prob = NaN(n, 3);

nsplit = 3;
prob_L = NaN(n, 3, nsplit);
prob_C = NaN(n, 3, nsplit);

prob_Lm = NaN(n, 3);
prob_Cm = NaN(n, 3);
num = NaN(n, 3);

close all
figure(1005); % clf
set(gcf, 'Name', 'Room->Arena->Room', 'Color', 'w', 'Position', [180   677  212   201])

% ax1= subplot(3,1,1);
% ax3= subplot(3,1,2);
% ax5= subplot(3,1,3);
% linkaxes([ax1,ax3,ax5],'x')
xinds_1 = [2167, 180, 3953];
xinds_2 = xinds_1+5400;
highlight_subj = 5;
spdconv = gausswin(30*2 + 1); spdconv = spdconv./sum(spdconv(:));
for i = 1:n
    %
    sessns = nnum(i):nnum(i)+2;
    for pn = 1:length(sessns)
%         figure
        for stype = 1:length(sesstype)
        fn = sprintf('%s%s245%d_%s%d.mat', ddir, expt_str, anum(i), sesstype{stype}, sessns(pn));
        if isfile(fn)
%             disp(fn)
            temp = load(fn);
            
            all_previous = true;
            [outvars] = SMART_useful_figs(temp.zonestruct, nsplit, 0, all_previous, false);
            prob_C(i, pn, :) = outvars.correct_runav_prop;
            prob_L(i, pn, :) = outvars.left_runav_prop;
            
            prob_Cm(i, pn) = outvars.correct_runav_prop(end);
            prob_Lm(i, pn) = outvars.left_runav_prop(end);
            num(i, pn) = length(temp.zonestruct.choices.isIDX);


            if i == highlight_subj
                y = rad2deg(temp.room.pol_theta); 
                spd = temp.arena.speed';
                spd = conv(spd, spdconv, 'same');
                path_lin = NaN(length(y), 1);
                path_lin_label = 0*y;
                z = temp.zonestruct.choices;
                for jj = 1:size(z.L_start_end)
                    inds = z.L_start_end(jj,1):z.L_start_end(jj,2);
                    s = spd(inds);
                    s = cumsum(s)./sum(s);
                    path_lin(inds) = s.*linspace(0,1,length(inds));
                    path_lin(inds(1)) = -1;
                    path_lin_label(inds) = 1;
                end
                for jj = 1:size(z.R_start_end)
                    inds = z.R_start_end(jj,1):z.R_start_end(jj,2);
                    s = spd(inds);
                    s = cumsum(s)./sum(s);
                    path_lin(inds) = s.*linspace(0,1,length(inds));
                    path_lin(inds(1)) = -1;
                    path_lin_label(inds) = 2;
                end
                plotinds = ~isnan(path_lin);
%                 figure(1004); % clf
                subplot_tight(3,1,(pn-1)*1 + 1, [.02, .1]); hold on;
                cla%                 x = temp.arena.x; y = temp.arena.y; t = temp.arena.timestamps./1000; 
                if pn==2
                    t = temp.arena.timestamps./1000; 
                    m = temp.arena.entrance_start_idx;
                else
                    t = temp.room.timestamps./1000; 
                    m = temp.room.entrance_start_idx;
                end
                
                path_lin_reward = path_lin_label>10;
                path_lin_reward(m) = true;

                path_lin = path_lin(plotinds);
                path_lin(path_lin==-1) = NaN;
                path_lin_label = path_lin_label(plotinds);
                path_lin_reward = path_lin_reward(plotinds);
                path_lin_L = path_lin*NaN; path_lin_L(path_lin_label==1) = path_lin(path_lin_label==1);
                path_lin_R = path_lin*NaN; path_lin_R(path_lin_label==2) = path_lin(path_lin_label==2);
                pl = find(path_lin_reward==true);
                
                plot(path_lin_L, 'b-', 'LineWidth', 2)
                plot(path_lin_R, 'r-', 'LineWidth', 2)
                clrs = {[.7 .7 1], [1 .7 .7]};
                for jjj = 1:length(pl)
                scatter(pl(jjj), path_lin(pl(jjj)), 10, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clrs{path_lin_label(pl(jjj))}, 'MarkerFaceAlpha', 1)
                end
                if pn==3; plot([4500 4500+find(t>=60, 1)], [-.25 -.25], 'k-', 'LineWidth', 2); end
                
                ylim([-.5 1.5])
                xlim([xinds_1(pn), xinds_2(pn)])
%                 xlim([0 12050])
%                 xlim([4233 9565])
%                 xlim([2000 12050])
                set(gca, 'YTick', [0 .5 1], 'XTick', [], 'TickDir', 'out', 'XColor', 'none', 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold')
                
                if pn==1; keep_arena = temp.arena;  keep_room = temp.room; keep_zone = temp.zonestruct; end
            end
            
        else
            fprintf('\n no file %s', fn)
        end
        end
    end
end
%%
if false
figure(10); clf; hold on
inds=19800:20465;
ax = keep_arena.x; ay = keep_arena.y;
axl = keep_arena.x*NaN; ayl = keep_arena.y*NaN;
axl(keep_zone.choices.paths.left) = ax(keep_zone.choices.paths.left); 
ayl(keep_zone.choices.paths.left) = ay(keep_zone.choices.paths.left); 
axr = keep_arena.x*NaN; ayr = keep_arena.y*NaN;
axr(keep_zone.choices.paths.right) = ax(keep_zone.choices.paths.right); 
ayr(keep_zone.choices.paths.right) = ay(keep_zone.choices.paths.right); 
r = intersect(keep_room.entrance_start_idx, inds);
% plot3(ax,ay, 1:length(ay), 'k', 'LineWidth', 2)
% plot3(axl,ayl, 1:length(ay), 'b', 'LineWidth', 3)
% plot3(axr,ayr, 1:length(ay), 'r', 'LineWidth', 3)
plot(ax(inds),ay(inds), 'k', 'LineWidth', 2)
plot(axl(inds),ayl(inds), 'b', 'LineWidth', 3)
plot(axr(inds),ayr(inds), 'r', 'LineWidth', 3)
plot([0 40], [-45 -45], 'k', 'LineWidth', 3)
scatter(ax(r),ay(r), 1000, 'k.')
axis([-50 50 -50 50])
end
% zlim([inds(1), inds(end)])
%%
figure(1004); clf
set(gcf, 'Name', 'Room->Arena->Room', 'Color', 'w', 'Position', [180   377   731   201])
clrs = [1 .4 .4; .4 .5 1; 1 .4 .4];
cmap_sub = cmap.*0;
offs = linspace(-.1, .1, nsplit);
% for ii = 1:3
%     ys = squeeze(prob_C(:,ii,:));
%     my = mean(ys, 1);
%     sy = std(ys, [], 1);
%     xs = (1:nsplit);
%     subplot(1,4,ii); hold on;
%     shadedErrorBar(xs, my, sy, 'LineProps', {'Color', clrs(ii,:), 'LineWidth', 2})
%     plot([0 nsplit], [.5 .5], 'k:', 'LineWidth', .5)
%     for iii = 1:n
%         plot(xs-offs, ys(iii, :), 'LineWidth', .1, 'Color', cmap_sub(aidx(iii), :))
%     end
%     iii = highlight_subj;
%     plot(xs-offs, ys(iii, :), 'LineWidth', .25, 'Color', 'm')
%     
%     ylim([-.1 1.1])
%     xlim([0.5  nsplit+.5])
%     set(gca, 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold')
%     if ii>1
%         set(gca, 'YTickLabel', '')
%     end
% end


subplot(2,4,4); hold on;
cmap_sub = cmap*0;
xs = 1:3;
my = mean(prob_Cm, 1);
offsets = (rand(8,1)-.5)/4;
offsets = offsets*ones(1,3);
sy = std(prob_Cm, [], 1);
   plot([0 4], [.5 .5], 'k:', 'LineWidth', .5)
shadedErrorBar(xs, my, sy, 'LineProps', {'Color', clrs(1,:), 'LineWidth', 2})
for ii = [1:n]%, 3]
    ys = squeeze(prob_Cm(ii, :));
    subplot(2,4,4); hold on;
    plot(xs+offsets(ii,:), ys, 'LineWidth', .1, 'Color', 'k')
%     scatter(xs, ys, 30, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap_sub(aidx(ii), :), 'MarkerFaceAlpha', .7)
    ylim([-.1 1.1])
    xlim([.5  3.5])
    set(gca, 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold')
    if false % ii==highlight_subj
        scatter(xs+offsets(ii,:), ys, 30, 'Marker', a_marker{ii}, 'MarkerEdgeColor', [1 .6 1], 'MarkerFaceColor', 'm', 'MarkerFaceAlpha', 1)
    else
        scatter(xs+offsets(ii,:), ys, 30, 'Marker', a_marker{ii}, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', cmap_sub(aidx(ii), :), 'MarkerFaceAlpha', 1)
    end
end
axis([.75 3.25 .2 .9])
set(gca, 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold', 'TickDir', 'out', 'YTick', [.25:.25:.75])
% subplot(2,4,8); hold on;
%%

D = readtable('C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\behav_summary.csv');

valid = D.swap_performance_win_win>0 & D.average_streak>0 & D.swap_performance_lose_win>0;
ww = 1-D.swap_performance_win_win(valid);
ll = D.average_streak(valid);
[c, p] = corr(1-ww,ll, 'type', 'Spearman');
[lreg, ~] = fit(ll, ww,'exp1','Normalize','on','Robust','Bisquare');

figure(1006); clf; 
set(gcf, 'Color', 'w', 'Position', [200 50 500 450]); 
xs = linspace(min(ll), max(ll), 100);
plot(xs, lreg(xs), 'LineWidth', 3, 'Color', [.4 .8 .8])
hold on;
scatter(ll, ww, 200, 'Marker', '.', 'MarkerEdgeColor', 'k')
% scatter(ll, ww, 50, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .2)
text(8, 0, sprintf('n = %d', sum(valid)))
set(gca, 'YTick', [0 .2 .4 .6 .8 1], 'YDir', 'reverse', 'YTickLabel', 1:-.2:0, 'XTick', [0 4 8 12], 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold')
xlabel('Average Streak')
ylabel('Win-Shift choices')
title(sprintf('Av Streak and Win-shift  -  corr=%1.2f, p=%1.5f', c, p))
axis([0 13 -.2 1.2])


% valid = D.swap_performance_win_win>0 & D.EntrPerMin>0;
ww = 1-D.swap_performance_win_win(valid);
ll = D.EntrPerMin(valid);
[c, p] = corr(1-ww,ll, 'type', 'Spearman');
[lreg, ~] = fit(ll, ww,'poly1','Normalize','on','Robust','Bisquare');

figure(1007); clf; 
set(gcf, 'Color', 'w', 'Position', [200 500 500 450]); 
xs = linspace(min(ll), max(ll), 100);
plot(xs, lreg(xs), 'LineWidth', 3, 'Color', [.4 .8 .8])
hold on;
scatter(ll, ww, 50, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .2)
text(5, 0, sprintf('n = %d', sum(valid)))
set(gca, 'YTick', [0 .2 .4 .6 .8 1], 'YDir', 'reverse', 'YTickLabel', 1:-.2:0, 'XTick', [0 3 6], 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold')
xlabel('Entries / Min')
ylabel('Win-Shift choices')
title(sprintf('Rew/Min and Win-shift  -  corr=%1.2f, p=%1.5f', c, p))
axis([-1 7 -.2 1.2])


% valid = D.swap_performance_lose_win>0 & D.average_streak>0;
ww = D.swap_performance_lose_win(valid);
ll = D.average_streak(valid);
[c, p] = corr(ww,ll, 'type', 'Spearman');
[lreg, ~] = fit(ll, ww,'exp1','Normalize','on','Robust','Bisquare');

figure(1008); clf; 
set(gcf, 'Color', 'w', 'Position', [800 50 500 450]); 
xs = linspace(min(ll), max(ll), 100);
plot(xs, lreg(xs), 'LineWidth', 3, 'Color', [.8 .8 .4])
hold on;
scatter(ll, ww, 50, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .2)
text(8, 0, sprintf('n = %d', sum(valid)))
set(gca, 'YTick', [0 .2 .4 .6 .8 1], 'YDir', 'normal', 'XTick', [0 4 8 12], 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold')
xlabel('Average Streak')
ylabel('Lose-Shift choices')
title(sprintf('Av Streak and Lose-shift  -  corr=%1.2f, p=%1.5f', c, p))
axis([0 13 -.2 1.2])


% valid = D.swap_performance_lose_win>0 & D.EntrPerMin>0;
ww = D.swap_performance_lose_win(valid);
ll = D.EntrPerMin(valid);
[c, p] = corr(ww,ll, 'type', 'Spearman');
[lreg, ~] = fit(ll, ww,'poly1','Normalize','on','Robust','Bisquare');

figure(1009); clf; 
set(gcf, 'Color', 'w', 'Position', [800 500 500 450]); 
xs = linspace(min(ll), max(ll), 100);
plot(xs, lreg(xs), 'LineWidth', 3, 'Color', [.8 .8 .4])
hold on;
scatter(ll, ww, 50, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .2)
text(5, 0, sprintf('n = %d', sum(valid)))
set(gca, 'YTick', [0 .2 .4 .6 .8 1], 'YDir', 'normal', 'XTick', [0 3 6], 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold')
xlabel('Entries / Min')
ylabel('Lose-Shift choices')
title(sprintf('Rew/Min and Lose-shift  -  corr= %1.2f, p= %1.5f', c, p))
axis([-1 7 -.2 1.2])


%%
% valid = D.arena_displacement_max>6*pi;
valid1 = ~isnan(D.probability_correct) & ~isnan(D.probability_left) & (D.probability_correct ~= D.probability_left);
valid2 = D.DistanceRan>80;
valid=valid1&valid2;
% valid = D.swap_performance_win_win>0;
p11 = D.swap_performance_win_win(valid);
p12 = -100+D.swap_performance_win_lose(valid);
p21 = D.swap_performance_lose_win(valid);
p22 = -100+D.swap_performance_lose_lose(valid);

p = NaN(2,2,sum(valid));
p(1,1,:) = D.swap_performance_win_win(valid);
p(1,2,:) = D.swap_performance_win_lose(valid);
p(2,1,:) = D.swap_performance_lose_win(valid);
p(2,2,:) = D.swap_performance_lose_lose(valid);


offs = (rand(sum(valid),4)-.5)/10;
offs = ones(sum(valid),1)*[1.1,1.9,3,4] + offs;
figure(99); clf; 
% set(gcf, 'Name', 'win-lose-perf', 'Color', 'w', 'Position', [283   506   445   243])
set(gcf, 'Name', 'win-lose-perf', 'Color', 'w', 'Position', [283   506   257   231])
vv = violinplot([p11 p21 p12 p22]);
rectangle('Position', [1 0 1 1], 'FaceColor', 'w', 'EdgeColor', 'w')
plot(offs', [p11 p21 NaN*p12 NaN*p22]', 'Color', 'k')
cm = viridis(50);
cmean = nanmean(squeeze([p(1,1,:), p(2,1,:), p(1,2,:), p(2,2,:)]), 2);
for ii = 1:2; 
    vv(ii).ViolinColor = cm(round(100*cmean(ii)), :);
    vv(ii).MedianColor = 'none';
    vv(ii).BoxColor = 'none';
    vv(ii).EdgeColor = 'k';
    vv(ii).ViolinAlpha = .8;
    vv(ii).ShowData = 0;
end
vv = violinplot([p11 p21 p12 p22]);
for ii = 1:2
    vv(ii).ViolinColor = 'none';
    vv(ii).MedianColor = light_colormap(cm(round(100*cmean(ii)), :), 2);
    vv(ii).MedianColor = 'w'; % cm(round(100*cmean(ii)), :)/1.2; % 'k';
    vv(ii).BoxWidth = .03; 
    vv(ii).BoxColor = 'k'; 
    vv(ii).MedianPlot.LineWidth = 2;
    vv(ii).MedianPlot.SizeData = 100;
    vv(ii).MedianPlot.MarkerEdgeColor = 'k';
    vv(ii).EdgeColor = 'none';
%     vv(ii).ViolinAlpha = 0;
%     vv(ii).ScatterPlot.MarkerFaceColor = cm(round(100*cmean(ii)), :)/1.2; % 'k';
    vv(ii).ScatterPlot.MarkerFaceColor = light_colormap(cm(round(100*cmean(ii)), :), 1.5); % 'k';
    vv(ii).ScatterPlot.MarkerFaceAlpha = 1;
    vv(ii).ScatterPlot.SizeData = 20;
    vv(ii).ScatterPlot.MarkerEdgeColor = 'k';
    vv(ii).ShowData = 1;
    vv(ii).ScatterPlot.XData = offs(:,ii) ;
end

set(gca, 'YTick', [0 .2 .4 .6 .8 1], 'YDir', 'normal', 'XTick', [1:4],...
    'XTickLabel', {'Win-Shift', 'Lose-Shift', 'W-L', 'L-L'}, 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold')
ylim([-.05 1])
text(2.1, .85, sprintf('n = %d', sum(valid)))
%%

x = D.probability_correct(valid);
% x = D.EntrPerMin(valid);
s = (D.average_streak(valid));
y = p11-p21;
s(s>7)=7;
figure(98); clf; hold on
xs = linspace(min(x), max(x), 100);
[lreg, ~] = fit(x, y,'poly1','Normalize','on','Robust','Bisquare');
plot(xs, lreg(xs), 'LineWidth', 3, 'Color', [.4 .8 .8])

% s(isnan(s))=0;
scatter(x, y, 200, s, '.')
colormap viridis
% axis([.38 1.01 -.6 .8])
axis([.38 1.01 -.6 .8])
colorbar
set(gca, 'YTick', [-.5:.25:.75], 'YTickLabel', 100*[-.5:.25:.75], 'YDir', 'normal', 'XTick', [.4:.1:1],...
    'XTickLabel', 40:10:100, 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold')
set(gcf, 'Name', 'WS-WL  x  probCorrect  x  av streak', 'Color', 'w', 'Position', [383   506   248   169])

[c, p] = corr(x, p11-p21);
title(sprintf('WS-LS x pCorr  -  corr= %1.2f, p= %1.5f', c, p))

% x = D.probability_correct(valid);
x = D.EntrPerMin(valid);
s = (D.average_streak(valid));
xs = linspace(min(x), max(x), 100);
s(s>7)=7;

figure(97); clf; hold on
[lreg, ~] = fit(x, y,'poly1','Normalize', 'on', 'Robust','Bisquare');
plot(xs, lreg(xs), 'LineWidth', 3, 'Color', [.4 .8 .8])
% s(isnan(s))=0;
scatter(x, y, 200, s, '.')
colormap viridis
% axis([.38 1.01 -.6 .8])
axis([-0 7.1 -.6 .8])
colorbar
set(gca, 'YTick', [-.5:.25:.75], 'YTickLabel', 100*[-.5:.25:.75], 'YDir', 'normal', 'XTick', [0:2:6],...
    'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold')
set(gcf, 'Name', 'WS-WL  x  probCorrect  x  av streak', 'Color', 'w', 'Position', [583   506   248   169])

[c, p] = corr(x, p11-p21);
title(sprintf('WS-LS x R/min  -  corr= %1.2f, p= %1.5f', c, p))



%%
% D = readtable('C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\behav_summary.csv');

valid1 = ~isnan(D.probability_correct) & ~isnan(D.probability_left) & (D.probability_correct ~= D.probability_left);
valid2 = D.DistanceRan>80;
valid=valid1&valid2;
x = D.EntrPerMin(valid);

figure; 
[c, xc] = histcounts(x, [0:.5:7]);
figure; 
plot(xc(1:end-1)+(xc(2)-xc(1)), c)


load('C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL5.mat')
% load('C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24513_IL8.mat')
%%
h = figure(102); clf; 
set(gcf, 'Position', [680   709   275   269])
hold on
x = room.x*-1;  y = room.y;  t = room.timestamps./1000;
% r = x*NaN;
xr = x*NaN;  yr = y*NaN;  tr = t*NaN;
for ii = 1:length(room.entrance_start_idx)
    ab = room.entrance_start_idx(ii)-5:room.entrance_start_idx(ii)+15;
    xr(ab) = x(ab); yr(ab) = y(ab); tr(ab) = t(ab); 
end
ds = 5;
plot3(x(1:ds:end),y(1:ds:end),t(1:ds:end), 'k-')
% for ii = 1:size(zonestruct.choices.Correct_start_end,1)
%     a = zonestruct.choices.Correct_start_end(ii,1):zonestruct.choices.Correct_start_end(ii,2);
%     b = find(x<-20 & y<20 & y>-20);
% %     b = find(abs(zonestruct.arena_zones_seq) == 3);
%     ab = intersect(b,a);
%     xr(ab) = x(ab); yr(ab) = y(ab); tr(ab) = t(ab); 
% end
plot3(xr(1:ds:end),yr(1:ds:end),tr(1:ds:end), 'g-')

x = (-1*arena.x)+100; y = arena.y; t = arena.timestamps./1000;
xr = x*NaN;  yr = y*NaN;  tr = t*NaN;
for ii = 1:length(room.entrance_start_idx)
    ab = room.entrance_start_idx(ii)-5:room.entrance_start_idx(ii)+15;
    xr(ab) = x(ab); yr(ab) = y(ab); tr(ab) = t(ab); 
end
plot3(x(1:ds:end),y(1:ds:end),t(1:ds:end))
% for ii = 1:size(zonestruct.choices.Correct_start_end,1)
%     a = zonestruct.choices.Correct_start_end(ii,1):zonestruct.choices.Correct_start_end(ii,2);
%     b = find(x<-20+100 & y<20 & y>-20);
% %     b = find(abs(zonestruct.arena_zones_seq) == 3);
%     ab = intersect(b,a);
%     xr(ab) = x(ab); yr(ab) = y(ab); tr(ab) = t(ab); 
% end
plot3(xr(1:ds:end),yr(1:ds:end),tr(1:ds:end), 'g-')

% x = arena.x+100; y = arena.y; t = arena.timestamps./1000;
% plot3(x(t<(60*10)),y(t<(60*10)),t(t<(60*10)), 'k-')
% plot3(x(t>(60*10)),y(t>(60*10)),t(t>(60*10)), 'm-')
% for ii = 1:size(zonestruct.choices.Correct_start_end,1)
%     a = zonestruct.choices.Correct_start_end(ii,1):zonestruct.choices.Correct_start_end(ii,2);
%     plot3(x(a),y(a),t(a), 'g-')
% end
    
plot3(-50 + 0*x(t>(60*10)),-30 + 0*y(t>(60*10)),t(t>(60*10)), 'r-')
set(gca, 'ZTick', 600:600:1800, 'XLim', [-55 155], 'YLim', [-55 55], 'View', [11.1000   38.7009], 'Box', 'off', 'ZLim', [0 1800])
plot2svg('C:\Users\gjb326\Desktop\TRACKER DOCS\figures\path_examplesess.svg', h, 'png') 

%%

valid1 = D.arena_displacement_max > 50;
valid2 = ~contains(D.SessName, 'ILCON');
valid3 = contains(D.Name, 'NAPA_24511');
valid = valid1&valid2&valid3;
% valid = D.swap_performance_win_win>0;
xs = D.TRNum(valid);
s1 = D.probability_correct(valid);
y1 = D.swap_performance_win_win(valid);
y2 = D.swap_performance_lose_win(valid);

figure; hold on;
plot(xs, y1);
plot(xs, y2);
ylim([0 1])
set(gca, 'YTick', [0:.5:1], 'XTick', 1:12)
% yyaxis('right')
% plot(xs, s1);
% ylim([.5 1])

