% load('2023_07_25_H17_44_28_TR7_@placecells_HPC_miniscope1.mat')
% params = [];
% APA_PKCZmouse_imaging_params_current;
% 
% params.min_spd_thresh = -1;
% params.pos_bins = [-22 -15:5:15 22];
% params.yaw_bins = -pi:pi/8:pi;
% params.rho_bins = [0 7 14 21];
% % params.rho_bin = [0 10 22];
% 
% mstemp = ms;
% r = ms.room;
% a = ms.arena;
% raw_traces = ms.neuron.C + ms.neuron.YrA;
% [dff_filt, raw_filt] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, false);
% ms.spks = normalize_rows(dff_filt);
% ms.spks(isnan(ms.spks)) = 0;
% 
% figure(10900); clf;
% raw_traces = normalize_rows(raw_traces);
% s = raw_traces;
% s(ms.spks==0) = NaN;
% stacked_traces(raw_traces(1:50,:), .9, {'k-', 'LineWidth', 1});
% stacked_traces(s(1:50,:), .9, {'m-', 'LineWidth', 3});
% 
% 
% % percent place cells, SFEP, cofiring correlations
% [ms] = APA_generate_polar_placemaps(ms, params);
% % [ms] = APA_generate_placemaps(ms, params);
% params.ipos_int_time = .25;
% % [ms.room.momentary_pos_info,  rtemp]  = Fenton_ipos(ms, params.ipos_int_time, 'room', params);
% % [ms.arena.momentary_pos_info, atemp]  = Fenton_ipos(ms, params.ipos_int_time, 'arena', params);
% [ms.room.momentary_pos_info,  rtemp]  = Fenton_ipos_testing(ms, params.ipos_int_time, 'room_polar', params);
% [ms.arena.momentary_pos_info, atemp]  = Fenton_ipos_testing(ms, params.ipos_int_time, 'arena_polar', params);
% [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
% 
% ipos = abs(ms.room.momentary_pos_info) - abs(ms.arena.momentary_pos_info);


run_version = 'v2.1';
this_ver = str2double(run_version(2:end));

fullstart = datetime();
fprintf('\n START time: \t\t%s  \n\n', fullstart)

animals = {'mHPC23454' 'mHPC23459'};
numAnimals = length(animals);

%% reprocess
for animalLoop = 1:numAnimals
    fprintf('\n ANIMAL :  %s\n', animals{animalLoop})
    ddir = ['E:\RecordingData\GarrettBlair\PKCZ_imaging\' animals{animalLoop} '\processed_files\'];
    pfiles = dir([ddir '*_@placecells_HPC*']);
    nfiles = length(pfiles);
    for sessind1 = 1:nfiles
        fullname = [pfiles(sessind1).folder '\' pfiles(sessind1).name];
        sessname = [pfiles(sessind1).name];
        
        load(fullname, 'analysis_version')
        file_ver = str2double(analysis_version(2:end));
        if file_ver < this_ver
            fprintf('REPROCESS %s\n', sessname)
            clearvars ms
            load(fullname, 'ms')
            params = [];
            APA_PKCZmouse_imaging_params_current;
            
            if sessind1==1
                plotthisshit = true;
            else
                plotthisshit = false;
            end
%             raw_traces = ms.neuron.C + ms.neuron.YrA;
%             [dff_filt, ~] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, plotthisshit);
%             
%             ms.spks = normalize_rows(dff_filt);
%             ms.spks(isnan(ms.spks)) = 0;
%             clearvars dff_filt raw_filt raw_traces
%             
%             fprintf('\t\t   [placefields]: \t\t%s  \n', datetime())
%             [ms] = APA_generate_polar_placemaps(ms, params);
%             
%             fprintf('\t\t   [ipos]: \t\t%s  \n', datetime())
%             ms.room.momentary_pos_info = []; ms.room.conjoint_ipos_min = []; ms.room.conjoint_ipos_av = [];
%             ms.arena.momentary_pos_info = []; ms.arena.conjoint_ipos_min = []; ms.arena.conjoint_ipos_av = [];
%             [ms.room.momentary_pos_info, ~, ms.room.conjoint_ipos_min, ms.room.conjoint_ipos_av] = ...
%                 Fenton_ipos_testing(ms, params.ipos_int_time, 'room_polar',  params);
%             [ms.arena.momentary_pos_info, ~, ms.arena.conjoint_ipos_min, ms.arena.conjoint_ipos_av] = ...
%                 Fenton_ipos_testing(ms, params.ipos_int_time, 'arena_polar', params);
%             ipos = abs(ms.room.momentary_pos_info) - abs(ms.arena.momentary_pos_info);
%             iposm = nanmean(ipos,1);
%             figure; imagesc(ipos, [-.1 .1]); colormap redblue;
%             fprintf('\t\t   [svm decode]: \t\t%s  \n', datetime())
%             [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
%             
%             
            fprintf('\t\t   [seg corr]: \t\t%s  \n', datetime())
            spks_1secbin = bin_spks_time(ms.spks>0, 1, ms.timestamps./1000, false);
            spks_1secbin(isnan(spks_1secbin)) = 0;
            corrtype = 'Kendall';
            [c,p] = corr(spks_1secbin','Type', corrtype);
            segs_corr = [];
            segs_corr.c = c;
            segs_corr.p = p;
            segs_corr.type = corrtype;
            segs_corr.tau = '1sec';
            segs_corr.binspks = spks_1secbin;
            
            spks1 = spks_1secbin;
            spks2 = spks1(:, floor(size(spks1,2)/2)+1:end);
            spks1 = spks1(:, 1:floor(size(spks1,2)/2));
            
            [c1,p1] = corr(spks1', 'Type', corrtype);
            [c2,p2] = corr(spks2', 'Type', corrtype);
            segs_corr.split1_corr = c1;
            segs_corr.split1_p = p1;
            segs_corr.split2_corr = c2;
            segs_corr.split2_p = p2;

            
            % bin spks to 1 sec, calc Tau corr
%             disp(size(segs_corr.c))
            analysis_version = run_version;
            save(fullname, 'segs_corr', 'analysis_version', '-append')
            fprintf('\t done time: \t\t%s  \n', datetime())
            clearvars ms analysis_version segs_corr
        else
            fprintf('SKIP %s\n', sessname)
        end
    end
    fprintf('\n\n')
end
fprintf('\n PROCESS time: \t\t%s  \n\n', datetime()-fullstart)


%% eval
time_origin = datetime(2020, 1, 1, 0, 0, 0); %

run_version = 'v2.06';
this_ver = str2double(run_version(2:end));

maxsess = 14;
exp_day = NaN(maxsess,numAnimals);
exp_type = cell(maxsess,numAnimals);
exp_label = cell(maxsess,numAnimals);
exp_num = NaN(maxsess,numAnimals);
rec_region = zeros(maxsess,numAnimals); % 1==HPC, 2==ACC
numEntr = NaN(maxsess,numAnimals);
sessTime = NaN(maxsess,numAnimals);
entrTime = NaN(maxsess,numAnimals);
cmTraveled = NaN(maxsess,numAnimals);

numCells = NaN(maxsess,numAnimals);
propPlaceR = NaN(maxsess,numAnimals);
propPlaceA = NaN(maxsess,numAnimals);
propPlace = NaN(maxsess,numAnimals);
SFEP = NaN(maxsess,numAnimals);
%
for animalLoop = 1:numAnimals
    fprintf('\n ANIMAL :  %s\n', animals{animalLoop})
    ddir = ['E:\RecordingData\GarrettBlair\PKCZ_imaging\' animals{animalLoop} '\processed_files\'];
    pfiles = dir([ddir '*_@placecells_HPC*']);
    nfiles = length(pfiles);
    for sessind1 = 1:nfiles
        %%
        fullname = [pfiles(sessind1).folder '\' pfiles(sessind1).name];
        sessname = [pfiles(sessind1).name];
        
        load(fullname, 'analysis_version')
        file_ver = str2double(analysis_version(2:end));
        if file_ver >= this_ver
            fprintf('EVAL %s\n', sessname)
            clearvars ms
            load(fullname, 'ms')
            [sessDate, trialname, trial_type, trial_num] = recdata_from_parentdir(ms.parentDir);
            exp_day(sessind1, animalLoop) = (days(sessDate - time_origin));
            exp_type{sessind1, animalLoop} = trial_type;
            exp_label{sessind1, animalLoop} = [ trial_type '_' num2str(trial_num)];
            exp_num(sessind1, animalLoop) = trial_num;%thisind;%
            sessTime(sessind1, animalLoop) = (ms.timestamps(end) - ms.timestamps(1))/1000;
            cmTraveled(sessind1, animalLoop) = sum(ms.arena.speed)/sessTime(sessind1, animalLoop);
            if strcmp(trial_type, 'HC') == false
                if isempty(ms.room.entranceTimes)
                entrTime(sessind1, animalLoop) = Inf;
                numEntr(sessind1, animalLoop) = 0;
                else
                entrTime(sessind1, animalLoop) = ms.room.entranceTimes(1);
                numEntr(sessind1, animalLoop) = length(ms.room.entranceTimes);
                end
            end
%             try
%                 entrTime(sessind, animalLoop) = ms.room.entranceTimes(1);
%             catch
%                 entrTime(sessind, animalLoop) = -1;
%             end
            numCells(sessind1, animalLoop) = size(ms.spks, 1);
            propPlaceR(sessind1, animalLoop) = nanmean(ms.room.pcell_stats.infoProb<=.05);
            propPlaceA(sessind1, animalLoop) = nanmean(ms.arena.pcell_stats.infoProb<=.05);
            propPlace(sessind1, animalLoop) = nanmean(ms.arena.pcell_stats.infoProb<=.05 | ms.room.pcell_stats.infoProb<=.05);
            
            ipos = abs(ms.room.momentary_pos_info) - abs(ms.arena.momentary_pos_info);
            iposm = nanmean(ipos,1);
            if any(iposm)
                SFEP(sessind1, animalLoop) = nansum(iposm>0) / (nansum(iposm>0 | iposm<0));
            end
        else
            disp('analysis version mismatch')
        end
    end
end
entrTimesec = entrTime./1000;
numEntrperMin = 60*numEntr./sessTime;
numEntrperMin = numEntrperMin;
%%
valid = strcmp(exp_type, 'TR') | strcmp(exp_type, 'RET') | strcmp(exp_type, 'CON');
valid = sum(valid,2)==2;

% valid = [0 1 1 0 1 1 1 0 0 0 0 0 0 0]' == 1;
valid = [1 1 1 1 1 1 1 0 0 0 0 0 0 0]' == 1;
% valid = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]' == 1;
retsplit = find(strcmp(exp_type(valid), 'RET'), 1)-.5;

figure(191); clf; hold on
plot(numCells(valid,:), '-o', 'MarkerFaceColor', 'w')
% plot(numCells./sessTime, '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-50 2000])
xlim([.5 sum(valid)+.5])
title('Recorded Cells')

figure(192); clf; 
subplot(2,3,1); hold on
plot(numEntrperMin(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)+1], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.5 4])
xlim([.5 sum(valid)+.5])
title('Entr / Min')

subplot(2,3,4); hold on
plot(entrTimesec(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)+1], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-100 1200])
xlim([.5 sum(valid)+.5])
title('First Entry (sec)')
    
    
subplot(2,3,3); hold on
plot(propPlaceA(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)+1], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.05 .25])
xlim([.5 sum(valid)+.5])
title('% place, Arena')
legend({'mHPC23454' 'mHPC23459'}, 'Location', 'northwest')

subplot(2,3,6); hold on
plot(propPlaceR(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)+1], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.05 .25])
xlim([.5 sum(valid)+.5])
title('% place, Room')

subplot(2,3,2); hold on
plot(propPlace(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)+1], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.05 .25])
xlim([.5 sum(valid)+.5])
title('% place, Either')
    
% figure(193); clf; hold on
subplot(2,3,5); hold on
plot([0 sum(valid)+2], [.5 .5], 'k-')
plot(SFEP(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)+1], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.05 1.1])
xlim([.5 sum(valid)+.5])
title('Room SFEP')

    
    
    
%% Cross session
% 
% cellreg_dir = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\matching_contours\manual_alignment_HPC\cellreg\';
% matching_fileName = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\matching_matrix.mat';
% cellreg_dir = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\matching_contours\manual_alignment_HPC\cellreg\';
% matching_fileName = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\matching_matrix.mat';
% Setup_matching_marix(cellreg_dir, matching_fileName);

% valid = [0 1 1 0 0 1 1 0 0 0 0 0 0 0]' == 1;
% valid = [0 1 1 0 0 1 1 0 0 0 0 0 0 0]' == 1;
valid = [1 1 1 1 1 1 1 1 0 0 0 0 0 0]' == 1;

cmap_filenames = {'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\matching_matrix.mat',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\matching_matrix.mat'};
cellsbynumdays = NaN(maxsess, 2);
% sortcorrcorr = NaN(maxsess, maxsess, 2);
% sortcorrp = NaN(maxsess, maxsess, 2);
pfieldcorr_arena = NaN(maxsess, 2);
pfieldcorr_room = NaN(maxsess, 2);
pfieldprob_arena = NaN(maxsess, 2);
pfieldprob_room = NaN(maxsess, 2);
nshared = NaN(maxsess, 2);
sortcorrcorr = NaN(maxsess, 2);
sortcorrp = NaN(maxsess, 2);
mscam_bg = [];
%
corrthresh = .01;
for animalLoop = 1:numAnimals
    fprintf('\n ANIMAL :  %s\n', animals{animalLoop})
    ddir = ['E:\RecordingData\GarrettBlair\PKCZ_imaging\' animals{animalLoop} '\processed_files\'];
    
    pfiles = dir([ddir '*_@placecells_HPC*']);
    nfiles = length(pfiles);
    matchdata = load( cmap_filenames{animalLoop});
    
    cellsbynumdays(:, animalLoop) = histcounts(sum(matchdata.cellmap>0, 2), [1:maxsess+1]);
%     for sessind1 = 1:nfiles-1
    for sessind1 = 6 % TR7 file
        if valid(sessind1)==true
            fullname1 = [pfiles(sessind1).folder '\' pfiles(sessind1).name];
            sessname1 = [pfiles(sessind1).name];
            s1 = load(fullname1, 'segs_corr', 'ms');
            for sessind2 = 1:nfiles % sessind1+1
                %%
                if valid(sessind2)==true && (sessind1~=sessind2)
                    fullname2 = [pfiles(sessind2).folder '\' pfiles(sessind2).name];
                    sessname2 = [pfiles(sessind2).name];
                    disp(' ')
                    disp(sessname1)
                    disp(sessname2)
                    s2 = load(fullname2, 'segs_corr', 'ms');
                    %                     s2 = load(fullname2, 'segs_corr', 'ms');
                    cm = matchdata.cellmap;
                    matched = cm(:, sessind1)>0 & cm(:, sessind2)>0;
                    nshared(sessind2, animalLoop) = sum(matched);
                    
                    s1segs = cm(matched, sessind1);
                    s2segs = cm(matched, sessind2);
                    isupper = triu( true(sum(matched)), 1) == 1;
                    p1 = s1.ms.room.pfields_smooth(s1segs, :, :);
                    p1 = reshape(p1, [size(p1,1), size(p1,2)*size(p1,3)]);
                    p2 = s2.ms.room.pfields_smooth(s2segs, :, :);
                    p2 = reshape(p2, [size(p2,1), size(p2,2)*size(p2,3)]);
                    c = NaN(sum(matched),1); p = NaN(sum(matched),1);
                    for i = 1:sum(matched)
                    [c(i), p(i)] = nancorr(p1(i,:)',p2(i,:)');
                    end
%                     pfieldcorr_room(sessind1, animalLoop) = median(c);
                    pfieldcorr_room(sessind2, animalLoop) = median(c);
                    pfieldprob_room(sessind2, animalLoop) = mean(p<=corrthresh);
                    
                    p1 = s1.ms.arena.pfields_smooth(s1segs, :, :);
                    p1 = reshape(p1, [size(p1,1), size(p1,2)*size(p1,3)]);
                    p2 = s2.ms.arena.pfields_smooth(s2segs, :, :);
                    p2 = reshape(p2, [size(p2,1), size(p2,2)*size(p2,3)]);
                    c = NaN(sum(matched),1); p = NaN(sum(matched),1);
                    for i = 1:sum(matched)
                        [c(i), p(i)] = nancorr(p1(i,:)',p2(i,:)');
                    end
%                     pfieldcorr_arena(sessind1, animalLoop) = median(c);
                    pfieldcorr_arena(sessind2, animalLoop) = median(c);
                    pfieldprob_arena(sessind2, animalLoop) = mean(p<=corrthresh);
                    
                    c1 = s1.segs_corr.c(s1segs, s1segs);
                    c2 = s2.segs_corr.c(s2segs, s2segs);
                    
                    c1 = c1(isupper);
                    c2 = c2(isupper);
                    [~, c1ord] = sort(c1, 'descend');
                    c1 = c1(c1ord);
                    c2 = c2(c1ord);
                    %                     [sortcorrcorr(sessind1, sessind2, animalLoop), sortcorrp(sessind1, sessind2, animalLoop)] = corr(c1, c2,'Type','Pearson');
%                     [sortcorrcorr(sessind1, animalLoop), sortcorrp(sessind1, animalLoop)] = corr(c1, c2,'Type','Pearson');
                    [sortcorrcorr(sessind2, animalLoop), sortcorrp(sessind2, animalLoop)] = corr(c1, c2,'Type','Pearson');
                    %                     [c,p] = corr(spks_1secbin','Type','Kendall');
                    if animalLoop == 1 && sessind1 == 6 && sessind2 == 7
                        c_ex1 = cat(2, c1, c2);
                    end
                    if animalLoop == 2 && sessind1 == 6 && sessind2 == 7
                        c_ex2 = cat(2, c1, c2);
                    end
%                     figure; scatter(c1, c2, 'k.')
                                        clearvars matched s2 cm c1 c2 p1 p2 c
                elseif valid(sessind2)==true && (sessind1==sessind2)
                    c = s1.ms.room.split_corr;
                    p = s1.ms.room.split_p;
                    pfieldcorr_room(sessind2, animalLoop) = median(c);
                    pfieldprob_room(sessind2, animalLoop) = mean(p<=corrthresh);
                    c = s1.ms.room.split_corr;
                    p = s1.ms.room.split_p;
                    
                    pfieldcorr_arena(sessind2, animalLoop) = median(c);
                    pfieldprob_arena(sessind2, animalLoop) = mean(p<=corrthresh);   
                    
                    nshared(sessind2, animalLoop) = size(s1.segs_corr.binspks, 1);
                    
                    spks1 = s1.segs_corr.binspks;
                    spks2 = spks1(:, floor(size(spks1,2)/2)+1:end);
                    spks1 = spks1(:, 1:floor(size(spks1,2)/2));
                    
%                     c1 = corr(spks1', 'Type', 'Kendal');
%                     c2 = corr(spks2', 'Type', 'Kendal');
                    c1 = corr(spks1', 'Type', 'Pearson');
                    c2 = corr(spks2', 'Type', 'Pearson');
                    c1 = c1(isupper);
                    c2 = c2(isupper);
                    [~, c1ord] = sort(c1, 'descend');
                    c1 = c1(c1ord);
                    c2 = c2(c1ord);

                    [sortcorrcorr(sessind2, animalLoop), sortcorrp(sessind2, animalLoop)] = corr(c1, c2,'Type','Pearson');
                    if animalLoop==2; mscam_bg = s1.ms.neuron.meanFrame; end
                end % valid check
            end % sessind1
        end % valid check
        clearvars s1
    end % sessind2
end

animalLoop=2;
ddir = ['E:\RecordingData\GarrettBlair\PKCZ_imaging\' animals{animalLoop} '\processed_files\'];
pfiles = dir([ddir '*_@placecells_HPC*']);
nfiles = length(pfiles);
sessind1=1;
fullname1 = [pfiles(sessind1).folder '\' pfiles(sessind1).name];
s1 = load(fullname1, 'ms');
mscam_bg1 = s1.ms.neuron.maxFrame;

sessind1=nfiles;
fullname1 = [pfiles(sessind1).folder '\' pfiles(sessind1).name];
s1 = load(fullname1, 'ms');
mscam_bg2 = s1.ms.neuron.maxFrame;

%%
% figure(193); clf; hold on

cmap = lbmap(256, 'RedBlue');
figure(194); clf; colormap(cmap)

% valid = [0 1 1 0 0 1 1 0 0 0 0 0 0 0]' == 1;
valid = [1 1 1 1 1 1 1 1 0 0 0 0 0 0]' == 1;
retsplit = find(strcmp(exp_type(valid), 'RET'), 1)-1;


figure(194); clf

subplot(1,4,1);
imagesc(c_ex1, [-.1 .1]); colorbar
title('Ex) mHPC23454')
axis off tight

subplot(1,4,2);
imagesc(c_ex2, [-.1 .1]); colorbar
title('Ex) mHPC23459')
axis off tight


subplot_tight(1,4,3:4, [.2 .15]); cla; hold on
set(gca, 'TickLabelInterpreter', 'none')
plot([1:sum(valid)], sortcorrcorr(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.05 .6])
xlim([.5 sum(valid)+1])
title('Population correlation similarity')
[cv] = corr(c_ex1(:,1), c_ex1(:,2));
text(sum(valid)+2, .3, sprintf('%1.4f', cv), 'Color', cmap(1,:))
 [cv] = corr(c_ex2(:,1), c_ex2(:,2));
text(sum(valid)+2, .2, sprintf('%1.4f', cv), 'Color', cmap(end,:))
    
figure(195); clf

% retsplit = find(strcmp(exp_type(valid), 'RET'), 1)-1; % to mark TR7
subplot(2,2,1); hold on
plot([1:sum(valid)], pfieldcorr_room(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.2 1.1])
xlim([.5 sum(valid)+1])
title('room corr')

subplot(2,2,3); hold on
plot([1:sum(valid)], pfieldcorr_arena(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.2 1.1])
xlim([.5 sum(valid)+1])
title('arena corr')


subplot(2,2,2); hold on
plot([1:sum(valid)], pfieldprob_room(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.1 1.1])
xlim([.5 sum(valid)+1])
title(sprintf('room prop<%0.2f', corrthresh))

subplot(2,2,4); hold on
plot([1:sum(valid)], pfieldprob_arena(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.1 1.1])
xlim([.5 sum(valid)+1])
title(sprintf('arena prop<%0.2f', corrthresh))

    figure(196); clf; hold on
% retsplit = find(strcmp(exp_type(valid), 'RET'), 1)-1; % to mark TR7
% subplot(2,2,1); hold on
plot([1:sum(valid)], nshared(valid,:), '-o', 'MarkerFaceColor', 'w')
plot([retsplit retsplit], [-10000 10000], 'k:')
set(gca, 'XTick', [1:sum(valid)], 'XTickLabel',exp_label(valid,1), 'XTickLabelRotation', 90)
ylim([-.2 1200])
xlim([.5 sum(valid)+1])
title('# of shared cells')
%%
cm = load("E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\matching_contours\manual_alignment_HPC\cellreg\cellRegistered_20231124_151816.mat");
% cs = load("E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\matching_contours\manual_alignment_HPC\cellreg\aligned_data_struct.mat");


cmap = cm.cell_registered_struct.cell_to_index_map;


inds = [1 6 7]; % HAB0 TR7 RET8 sessions
shared = sum(cmap(:, inds)>0, 2) == length(inds);
im = cell(3,1);
px = NaN(sum(shared), length(inds));
py = NaN(sum(shared), length(inds));
for i = 1:3
c1 = cm.cell_registered_struct.spatial_footprints_corrected{inds(i)}; % 
c1 = c1(cmap(shared, inds(i)), :, :);
c2 = reshape(c1, [size(c1,1) size(c1,2)*size(c1,3)]);
c2 = normalize_rows(c2);
c2(c2<.7) = 0; c2 = normalize_rows(c2);
c3 = reshape(c2, [size(c1,1) size(c1,2) size(c1,3)]);
imr = squeeze(sum(c3,1));
im{i} = imr;

pall = cm.cell_registered_struct.centroid_locations_corrected{inds(i)}; % 
pall = pall(cmap(shared, inds(i)), :);
px(:,i) = pall(:,1);
py(:,i) = pall(:,2);
end


% imc = cat(3, im{1}+im{2}, im{2}+im{3}, im{1}+im{3})/2;
%%
imc = cat(3, im{1}, im{2}, im{3})>0;
figure(196); clf;
subplot_tight(2,2,1)
imshow(mscam_bg1./200);
title('Max Projection, Day 1')
subplot_tight(2,2,2)
imshow(mscam_bg2./200);
title(['Max Projection, Day ' num2str(round(exp_day(end,2) - exp_day(1,2)))])
axis image
subplot_tight(2,2,3)
image(imc/1.2); axis image
title('Contours: \color{red}HAB0, \color{green}TR7, \color{blue}RET8')
%
subplot_tight(2,2,4); hold on
cm = jet(sum(shared));
[~, ord] = sort(rand(sum(shared),1));
cm = cm(ord,:);
for i = 1:1:sum(shared)
scatter(px(i,1), py(i,1), 20, 'o', 'MarkerFaceColor', cm(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .7);
scatter(px(i,2), py(i,2), 20, 'o', 'MarkerFaceColor', cm(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .7);
scatter(px(i,3), py(i,3), 20, 'o', 'MarkerFaceColor', cm(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .7);
end
set(gca, 'Color', 'k', 'YDir', 'reverse')
axis image
axis([0 size(imc,2) 0 size(imc,1)])
title('Matched centroids')


    