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


run_version = 'v2.01';
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
    for sessind = 1:nfiles
        fullname = [pfiles(sessind).folder '\' pfiles(sessind).name];
        sessname = [pfiles(sessind).name];
        
        load(fullname, 'analysis_version')
        file_ver = str2double(analysis_version(2:end));
        if file_ver < this_ver
            fprintf('REPROCESS %s\n', sessname)
            clearvars ms
            load(fullname, 'ms')
            params = [];
            APA_PKCZmouse_imaging_params_current;
            
            if sessind==1
                plotthisshit = true;
            else
                plotthisshit = false;
            end
            raw_traces = ms.neuron.C + ms.neuron.YrA;
            [dff_filt, ~] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, plotthisshit);
            
            ms.spks = normalize_rows(dff_filt);
            ms.spks(isnan(ms.spks)) = 0;
            clearvars dff_filt raw_filt raw_traces
            
            fprintf('\t\t   [placefields]: \t\t%s  \n', datetime())
            [ms] = APA_generate_polar_placemaps(ms, params);
            
            fprintf('\t\t   [ipos]: \t\t%s  \n', datetime())
            [ms.room.momentary_pos_info, ~, ms.room.conjoint_ipos_min, ms.room.conjoint_ipos_av] = ...
                Fenton_ipos_testing(ms, params.ipos_int_time, 'room_polar',  params);
%             [ms.arena.momentary_pos_info, ~, ms.arena.conjoint_ipos_min, ms.arena.conjoint_ipos_av] = ...
%                 Fenton_ipos_testing(ms, params.ipos_int_time, 'arena_polar', params);
            
            fprintf('\t\t   [svm decode]: \t\t%s  \n', datetime())
            [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
            
            
            fprintf('\t\t   [seg corr]: \t\t%s  \n', datetime())
            spks_1secbin = bin_spks_time(ms.spks>0, 1, ms.timestamps./1000, false);
            spks_1secbin(isnan(spks_1secbin)) = 0;
            [c,p] = corr(spks_1secbin','Type','Kendall');
            segs_corr = [];
            segs_corr.c = c;
            segs_corr.p = p;
            segs_corr.tau = '1sec';
            segs_corr.binspks = spks_1secbin;
            % bin spks to 1 sec, calc Tau corr
            
            analysis_version = run_version;
            save(fullname, 'ms', 'analysis_version', 'segs_corr', '-append')
            fprintf('\t done time: \t\t%s  \n', datetime())
        else
            fprintf('SKIP %s\n', sessname)
        end
    end
    fprintf('\n\n')
end
fprintf('\n PROCESS time: \t\t%s  \n\n', datetime()-fullstart)


%% eval
% need to rerun arena ipos !
time_origin = datetime(2020, 1, 1, 0, 0, 0); %

maxsess = 14;
exp_day = NaN(maxsess,numAnimals);
exp_type = cell(maxsess,numAnimals);
exp_num = NaN(maxsess,numAnimals);
rec_region = zeros(maxsess,numAnimals); % 1==HPC, 2==ACC
numEntr = NaN(maxsess,numAnimals);
sessTime = NaN(maxsess,numAnimals);
entrTime = NaN(maxsess,numAnimals);

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
    for sessind = 1:nfiles
        fullname = [pfiles(sessind).folder '\' pfiles(sessind).name];
        sessname = [pfiles(sessind).name];
        
        load(fullname, 'analysis_version')
        file_ver = str2double(analysis_version(2:end));
        if file_ver == this_ver
            fprintf('EVAL %s\n', sessname)
            clearvars ms
            load(fullname, 'ms')
            [sessDate, trialname, trial_type, trial_num] = recdata_from_parentdir(ms.parentDir);
            exp_day(sessind, animalLoop) = (days(sessDate - time_origin));
            exp_type{sessind, animalLoop} = trial_type;
            exp_num(sessind, animalLoop) = trial_num;%thisind;%
            sessTime(sessind, animalLoop) = ms.timestamps(end) - ms.timestamps(1);
            if strcmp(trial_type, 'HC') == False
                entrTime(sessind, animalLoop) = ms.room.entranceTimes(1);
                numEntr(sessind, animalLoop) = length(ms.room.entranceTimes);
            end
%             try
%                 entrTime(sessind, animalLoop) = ms.room.entranceTimes(1);
%             catch
%                 entrTime(sessind, animalLoop) = -1;
%             end
            numCells(sessind, animalLoop) = size(ms.spks, 1);
            propPlaceR(sessind, animalLoop) = nanmean(ms.room.pcell_stats.infoProb<=.05);
            propPlaceA(sessind, animalLoop) = nanmean(ms.arena.pcell_stats.infoProb<=.05);
            propPlace(sessind, animalLoop) = nanmean(ms.arena.pcell_stats.infoProb<=.05 | ms.room.pcell_stats.infoProb<=.05);
            
            ipos = abs(ms.room.momentary_pos_info) - abs(ms.arena.momentary_pos_info);
            iposm = nanmean(ipos,1);
            if any(iposm)
                SFEP(sessind, animalLoop) = nansum(iposm>0) / (nansum(iposm>0)) || nansum(iposm<0);
            end
        else
            disp('analysis version mismatch')
        end
    end
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    