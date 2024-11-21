% behav var to look into - total path length, 1st/2nd half in water sess,
% inter-entrance times
clear
params = [];
APA_rat_imaging_params_current;
% APA_PKCZmouse_imaging_params_current;

% C:\Users\gjb326\Documents\MATLAB\GJB_fixbanding
% gjb_RUN_fixbanding_demoversion.mat

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

%
% animals = {'Hipp16942'};
% animals = {'Hipp18240'};
% animals = {'HPCACC24500', 'HPCACC24502', 'HPCACC24503'};% animals = {'Acc20832', 'Acc19947'};
animals = {'HPCACC24514'};% animals = {'Acc20832', 'Acc19947'};
experiment_folder = 'E:/GarrettBlair/APA/';
% animals = {'HPCACC24504', 'HPCACC24505'};% animals = {'Acc20832', 'Acc19947'};
% experiment_folder = 'F:/GarrettBlair/APA/';
numAnimals = length(animals);
analysis_version = 'v1.42';
dir_list_fname = '_directory_list.csv';
% experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/';
% experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA/';
% experiment_folder = 'Z:/f/fentonlab/RAWDATA/CaImage/GarrettBlair/ImagingData/APA_HPCACC/';
% experiment_folder = 'D:/APA recordings/';

DAT_Dir             = sprintf('%sDAT_files/', experiment_folder);

rerun_behav         = false;
rerun_place         = false;
resave_proccessed   = false;
resave_contours     = false;
rerun_processed     = false;
resave_simple       = false; %%%%%%%%%%
fit_contours_fullFOV= false;


PRINT_ONLY = false;

%% to add
% speed threshold svm training
% apply svm split by occ to pfields split
% overdispersion calc
%%

for animalLoop = 1:numAnimals
    %%
    animal_name = animals{animalLoop};
    dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
    for cameraLoop = 1:length(params.cameraName)
        cameraName = params.cameraName{cameraLoop};
        AnimalDir = setup_imaging_Sessionfiles(animal_name, dir_file, experiment_folder, cameraName);
        fprintf('\n\nSTART ANALYSIS FOR ANIMAL:       %s\n', animal_name)
        % read in behavior data from tracker
        for sessionLoop = 1:AnimalDir.numSess
            %%
            behaviorFile        = AnimalDir.behaviorFile{sessionLoop};
            camFolderExists = isfolder([AnimalDir.SessionDir{sessionLoop} cameraName]);
            run_this_sess = isempty(dir(behaviorFile)) == true || rerun_behav==true;
            if run_this_sess==true && camFolderExists==true
                %%%%%%%%%%%%%%%%%%%
                fprintf('~~~BEHAV analysis beginning: %s...', AnimalDir.Sessname{sessionLoop})
                fprintf('\n%s\n...', AnimalDir.SessionDir{sessionLoop})
                %
                if PRINT_ONLY==false
%                     disp('a')
                    [ms, behav, params_sub] = APA_generate_data_struct(AnimalDir, sessionLoop, params);
                    ms.params               = params_sub;
                    ms.cameraName           = cameraName;
                    [ms]                    = extract_caiman_data(ms, params, cameraName);
                    save(behaviorFile, 'ms', 'behav', 'analysis_version', 'AnimalDir', '-v7.3');
                end
                fprintf('Done! ->  %s\n', behaviorFile)
                %%%%%%%%%%%%%%%%%%%
            elseif run_this_sess==false && camFolderExists==true
                fprintf('~~~BEHAV analysis skipped: %s\n', AnimalDir.Sessname{sessionLoop})
            elseif run_this_sess==true && camFolderExists==false
                warning('~~~BEHAV analysis skipped, no cam found: \n\t%s\n', [AnimalDir.SessionDir{sessionLoop} cameraName])
            else
                fprintf('~~~BEHAV analysis skipped: %s\n', AnimalDir.Sessname{sessionLoop})
            end
        end % sessionLoop 
        %
    end % cameraLoop
end
%%
for animalLoop = 1:numAnimals
    animal_name = animals{animalLoop};
    dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
    for cameraLoop = 1:length(params.cameraName)
        %%
        cameraName = params.cameraName{cameraLoop};
        AnimalDir = setup_imaging_Sessionfiles(animal_name, dir_file, experiment_folder, cameraName);%DAT_Dir, processedDir, contourDir);
        fprintf('\n\nSTART ANALYSIS FOR ANIMAL:       %s\n', animal_name)
        qq = tic;
        % read in caiman data and construct rate maps
        for sessionLoop = 1:AnimalDir.numSess
            placecellFile       = AnimalDir.placecellFile{sessionLoop};
            camFolderExists = isfolder([AnimalDir.SessionDir{sessionLoop} cameraName]);
            run_this_sess = isempty(dir(placecellFile)) == true || rerun_place==true;
            if run_this_sess==true && camFolderExists==true
                fprintf('~~~PCELL analysis beginning: %s...', AnimalDir.Sessname{sessionLoop})
                behaviorFile        = AnimalDir.behaviorFile{sessionLoop};
                if PRINT_ONLY==false
                    %
                    load(behaviorFile, 'ms');
                    if isfield(ms, 'msBadFramesFile')
                        % remove the bade frames
                        
                    end
                    raw_traces = ms.neuron.C + ms.neuron.YrA;
                    [dff_filt, ~] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, true);
                    
                    ms.spks = normalize_rows(dff_filt);
                    ms.spks(isnan(ms.spks)) = 0;
%                     ms.spks = normalize_rows(ms.neuron.S_mat);
%                     ms.spks(isnan(ms.spks)) = 0;
                    tic
                    [ms] = APA_generate_placemaps(ms, params);
                    toc
%                     [ms.room.momentary_pos_info,  rtemp]  = Fenton_ipos(ms, params.ipos_int_time, 'room', params);
%                     [ms.arena.momentary_pos_info, atemp]  = Fenton_ipos(ms, params.ipos_int_time, 'arena', params);
%                     [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
                    
                    save(placecellFile, 'ms', 'params', 'analysis_version', 'AnimalDir', '-v7.3');
                end
                %
                fprintf(' \nDone!\n\t\t%s\n', placecellFile)
                toc(qq)
                fprintf('\n\n')
            else
                fprintf('~~~PCELL analysis skipped: %s\n', AnimalDir.Sessname{sessionLoop})
            end
        end % sessionLoop
        
        % resave processed files to single directory
        %%
        for sessionLoop = 1:AnimalDir.numSess
            processedFile   = AnimalDir.processedFile{sessionLoop};
            if isempty(dir(processedFile)) == true || resave_proccessed==true
                fprintf('~~~PROCESSED file creation: %s...', AnimalDir.Sessname{sessionLoop})
                placecellFile       = AnimalDir.placecellFile{sessionLoop};
                %
                if PRINT_ONLY==false && isfile(placecellFile)
                    load(placecellFile, 'ms');
                    save(processedFile, 'ms', 'params', 'analysis_version', 'AnimalDir', '-v7.3');
                end
                %
                fprintf(' Done!\n\t\t%s\n', processedFile)
            else
                fprintf('~~~PROCESSED file creation skipped: %s\n', AnimalDir.Sessname{sessionLoop})
            end % resave processed files to single directory
            
%             spikeFile       = AnimalDir.spikeFile{sessionLoop};
%             if isempty(dir(spikeFile)) == true || resave_proccessed==true
%                 fprintf('~~~SPIKES file creation: %s...', AnimalDir.Sessname{sessionLoop})
%                 %
%                 if PRINT_ONLY==false && isfile(processedFile)
%                     load(processedFile, 'ms');
%                     spks = ms.neuron.S_matw;
%                     time = ms.timestamps;
%                     room_x = ms.room.x;
%                     room_y = ms.room.y;
%                     arena_x = ms.arena.x;
%                     arena_y = ms.arena.y;
%                     parentDir = ms.parentDir;
%                     save(spikeFile, 'spks', 'time', 'room_x', 'room_y', 'arena_x', 'arena_y', 'parentDir')
%                 end
%                 %             spks = ms.neuron.S_mat;
%                 %
%                 fprintf(' Done!\n\t\t%s\n', processedFile)
%             else
%                 fprintf('~~~SPIKE file creation skipped: %s\n', AnimalDir.Sessname{sessionLoop})
%             end
        end
        % saving contour files
        for sessionLoop = 1:AnimalDir.numSess
            contourFile = AnimalDir.contourFile{sessionLoop};
            processedFile = AnimalDir.processedFile{sessionLoop};
            valid_resave = isfile(contourFile) == false || resave_contours==true;
            if  valid_resave && isfile(processedFile)
                fprintf('~~~RESAVE CONTOURS beginning: %s...', AnimalDir.Sessname{sessionLoop})
                if PRINT_ONLY==false
                    load(processedFile, 'ms');
                    contours = gbContours(ms.neuron.fullA, ms.neuron.dims, [], .6);
                    % Can load in meta data to refit contours into the camera frame
                    if fit_contours_fullFOV == true
                        recordingDir = AnimalDir.SessionDir{sessionLoop};
                        metaDataFile = sprintf('%sMiniLFOV/metaData.json', recordingDir);
                        [camData] = readJSON(metaDataFile);
                        roiSize = [camData.ROI.height, camData.ROI.width];
                        contour_bounds = ms.params.crop_params.cropROI;
                        [contours] = fit_contous_in_template(contours, contour_bounds, roiSize);
                    end
                    %
                    save(contourFile, 'contours');
                    clearvars ms
                    %
                end
                fprintf(' Done!\n\t\t%s\n', processedFile)
            elseif valid_resave && ~exist(processedFile, 'file')==2
                fprintf('~~~RESAVE CONTOURS skipped: processed data not found!%s\n', AnimalDir.Sessname{sessionLoop})
            else
                fprintf('~~~RESAVE CONTOURS skipped: %s\n', AnimalDir.Sessname{sessionLoop})
            end
        end
        
        fprintf('~~~MAKE SIMPLE FILES beginning: \n')
        for sessionLoop = 1:AnimalDir.numSess
            %%
            processedFile = AnimalDir.processedFile{sessionLoop};
            simpleFile = AnimalDir.simpleFile{sessionLoop};
            if isempty(dir(simpleFile)) == true || resave_simple==true
                fprintf('\t %s...\n', simpleFile)
                simple_file_subfcn(processedFile, simpleFile)
            else
                fprintf(' skipped %s...\n', simpleFile)
            end
        end
        fprintf(' Done!\n\t\t%s\n', processedFile)
        
        if false % USE APA_analysis_batch_reprocess_files.m instead %% rerun_processed == true
        fprintf('~~~PROCESSED file rerunning for: %s...', AnimalDir.Sessname{sessionLoop})
        fprintf(' analysis_version : %s', analysis_version)
        for sessionLoop = 1:AnimalDir.numSess
            processedFile   = AnimalDir.processedFile{sessionLoop};
%             placecellFile       = AnimalDir.placecellFile{sessionLoop};
            if isfile(processedFile)
                fprintf('~~~PROCESSED file rerunning for analysis: %s...', AnimalDir.Sessname{sessionLoop})
                prev_version = load(processedFile, 'analysis_version');
                %
                if strcmp(prev_version.analysis_version, analysis_version) == false
                    clearvars ms
                    load(processedFile, 'ms');
                    
                    raw_traces = ms.neuron.C + ms.neuron.YrA;
                    [dff_filt, ~] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, true);
                    
                    ms.spks = normalize_rows(dff_filt);
                    ms.spks(isnan(ms.spks)) = 0;
                    [ms] = APA_generate_placemaps(ms, params);
                    
                    [ms.room.momentary_pos_info,  rtemp]  = Fenton_ipos(ms, params.ipos_int_time, 'room', params);
                    [ms.arena.momentary_pos_info, atemp]  = Fenton_ipos(ms, params.ipos_int_time, 'arena', params);
                    [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
                                        
                    save(processedFile, 'ms', 'params', 'analysis_version', 'AnimalDir');
                else
                    fprintf('~~~PROCESSED file RERUN skipped: %s\n', AnimalDir.Sessname{sessionLoop})
                    fprintf('\tsame analysis_version : %s\n', prev_version.analysis_version)
                end
                %
                fprintf(' Done!\n\t\t%s\n', processedFile)
            else
                fprintf('~~~PROCESSED file RERUN skipped: %s\n', AnimalDir.Sessname{sessionLoop})
                fprintf('\t(no processed file found) \n')
            end % resave processed files to single directory
        end % % RERUN PROCESSED FILES
        end
    end % cameraLoop
    if false
%%%%%%%         APA_analysis_batch_reprocess_files;
        APA_LFOV_matching;
        cellreg_dir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\matching_contours\manual_alignment\Acc19947\cellreg\';
        matching_fileName = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\matching_contours\manual_alignment\matching_matrix.mat';
        cellreg_dir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\matching_contours\manual_alignment\Acc20832\cellreg\';
        matching_fileName = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\matching_contours\manual_alignment\matching_matrix.mat';
        cellreg_dir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\OCAM imaging\Hipp18240\matching_contours\cellreg\';
        matching_fileName = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\OCAM imaging\Hipp18240\matching_contours\matching_matrix.mat';
        %         cellreg_dir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\matching_contours\cell_reg\';
        %         matching_fileName = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\matching_contours\matching_matrix.mat';
        %         cellreg_dir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\matching_contours\cellreg\';
        %         matching_fileName = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\matching_contours\matching_matrix.mat';
        Setup_matching_marix(cellreg_dir, matching_fileName);
    end
end

function simple_file_subfcn(filename, newname)
%%
% ind1 = strfind(filename, 'processed_files')+length('processed_files') + 1;
% ind2 = strfind(filename, '.mat')+length('.mat') - 1;
% newname = ['D:\GarrettBlair\APA\HPCACC24500\simple_files\' filename(ind1:ind2)];
% load(filename, 'ms')
load(filename, 'ms')

vars2save = {'parentDir' 'filename' 'time_ms' 'dt_ms' 'roomx' 'roomy' 'arenax' 'arenay'...
    'spks' 'craw' 'roomEntranceTimes' 'roomShockTimes' 'arenaEntranceTimes' 'arenaShockTimes'...
    'thisregion' 'sessDate' 'trialname' 'trial_type' 'trial_num'...
    'shock_zone_center' 'shock_zone_size', ...
    'avoid_ind' 'is_entrance' 'shock_zone_dist'...
    'shk_ind' 'entr_ind' 'shock_approach_ind' 'distance_entrance_size'};

parentDir = ms.parentDir;
time_ms = single(ms.timestamps);
dt_ms   = single(ms.dt_corrected/.1000);
roomx   = single(ms.room.x);
roomy   = single(ms.room.y);
arenax  = single(ms.arena.x);
arenay  = single(ms.arena.y);
spks    = single(ms.spks);
craw    = single(normalize_rows(ms.neuron.C+ms.neuron.YrA));


roomEntranceTimes    = single(ms.room.entranceTimes);
roomShockTimes       = single(ms.room.shockTimes);
try
arenaEntranceTimes   = single(ms.arena.entranceTimes);
arenaShockTimes      = single(ms.arena.shockTimes);
catch
arenaEntranceTimes   = [];
arenaShockTimes      = [];
end
% Get some session info
[sessDate, trialname, trial_type, trial_num] = recdata_from_parentdir(ms.parentDir);
if ~isempty(strfind(ms.parentDir, 'Hipp18240/2022_09_08/16_09_56_HAB/'))
    trial_num = 0;
end
if ~isfield(ms, 'cameraName')
    ms.cameraName = 'MiniLFOV';
end
if contains(ms.cameraName, 'ACC_miniscope')
    thisregion='ACC';
elseif contains(ms.cameraName, 'HPC_miniscope')
    thisregion='HPC';
elseif contains(ms.cameraName, 'MiniLFOV') && contains(ms.parentDir, 'Acc')
    thisregion='ACC';
elseif contains(ms.cameraName, 'MiniLFOV') && contains(ms.parentDir, 'Hipp')
    thisregion='HPC';
elseif contains(ms.cameraName, 'MiniLFOV') && contains(ms.parentDir, 'HPC')
    thisregion='HPC';
end

% Calculate shock zone distances
% if strcmp(trial_type, 'TR')
if strcmp(trial_type, 'CON')
    shock_zone_center = 3*pi/2; % typical room shock configuration
    shock_zone_size = pi/6; % size in rad from center to edge
else
    shock_zone_center = pi/2; % typical room shock configuration
    shock_zone_size = pi/6; % size in rad from center to edge
end
distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure

% [shock_approach_ind, is_entrance, ~, ~, shock_zone_dist] = shock_zone_dist_peakfinding(ms.room, ms.timestamps./1000,...
%     shock_zone_center, shock_zone_size, distance_entrance_size, true);
[avoid_ind, is_entrance, ~, ~, shock_zone_dist, shk_ind, entr_ind, shock_approach_ind] = shock_zone_dist_peakfinding(ms.room, time_ms./1000,...
    shock_zone_center, shock_zone_size, distance_entrance_size, false);
% [avoid_ind, is_entrance, ~, ~, shock_zone_dist, shk_ind, entr_ind, shock_approach_ind] = shock_zone_dist_peakfinding(ms.room, time_ms./1000,...
%     shock_zone_center, shock_zone_size, distance_entrance_size, false);
trial_num = single(trial_num);
shock_zone_size = single(shock_zone_size);
shock_zone_center = single(shock_zone_center);
% distance_entrance_size = single(distance_entrance_size);
% shock_zone_dist = single(shock_zone_dist);
% shock_approach_ind = single(shock_approach_ind);
sessDate = char(sessDate);


avoid_ind = single(avoid_ind);
is_entrance = logical(is_entrance);
shock_zone_dist = single(shock_zone_dist);
shk_ind = single(shk_ind);
entr_ind = single(entr_ind);
shock_approach_ind = single(shock_approach_ind);
distance_entrance_size = single(distance_entrance_size);

save(newname, vars2save{:})
end