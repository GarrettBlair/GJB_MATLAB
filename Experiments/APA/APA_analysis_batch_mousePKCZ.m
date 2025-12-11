% behav var to look into - total path length, 1st/2nd half in water sess,
% inter-entrance times

params = [];
APA_PKCZmouse_imaging_params_current

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

%
% animals = {'Hipp16942'};
% animals = {'Hipp18240'};
% animals = {'mHPC23454', 'mHPC23459'};% animals = {'Acc20832', 'Acc19947'};
% animals = {'mHPC23459'};% animals = {'Acc20832', 'Acc19947'};
animals = {'mHPC24458'};%, 'mHPC24459'};% animals = {'Acc20832', 'Acc19947'};
% animals = {'mHPC24459'};% animals = {'Acc20832', 'Acc19947'};
% animals = {'mHPC24458'};% animals = {'Acc20832', 'Acc19947'};
numAnimals = length(animals);
analysis_version = 'v1.61';
dir_list_fname = '_directory_list.csv';
% experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/';
% experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_water/';
% experiment_folder = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\PKCZ_imaging\';
% experiment_folder = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\';
% experiment_folder = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\';

experiment_folder = 'F:\GarrettBlair\APA\PCKZ_imaging\';
% experiment_folder = '\\sshfs.r/garrettb@monk.cns.nyu.edu/f/fentonlab/RAWDATA/CaImage/GarrettBlair/ImagingData/PKCZ_imaging/'
DAT_Dir             = sprintf('%sDAT_files/', experiment_folder);

rerun_behav         = false;
rerun_place         = true;
resave_proccessed   = true;
resave_contours     = false;
fit_contours_fullFOV= false;

% params.parfor_progbar           = true; % whether to show the parfor progress figure

% params.plotting                 = false;
params.reuse_contour_crop = '';

PRINT_ONLY = true;

%% to add
%
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
            behaviorFile        = AnimalDir.behaviorFile{sessionLoop};
            camFolderExists = isfolder([AnimalDir.SessionDir{sessionLoop} cameraName]);
            run_this_sess = isempty(dir(behaviorFile)) == true || rerun_behav==true;
            if run_this_sess==true && camFolderExists==true
                %%%%%%%%%%%%%%%%%%%
                fprintf('~~~BEHAV analysis beginning: %s... [%d]\n', AnimalDir.Sessname{sessionLoop}, sessionLoop)
                fprintf('\n%s\n...', AnimalDir.SessionDir{sessionLoop})
                %
                if PRINT_ONLY==false
                    [ms, behav, params_sub] = APA_generate_data_struct(AnimalDir, sessionLoop, params);
                    ms.params               = params_sub;
                    ms.cameraName           = cameraName;
                    [ms]                    = extract_caiman_data(ms, params, cameraName);
                    save(behaviorFile, 'ms', 'behav', 'analysis_version', 'AnimalDir', '-v7.3');
                else
                    
                end
                fprintf('Done! ->  %s\n', behaviorFile)
                %%%%%%%%%%%%%%%%%%%
            elseif run_this_sess==false && camFolderExists==true
                fprintf('~~~BEHAV analysis skipped: %s\n', AnimalDir.Sessname{sessionLoop})
            elseif run_this_sess==true && camFolderExists==false
                warning('~~~BEHAV analysis skipped, no cam found: \n\t%s\n', [AnimalDir.SessionDir{sessionLoop} cameraName])
            else
                fprintf('~~~BEHAV analysis skipped: %s [%d]\n', AnimalDir.Sessname{sessionLoop}, sessionLoop)
            end
        end % sessionLoop 
        %
    end % cameraLoop
end
%%
for animalLoop = 2%:numAnimals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    animal_name = animals{animalLoop};
    dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
    for cameraLoop = 1:length(params.cameraName)
        %%
        cameraName = params.cameraName{cameraLoop};
        AnimalDir = setup_imaging_Sessionfiles(animal_name, dir_file, experiment_folder, cameraName);%DAT_Dir, processedDir, contourDir);
        fprintf('\n\nSTART ANALYSIS FOR ANIMAL:       %s\n', animal_name)
        qq = tic;
        % read in caiman data and construct rate maps
        for sessionLoop = 11:12%1:AnimalDir.numSess %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            placecellFile       = AnimalDir.placecellFile{sessionLoop};
            camFolderExists = isfolder([AnimalDir.SessionDir{sessionLoop} cameraName]);
            run_this_sess = isempty(dir(placecellFile)) == true || rerun_place==true;
            if run_this_sess==true && camFolderExists==true
                fprintf('~~~PCELL analysis beginning: %s... [%d]\n', AnimalDir.Sessname{sessionLoop}, sessionLoop)
                behaviorFile        = AnimalDir.behaviorFile{sessionLoop};
                if PRINT_ONLY==false
                    %
                    load(behaviorFile, 'ms');
                    
                    raw_traces = ms.neuron.C + ms.neuron.YrA;
                    [dff_filt, ~] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, false);

                    ms.spks = normalize_rows(dff_filt);
                    ms.spks(isnan(ms.spks)) = 0;
                    clearvars dff_filt raw_filt raw_traces
%                     ms.spks = normalize_rows(ms.neuron.S_mat);
%                     ms.spks(isnan(ms.spks)) = 0;

                    [ms] = APA_generate_placemaps(ms, params);
                    
                    [ms.room.momentary_pos_info,  rtemp]  = Fenton_ipos(ms, params.ipos_int_time, 'room', params);
                    [ms.arena.momentary_pos_info, atemp]  = Fenton_ipos(ms, params.ipos_int_time, 'arena', params);
                    [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
                    
                    save(placecellFile, 'ms', 'params', 'analysis_version', 'AnimalDir', '-v7.3');
                end
                %
                fprintf(' \nDone!\n\t\t%s\n', placecellFile)
                toc(qq)
            else
                fprintf('~~~PCELL analysis skipped: %s [%d]\n', AnimalDir.Sessname{sessionLoop}, sessionLoop)
            end
            fprintf('\t\ttime: %s  \n\t\tsession process time -- %.0f seconds\n\n', datetime(), toc(qq))
        end % sessionLoop
        
        % resave processed files to single directory
        %%
%         if resave_proccessed == true
        for sessionLoop = 11:12%1:AnimalDir.numSess
            processedFile   = AnimalDir.processedFile{sessionLoop};
            if isempty(dir(processedFile)) == true || resave_proccessed==true
                fprintf('~~~PROCESSED file creation: %s... [%d]\n', AnimalDir.Sessname{sessionLoop}, sessionLoop)
                placecellFile       = AnimalDir.placecellFile{sessionLoop};
                %
                if PRINT_ONLY==false && isfile(placecellFile)
                    load(placecellFile, 'ms');
                    save(processedFile, 'ms', 'params', 'analysis_version', 'AnimalDir', '-v7.3');
                end
                %
                fprintf('\t%s\n\tDone!\n\n', processedFile)
            else
                fprintf('~~~PROCESSED file creation skipped: %s [%d]\n', AnimalDir.Sessname{sessionLoop}, sessionLoop)
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
%         end
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
                    [ms] = APA_generate_placemaps(ms, params);
                    
                    [ms.room.momentary_pos_info,  rtemp]  = Fenton_ipos(ms, params.ipos_int_time, 'room', params);
                    [ms.arena.momentary_pos_info, atemp]  = Fenton_ipos(ms, params.ipos_int_time, 'arena', params);
                    [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
                                        
                    save(processedFile, 'ms', 'params', 'analysis_version', 'AnimalDir', '-v7.3');
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
        APA_analysis_batch_reprocess_files;
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
