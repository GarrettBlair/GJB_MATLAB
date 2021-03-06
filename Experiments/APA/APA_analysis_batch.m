
params = [];
params.arena_radius             = 40; % in cm
params.pos_bins                 = -40:4:40; % in cm, x and y
params.yaw_bin                  = -pi:pi/8:pi;
params.behav_smoothing_interval = .5; % in seconds, length of smoothing kernel
params.occupancy_thresh         = 1; % in seconds, minimum time in each bin to be counted for place map
params.pfield_kernel_radius     = 3; % kernel ends up being [n*2 + 1] in bins
params.smin_vals                = -50:5:-5; % smin values used to create the 'deconv_sweep.mat'
% params.speed_thresh             = 5; % speed thresh in cm/sec
params.num_partitions           = 2;
params.max_spd_thresh           = 100;
params.min_spd_thresh           = 5;
params.min_samples              = 10;

params.rotate_behav             = true;
params.nan_interp               = true;
params.correct_dt               = true; % correct for large jumps in timestamp file when constructing vmap
params.plotting                 = false;
params.reuse_contour_crop       = 'Crop_params.mat'; % use the previous ms file contour crop, unless empty



warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

%%
animals = {'Hipp16942'};
numAnimals = length(animals);
analysis_version = 'v1.0';
dir_list_fname = '_directory_list.csv';
experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/';

rerun_behav         = false;
rerun_place         = false;
resave_proccessed   = false;
resave_contours     = true;

fit_contours_fullFOV= false;

for animalLoop = 1:numAnimals
    %%
    animal_name = animals{animalLoop};
    dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
    processedDir        = sprintf('%s%s/processed_files/', experiment_folder, animal_name);
    contourDir        = sprintf('%s%s/matching_contours/', experiment_folder, animal_name);
    
    AnimalDir = setup_imaging_Sessionfiles(animal_name, dir_file, processedDir, contourDir);
    %
    fprintf('\n\nSTART ANALYSIS FOR ANIMAL:       %s\n', animal_name)
    tic
    
%     AnimalDir.numSess = 1
    for sessionLoop = 1:AnimalDir.numSess
        behaviorFile        = AnimalDir.behaviorFile{sessionLoop};
        if isempty(dir(behaviorFile)) == true || rerun_behav==true
            fprintf('~~~BEHAV analysis beginning: %s...', AnimalDir.Sessname{sessionLoop})
            recordingDir   = AnimalDir.SessionDir{sessionLoop};
            %
            [ms, behav, params_sub] = APA_generate_data_struct(recordingDir, params);
            ms.params = params_sub;
            save(behaviorFile, 'ms', 'behav', 'analysis_version', 'AnimalDir');
            %
            fprintf(' Done!\n\t\t%s\n', behaviorFile)
        else
            fprintf('~~~BEHAV analysis skipped: %s\n', AnimalDir.Sessname{sessionLoop})
        end
    end
    toc 
    tic
    for sessionLoop = 1:AnimalDir.numSess
        placecellFile       = AnimalDir.placecellFile{sessionLoop};
        if isempty(dir(placecellFile)) == true || rerun_place==true
            fprintf('~~~PCELL analysis beginning: %s...', AnimalDir.Sessname{sessionLoop})
            behaviorFile        = AnimalDir.behaviorFile{sessionLoop};
            %
            load(behaviorFile, 'ms');
            [ms]            = APA_generate_placemaps(ms, params);
            save(placecellFile, 'ms', 'params', 'analysis_version', 'AnimalDir');
            %
            fprintf(' Done!\n\t\t%s\n', placecellFile)
        else
            fprintf('~~~PCELL analysis skipped: %s\n', AnimalDir.Sessname{sessionLoop})
        end
    end
    toc
    tic
    for sessionLoop = 1:AnimalDir.numSess
        processedFile   = AnimalDir.processedFile{sessionLoop};
        if isempty(dir(processedFile)) == true || resave_proccessed==true
            fprintf('~~~PROCESSED file creation: %s...', AnimalDir.Sessname{sessionLoop})
            placecellFile       = AnimalDir.placecellFile{sessionLoop};
            %
            load(placecellFile, 'ms');
            save(processedFile, 'ms', 'params', 'analysis_version', 'AnimalDir');
            %
            fprintf(' Done!\n\t\t%s\n', processedFile)
        else
            fprintf('~~~PROCESSED file creation skipped: %s\n', AnimalDir.Sessname{sessionLoop})
        end
    end
    for sessionLoop = 1:AnimalDir.numSess
        contourFile = AnimalDir.contourFile{sessionLoop};
        processedFile = AnimalDir.processedFile{sessionLoop};
        valid_resave = isempty(dir(contourFile)) == true || resave_contours==true;
        if  valid_resave && exist(processedFile, 'file')==2
            fprintf('~~~RESAVE CONTOURS beginning: %s...', AnimalDir.Sessname{sessionLoop})
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
            fprintf(' Done!\n\t\t%s\n', processedFile)
        elseif valid_resave && ~exist(processedFile, 'file')==2
            fprintf('~~~RESAVE CONTOURS skipped: processed data not found!%s\n', AnimalDir.Sessname{sessionLoop})
        else
            fprintf('~~~RESAVE CONTOURS skipped: %s\n', AnimalDir.Sessname{sessionLoop})
        end
    end
    if false
        APA_LFOV_matching
    end
    toc
end
