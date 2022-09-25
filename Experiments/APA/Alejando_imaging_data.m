
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
animals = {'TestMouse1', 'TestMouse2'};
numAnimals = length(animals);
analysis_version = 'v1.0';
dir_list_fname = '_directory_list.csv';
experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/AlejandroGrau/';

% rerun_behav         = false;
% rerun_place         = true;
resave_proccessed   = false;
resave_contours     = false;
for animalLoop = 1:numAnimals
    %%
    animal_name = animals{animalLoop};
    dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
    processedDir        = sprintf('%s%s/processed_files/', experiment_folder, animal_name);
    contourDir          = sprintf('%s%s/matching_contours/', experiment_folder, animal_name);
    
    AnimalDir = setup_imaging_Sessionfiles(animal_name, dir_file, processedDir, contourDir);
    %
    fprintf('\n\nSTART ANALYSIS FOR ANIMAL:       %s\n', animal_name)
    tic    
    for sessionLoop = 1:AnimalDir.numSess
        processedFile = AnimalDir.processedFile{sessionLoop};
        if isempty(dir(processedFile)) == true || resave_proccessed==true
            fprintf('~~~CAIMAN Extraction beginning: %s...', AnimalDir.Sessname{sessionLoop})
            recording_dir   = AnimalDir.SessionDir{sessionLoop};
            if AnimalDir.ispath(sessionLoop) == 1
            %
            [ms, behav, params] = APA_generate_data_struct(recording_dir, params);
            %
            caiman_data = load(ms.caimanFilename);
            [nsegs,nframes] = size(caiman_data.C);
            sminSweepFilename = sprintf('%s/MiniLFOV/deconv_sweep.mat', ms.parentDir);
            [smat, smat_weighted, good_idx, ~] = deconv_sweep_read(sminSweepFilename, params.smin_vals);
            all_good_idx = find(sum(good_idx,1)==size(good_idx,1));
            bad_idx = setdiff(1:size(caiman_data.C,1), all_good_idx);
            caiman_data.idx_components = all_good_idx;
            caiman_data.idx_components_bad = bad_idx;
            temp = sum(smat, 1);
            caiman_data.S_mat = reshape(temp, [nsegs, nframes]);
            temp = sum(smat_weighted, 1);
            caiman_data.S_matw = reshape(temp, [nsegs, nframes]);
            [~, bad_inds, ~, valid_contour_bounds] = Draw_contour_bounding(caiman_data.fullA, ...
                caiman_data.dims, caiman_data.maxFrame, caiman_data.idx_components, false);
            allbad = unique([caiman_data.idx_components_bad, bad_inds']);
            fprintf('\nRemoving %d bad components\n', length(allbad))
            neuron = remove_segments(caiman_data, allbad, false);
            ms.neuron = neuron;
            ms.valid_contour_bounds = valid_contour_bounds;
            
            save(processedFile, 'ms', 'params', 'analysis_version', 'AnimalDir');
            clearvars ms neuron caiman_data
            %
            fprintf(' Done!\n\t\t%s\n', processedFile)
            else
               warning('\n Could not find specified session path\n')
            end
        else
            fprintf('~~~CAIMAN Extraction skipped: %s\n', AnimalDir.Sessname{sessionLoop})
        end
    end
    for sessionLoop = 1:AnimalDir.numSess
        contourFile = AnimalDir.contourFile{sessionLoop};
        processedFile = AnimalDir.processedFile{sessionLoop};
        valid_resave = isempty(dir(contourFile)) == true || resave_contours==true;
        if  valid_resave && exist(processedFile, 'file')==2
            fprintf('~~~RESAVE CONTOURS beginning: %s...', AnimalDir.Sessname{sessionLoop})
            load(processedFile, 'ms');
            %
            contours = gbContours(ms.neuron.fullA, ms.neuron.dims, [], .6);
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
%     APA_LFOV_matching

    toc 
end