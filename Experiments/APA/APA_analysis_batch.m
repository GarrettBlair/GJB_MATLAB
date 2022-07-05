
params = [];
params.arena_radius             = 40; % in cm
params.pos_bins                 = -40:4:40; % in cm, x and y
params.yaw_bin                  = -pi:pi/8:pi;
params.behav_smoothing_interval = .5; % in seconds, length of smoothing kernel
params.occupancy_thresh         = 1; % in seconds, minimum time in each bin to be counted for place map
params.pfield_kernel_radius     = 3; % kernel ends up being [n*2 + 1] in bins
% params.speed_thresh             = 5; % speed thresh in cm/sec
params.num_partitions           = 2;
params.max_spd_thresh           = 100;
params.min_spd_thresh           = 5;
params.min_samples              = 10;

params.rotate_behav             = true;
params.nan_interp               = true;
params.correct_dt               = true; % correct for large jumps in timestamp file when constructing vmap
params.plotting                 = false;



warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

%%
animals = {'Hipp16942'};
numAnimals = length(animals);
analysis_version = 'v1.0';
dir_list_fname = '_directory_list.csv';
experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/';

for animalLoop = 1:numAnimals
    %%
    animal_name = animals{animalLoop};
    dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
    processedDir        = sprintf('%s%s/processed_files/', experiment_folder, animal_name);
    
    AnimalDir = setupSessionfiles(animal_name, dir_file, processedDir);
    %
    fprintf('\n\nSTART ANALYSIS FOR ANIMAL:       %s\n', animal_name)
    for sessionLoop = 1:AnimalDir.numSess
        behaviorFile        = AnimalDir.behaviorFile{sessionLoop};
        if isempty(dir(behaviorFile)) == true
            fprintf('~~~BEHAV analysis beginning: %s...', AnimalDir.Sessname{sessionLoop})
            recording_dir   = AnimalDir.SessionDir{sessionLoop};
            %
            [ms, behav, params] = APA_generate_data_struct(recording_dir, params);
            save(behaviorFile, 'ms', 'behav', 'params', 'analysis_version', 'AnimalDir');
            %
            fprintf(' Done!\n\t\t%s\n', behaviorFile)
        else
            fprintf('~~~BEHAV analysis skipped: %s\n', AnimalDir.Sessname{sessionLoop})
        end
    end
    for sessionLoop = 1:AnimalDir.numSess
        placecellFile       = AnimalDir.placecellFile{sessionLoop};
        if isempty(dir(placecellFile)) == true
            fprintf('~~~PCELL analysis beginning: %s...', AnimalDir.Sessname{sessionLoop})
            behaviorFile        = AnimalDir.behaviorFile{sessionLoop};
            %
            load(behaviorFile, 'ms');
            [ms]            = APA_generate_placmaps(ms, params);
            save(placecellFile, 'ms', 'params', 'analysis_version', 'AnimalDir');
            %
            fprintf(' Done!\n\t\t%s\n', placecellFile)
        else
            fprintf('~~~PCELL analysis skipped: %s\n', AnimalDir.Sessname{sessionLoop})
        end
    end
    for sessionLoop = 1:AnimalDir.numSess
        processedFile   = AnimalDir.processedFile{sessionLoop};
        if isempty(dir(processedFile)) == true
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
end

function out_struct = setupSessionfiles(animal_name, dir_file, processedDir)
    %%
    output_behav_file   = 'experiment/ms_behav_data.mat';
    pcellname = 'ms_placecells_data.mat';
    output_pcell_file   = strcat('experiment/', pcellname);
    if isempty(dir(processedDir)) == true
        mkdir(processedDir)
    end
    
    temp = readtable(dir_file, 'Delimiter', '&&');
    numSess = size(temp,1);
    out_struct = [];
    out_struct.SessionDir   = temp.SessionDir;
    out_struct.numSess      = numSess;
    for sessionLoop = 1:numSess
        recording_dir = out_struct.SessionDir{sessionLoop};
        Sessname      = recording_dir(strfind(recording_dir, animal_name) + length(animal_name) + 1 : end-1);
        s1            = strfind(Sessname, '/');
        Sessname      = strcat(Sessname(1:s1-1), '___', Sessname(s1+1:end));
        alt_pcellName = strcat(Sessname, '_', pcellname);

        behaviorFile  = sprintf('%s%s', recording_dir, output_behav_file);
        placecellFile  = sprintf('%s%s', recording_dir, output_pcell_file);
        processedFile = sprintf('%s%s', processedDir,  alt_pcellName);
        out_struct.Sessname{sessionLoop, 1}         = Sessname;
        out_struct.behaviorFile{sessionLoop, 1}     = behaviorFile;
        out_struct.placecellFile{sessionLoop, 1}    = placecellFile;
        out_struct.processedFile{sessionLoop, 1}    = processedFile;
    end
end