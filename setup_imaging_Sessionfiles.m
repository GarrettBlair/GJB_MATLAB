function out_struct = setup_imaging_Sessionfiles(animal_name, dir_file, experiment_folder, miniscopeName)
    %%
%     output_behav_file   = 'experiment/ms_behav_data.mat';
    DAT_Dir             = sprintf('%sDAT_files/', experiment_folder);
    animal_DAT_str   = sprintf('%s%s_', DAT_Dir, animal_name);
    if ~isempty(miniscopeName)
        pcellString     = sprintf('@placecells_%s.mat', miniscopeName);
        spikeFileString = sprintf('@spikes_%s.mat', miniscopeName);
        contourString   = sprintf('@contours_%s.mat', miniscopeName);
    else
        miniscopeName = 'MiniLFOV';
        pcellString = '@placecells.mat';
        spikeFileString = '@spikes.mat';
        contourString = 'contours.mat';
    end
    processedDir        = sprintf('%s%s/processed_files/', experiment_folder, animal_name);
    if isempty(dir(processedDir)) == true
        mkdir(processedDir)
    end
    spikeFileDir        = sprintf('%s%s/processed_files/spikeFiles/', experiment_folder, animal_name);
    if isempty(dir(spikeFileDir)) == true
        mkdir(spikeFileDir)
    end
    contourDir          = sprintf('%s%s/matching_contours/', experiment_folder, animal_name);
    if isempty(dir(contourDir)) == true
        mkdir(contourDir)
    end
    
    temp = readtable(dir_file, 'Delimiter', '&&');
    numSess = size(temp,1);
    out_struct = [];
    out_struct.numSess = 0;
    for sessionLoop = 1:numSess
        % extract the recording session name - DATE___TIME
        % (YEAR_MONTH_DAY___HOUR_MIN_SEC_TRIALTYPE)
        recording_dir   = temp.SessionDir{sessionLoop};
        if isfile([recording_dir 'OMIT_REC.txt']) == false
            Sessname        = recording_dir(strfind(recording_dir, animal_name) + length(animal_name) + 1 : end-1);
            s1              = strfind(Sessname, '/');
            if length(Sessname) > 19
                SessType        = Sessname(s1+10:end);
            else
                SessType        = 'None';
            end
            Sessname        = strcat(Sessname(1:s1-1), '_H', Sessname(s1+1:end));
            alt_pcellName   = strcat(Sessname, '_', pcellString);
            spikeFileName   = strcat(Sessname, '_', spikeFileString);
            contourFileName = strcat(Sessname, '_', contourString);
            % setup the filenames, full path
            tracking_room    = sprintf('%s%s_Room.dat', animal_DAT_str, SessType);
            tracking_arena    = sprintf('%s%s_Arena.dat', animal_DAT_str, SessType);
            behavior_saveFile    = sprintf('%sexperiment/%s_@behavior_%s.mat', recording_dir, Sessname, miniscopeName);
            placecellFile   = sprintf('%sexperiment/%s', recording_dir, alt_pcellName);
            processedFile   = sprintf('%s%s', processedDir,  alt_pcellName);
            spikeFile       = sprintf('%s%s', spikeFileDir,  spikeFileName);
            contourFile     = sprintf('%s%s', contourDir,  contourFileName);
            % pass the filenames to the data structure
            out_struct.SessionDir{sessionLoop, 1}           = temp.SessionDir{sessionLoop, 1};
            out_struct.numSess                              = out_struct.numSess + 1;
            out_struct.miniscopeName                        = miniscopeName;
            out_struct.Sessname{sessionLoop, 1}             = Sessname;
            out_struct.SessType{sessionLoop, 1}             = SessType;
            out_struct.tracking_room{sessionLoop, 1}        = tracking_room;
            out_struct.tracking_arena{sessionLoop, 1}       = tracking_arena;
            out_struct.behaviorFile{sessionLoop, 1}         = behavior_saveFile;
            out_struct.placecellFile{sessionLoop, 1}        = placecellFile;
            out_struct.processedFile{sessionLoop, 1}        = processedFile;
            out_struct.spikeFile{sessionLoop, 1}            = spikeFile;
            out_struct.contourFile{sessionLoop, 1}          = contourFile;
            if isfolder(recording_dir)==1
                out_struct.ispath(sessionLoop, 1)           = true;
            else
                out_struct.ispath(sessionLoop, 1)           = false;
            end
        else
            fprintf('Session dir skipped: \n\t%s\n', recording_dir)
        end
    end
end