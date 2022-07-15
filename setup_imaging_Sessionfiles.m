function out_struct = setup_imaging_Sessionfiles(animal_name, dir_file, processedDir, contourDir)
    %%
    output_behav_file   = 'experiment/ms_behav_data.mat';
    pcellString = 'ms_placecells_data.mat';
    contourString = 'contours.mat';
    output_pcell_file   = strcat('experiment/', pcellString);
    if isempty(dir(processedDir)) == true
        mkdir(processedDir)
    end
    if isempty(dir(contourDir)) == true
        mkdir(contourDir)
    end
    
    temp = readtable(dir_file, 'Delimiter', '&&');
    numSess = size(temp,1);
    out_struct = [];
    out_struct.SessionDir   = temp.SessionDir;
    out_struct.numSess      = numSess;
    for sessionLoop = 1:numSess
        % extract the recording session name - DATE___TIME
        % (YEAR_MONTH_DAY___HOUR_MIN_SEC)
        recording_dir   = out_struct.SessionDir{sessionLoop};
        Sessname        = recording_dir(strfind(recording_dir, animal_name) + length(animal_name) + 1 : end-1);
        s1              = strfind(Sessname, '/');
        Sessname        = strcat(Sessname(1:s1-1), '___', Sessname(s1+1:end));
        alt_pcellName   = strcat(Sessname, '_', pcellString);
        contourFileName = strcat(Sessname, '_', contourString);
        % setup the filenames, full path
        behaviorFile    = sprintf('%s%s', recording_dir, output_behav_file);
        placecellFile   = sprintf('%s%s', recording_dir, output_pcell_file);
        processedFile   = sprintf('%s%s', processedDir,  alt_pcellName);
        contourFile     = sprintf('%s%s', contourDir,  contourFileName);
        % pass the filenames to the data structure
        out_struct.Sessname{sessionLoop, 1}         = Sessname;
        out_struct.behaviorFile{sessionLoop, 1}     = behaviorFile;
        out_struct.placecellFile{sessionLoop, 1}    = placecellFile;
        out_struct.processedFile{sessionLoop, 1}    = processedFile;
        out_struct.contourFile{sessionLoop, 1}      = contourFile;
    end
end