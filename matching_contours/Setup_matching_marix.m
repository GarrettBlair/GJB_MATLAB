function Setup_matching_marix(cellreg_dir, matching_fileName, overwrite_flag)
%%
% cellreg_dir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\matching_contours\manual_cellreg\';
if ~isfile(matching_fileName) || overwrite_flag==true
    dd = dir([cellreg_dir 'cellRegistered_*']);
    if ~isempty(dd)
        if length(dd)>1
            currdir = pwd;
            cd(cellreg_dir)
            fn = uigetfile('cellRegistered_*', 'Select the CellReg .mat file');
            cd(currdir)
        else
            fn = dd.name;
        end
        fname = [cellreg_dir fn];
        X = load(fname);
        A = load([cellreg_dir 'aligned_data_struct.mat']);
        cellmap = X.cell_registered_struct.cell_to_index_map;
        file_names = A.aligned_data_struct.file_names';
        session_names = A.aligned_data_struct.sessions_list';
        for i = 1:length(session_names)
            seps = strfind(session_names{i}, '\');
            session_names{i} = session_names{i}(seps(end)+1 : end);
        end
        
        fprintf('\nSaved to %s\n\t%s\n\n', matching_fileName, cellreg_dir)
        save(matching_fileName, 'cellmap', 'file_names', 'session_names');
    else
        error('No cell reg file found!')
    end
else
    fprintf('\nSkipping - %s...\n\n', matching_fileName)
end
