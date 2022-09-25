function Setup_matching_marix(cellreg_dir, matching_fileName)
%%
% cellreg_dir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\matching_contours\manual_cellreg\';
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
    file_names = A.aligned_data_struct.file_names;
    session_names = A.aligned_data_struct.sessions_list;
    answer = 'Yes';
    if isfile(matching_fileName)
        answer = questdlg('Would you like to overwrite?', ...
            matching_fileName, ...
            'Yes','No','No');
    end
    if strcmp(answer, 'Yes')
        fprintf('\nSaved to %s\n\t%s\n\n', matching_fileName, cellreg_dir)
        save(matching_fileName, 'cellmap', 'file_names', 'session_names');
    else
        fprintf('\nAborting...\n\n')
    end
else
    error('No cell reg file found!')
end
