%%
% topdir = '\\sshfs.r/garrettb@monk.cns.nyu.edu/f/fentonlab/RAWDATA/CaImage/GarrettBlair/ImagingData/'; %PKCZ_imaging/'
topdir = 'E:\RecordingData\GarrettBlair\';
expt_name = 'PKCZ_imaging';
aname = 'mHPC24458';
% aname = 'mHPC23459';
dataDir = [topdir expt_name '\' aname '\processed_files'];
region = 'HPC';
saveDir = [topdir expt_name '\' aname '\matching_contours\manual_alignment_' region '\'];
file_str = ['*_' region '_miniscope*'];

manual_contour_matching(dataDir, saveDir, aname, file_str, [], [4, 12, 14], [2])
cd(saveDir)
% next do cell reg
Setup_matching_marix([saveDir '\cellreg\'], [dataDir '\cellmatching_' region '.mat'])

%% 4/4/2024  -  4/11/2024 HPCACC 24504 24505 24514
expt_name = 'APA';
% aname = 'HPCACC24500';
aname = 'HPCACC24505'; % 'HPCACC24514';
region = 'ACC'; % 'HPC';
% dataDir = ['D:\GarrettBlair\' expt_name '\' aname '\processed_files_OFCON'];
% saveDir = ['D:\GarrettBlair\' expt_name '\' aname '\matching_contours_OFCON\manual_alignment_' region '\'];
dataDir = ['D:\GarrettBlair\' expt_name '\' aname '\processed_files'];
saveDir = ['D:\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment_' region '\'];
% region = 'HPC';
file_str = ['*_' region '_miniscope*'];

manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])
cd(saveDir)
% do this after CellReg alignment
% Setup_matching_marix([saveDir '\cellreg\'], [dataDir '\cellmatching_' region '.mat'])

%%
expt_name = 'PKCZ_imaging';
aname = 'mHPC23454';
% aname = 'mHPC23459';
dataDir = ['E:\RecordingData\GarrettBlair\' expt_name '\' aname '\processed_files'];
region = 'HPC';
saveDir = ['E:\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment_' region '\'];
file_str = ['*_' region '_miniscope*'];

manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], 6)
cd(saveDir)

Setup_matching_marix([saveDir '\cellreg\'], [dataDir '\cellmatching_' region '.mat'])

%%
expt_name = 'APA';
% aname = 'HPCACC24500';
aname = 'HPCACC24503';
% dataDir = ['D:\GarrettBlair\' expt_name '\' aname '\processed_files_OFCON'];
% saveDir = ['D:\GarrettBlair\' expt_name '\' aname '\matching_contours_OFCON\manual_alignment_' region '\'];
dataDir = ['D:\GarrettBlair\' expt_name '\' aname '\processed_files'];
saveDir = ['D:\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment_' region '\'];
% region = 'HPC';
region = 'HPC';
file_str = ['*_' region '_miniscope*'];

manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])
cd(saveDir)

Setup_matching_marix([saveDir '\cellreg\'], [dataDir '\cellmatching_' region '.mat'])
%%
expt_name = 'APA';
% aname = 'HPCACC24500';
aname = 'HPCACC24502';
dataDir = ['D:\GarrettBlair\' expt_name '\' aname '\processed_files'];
file_str = '*_HPC_miniscope1.mat';
% file_str = '*_ACC_miniscope2.mat';

saveDir = ['D:\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment_ACC\'];
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])

%%
expt_name = 'APA';
aname = 'Hipp18240';
dataDir = ['D:\GarrettBlair\' expt_name '\' aname '\processed_files'];
file_str = '*_@placecells.mat';

saveDir = ['D:\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\'];
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])

%%
% expt_name = 'APA_water';
expt_name = 'APA';
aname = 'Acc20832';
% dataDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\processed_files'];
dataDir = ['D:\GarrettBlair\' expt_name '\' aname '\temppro_ret'];
file_str = '*_@placecells.mat';

% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])
% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\hab\'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [1:8], [])
% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\learn\'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [3, 5:13], [])
% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\newlearn\'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [11:22], [])
saveDir = ['D:\GarrettBlair\' expt_name '\' aname '\temppro_ret\manual_alignment\'];
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])
cd(saveDir)

%%
expt_name = 'APA_water';
aname = 'Acc19947';
dataDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\processed_files'];
% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment'];
file_str = '*_@placecells.mat';

% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])
% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\hab\'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [1:8], [])
% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\learn\'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [3, 5:13], [])
saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\newlearn\'];
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [11:23], [])
% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\all\'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])


%%
% dataDir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\processed_files';
% saveDir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\matching_contours\manual';
% aname = 'Hipp16942';
% file_str = '*_ms_placecells_data.mat';

expt_name = 'APA_water';
aname = 'Hipp18240';
dataDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\processed_files'];
saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment'];
file_str = '*_@placecells.mat';

%     sess = sesscomp(i, :); % sessns{i};
% refsess = find(sess==sess(1));
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])

%% Matching cells for Hannah's stress expt, processed files and contours copied from APA_water
expt_name = 'OCAM imaging';
aname = 'Hipp18240';
dataDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\processed_files'];
saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment'];
file_str = '*_@placecells.mat';

%     sess = sesscomp(i, :); % sessns{i};
% refsess = find(sess==sess(1));
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])

%% Matching cells for water manipulation comparisons, processed files and contours copied from APA_water
aname = 'Hipp18240';
dataDir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files';
saveDir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\matching_contours\manual_alignment';
file_str = '*_@placecells.mat';

%     sess = sesscomp(i, :); % sessns{i};
% refsess = find(sess==sess(1));
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])

%%
aname = 'TestMouse2';
dd = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\';
dataDir = sprintf('%s%s\\processed_files', dd, aname);
saveDir = sprintf('%s%s\\matching_contours', dd, aname);
file_str = '*_ms_placecells_data.mat';
%     sess = sesscomp(i, :); % sessns{i};
% refsess = find(sess==sess(1));
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])
%%
aname = 'TestMouse1';
dd = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\';
dataDir = sprintf('%s%s\\processed_files', dd, aname);
saveDir = sprintf('%s%s\\matching_contours', dd, aname);
file_str = '*_ms_placecells_data.mat';
%     sess = sesscomp(i, :); % sessns{i};
% refsess = find(sess==sess(1));
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])

