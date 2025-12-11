%% HPCACC34990
% topdir = '\\sshfs.r/garrettb@monk.cns.nyu.edu/f/fentonlab/RAWDATA/CaImage/GarrettBlair/ImagingData/'; %PKCZ_imaging/'
topdir = 'D:\GarrettBlair\';
expt_name = 'APA';
aname ='HPCACC34990'; % 'mHPC24458'; % 
% aname = 'mHPC23459';
dataDir = [topdir expt_name '\' aname '\processed_files'];
region = 'ACC';
saveDir = [topdir expt_name '\' aname '\matching_contours\manual_alignment_' region '\'];
file_str = ['*' region '_miniscope*'];

% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])
% cd(saveDir)
% next do cell reg
% CellReg;
% name_name = 'subIL';
name_name = 'RAR';
Setup_matching_marix([saveDir '\cellreg_' name_name '\'], [dataDir '\' name_name '_cellmatching_' region '.mat'], true)
% Setup_matching_marix([saveDir '\cellreg\'], [dataDir '\matching_matrix.mat'], true)

