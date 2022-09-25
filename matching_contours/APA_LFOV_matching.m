% dataDir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\processed_files';
% saveDir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\matching_contours\manual';
% aname = 'Hipp16942';
% file_str = '*_ms_placecells_data.mat';

expt_name = 'APA_water';
aname = 'Hipp18240';
dataDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\processed_files'];
saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual'];
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

