dataDir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\processed_files';
saveDir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\matching_contours\manual';

aname = 'Hipp16942';
file_str = '*_ms_placecells_data.mat';
%     sess = sesscomp(i, :); % sessns{i};
% refsess = find(sess==sess(1));
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])

%%
% dataDir = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\TestMouse1\processed_files';
% saveDir = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\TestMouse1\matching_contours';
% 
% aname = 'TestMouse1';
% file_str = '*_ms_placecells_data.mat';
% %     sess = sesscomp(i, :); % sessns{i};
% % refsess = find(sess==sess(1));
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])
