expt_name = 'APA_water';
aname = 'Acc20832';
dataDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\processed_files'];
file_str = '*_@placecells.mat';

% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])
% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\hab\'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [1:8], [])
% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\learn\'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [3, 5:13], [])
% saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\newlearn\'];
% manual_contour_matching(dataDir, saveDir, aname, file_str, [], [11:22], [])
saveDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\' expt_name '\' aname '\matching_contours\manual_alignment\all\'];
manual_contour_matching(dataDir, saveDir, aname, file_str, [], [], [])

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

