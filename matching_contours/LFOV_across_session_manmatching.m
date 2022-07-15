clear;
% %%
% top_dir = 'F:\ScaryData2\Hipp12to18\';
% save_dir = 'F:\ScaryData2\Hipp12to18\matching_v3\shock\';
% 
% animalNums = [6 7 8 9 12 13 15 18];
% sessns = [];
% sessns{1} = [3:9];
% sessns{2} = [4:10];
% sessns{3} = [1:7];
% sessns{4} = [2:8];
% sessns{5} = [7:13];
% sessns{6} = [7:13];
% sessns{7} = [6:13];
% sessns{8} = [9:17];
% shockSess = [6 7 4 5 10 10 10 13];
% 
% for i = 1:8
%     aname = sprintf('Hipp%d', animalNums(i));
%     file_str = sprintf('%s_linear*', aname);
%     sess = sessns{i};
%     refsess = find(sess==shockSess(i));
%     manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
% end
% 
% %%
% save_dir = 'F:\ScaryData2\Hipp12to18\matching_v3\barier\';
% 
% animalNums = [6 7 8 9 12 13 15 18];
% sessns = [];
% sessns{5} = [3:9];
% sessns{6} = [3:9];
% sessns{7} = [11:17];
% sessns{8} = [4:10];
% barrierSess = [0 0 0 0 6 6 14 7];
% 
% for i = 5:8
%     aname = sprintf('Hipp%d', animalNums(i));
%     file_str = sprintf('%s_linear*', aname);
%     sess = sessns{i};
%     refsess = find(sess==barrierSess(i));
%     manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
% end
% 
% %%
% save_dir = 'F:\ScaryData2\Hipp12to18\matching_v3\scopo\';
% 
% animalNums = [6 7 8 9 12 13 15 18];
% sessns = [];
% sessns{5} = [19:26];
% sessns{6} = [19:26];
% sessns{7} = [19:26];
% sessns{8} = [19:26];
% barrierSess = [0 0 0 0 23 23 23 23];
% 
% for i = 5:8
%     aname = sprintf('Hipp%d', animalNums(i));
%     file_str = sprintf('%s_linear*', aname);
%     sess = sessns{i};
%     refsess = find(sess==barrierSess(i));
%     manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
% end

%% Aligning sessions for the Guo Blair LFOV methods paper
% sesscomp = [ 1 2 7; 1 2 9; 1 2 11; 2 3 8; 4 5 26; 4 5 19; 3 4 14; 15 16 25];
% top_dir = 'F:\ScaryData2\Hipp12to18\';
top_dir = 'F:\GB Minidata\linear_data';
save_dir = 'F:\ScaryData2\Guo Blair Methods data\matching';
% sesscomp = [ 1 2 3; 1 2 3; 1 2 3; 2 3 4; 4 5 6; 4 5 6; 3 4 5; 15 16 17];
sesscomp = [ 8 9 10; 10 11 12 ; 10 11 12; 7 8 9; 15 16 17; 16 17 18; 16 17 18; 18 19 20];
animalNums = [6 7 8 9 12 13 15 18];

for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    file_str = sprintf('%s_linear*', aname);
    sess = sesscomp(i, :); % sessns{i};
    refsess = find(sess==sess(1));
    manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
end

%%
top_dir = 'F:\GB Minidata\linear_data';
save_dir = 'F:\ScaryData2\Guo Blair Methods data\matching';
sesscomp = [ 7 8 9; 13 14 15; 8 9 10; 9 10 11; 9 10 11; 9 10 11; 24 25 26];
animalNums = [30:32, 34:37];

for i = length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    file_str = sprintf('%s_linear*', aname);
    sess = sesscomp(i, :); % sessns{i};
    refsess = find(sess==sess(1));
    manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
end

for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    cmap_name = sprintf('%s\\%s_cmap.mat', save_dir, aname);
    if exist(cmap_name, 'file')~=2
        cellregfname = dir(sprintf('%s\\%s\\cellreg\\cellRegistered*', save_dir, aname));
        temp = load([cellregfname.folder '\' cellregfname.name]);
        cmap = temp.cell_registered_struct.cell_to_index_map;
        
        alignedfname = dir(sprintf('%s\\%s\\cellreg\\aligned_data*', save_dir, aname));
        temp = load([alignedfname.folder '\' alignedfname.name]);
        sessionList = temp.aligned_data_struct.sessions_list;
        sessionNums = NaN(length(sessionList),1);
        for ii = 1:length(sessionList)
            aa = sessionList{ii}(strfind(sessionList{ii}, 'linear')+6:end);
            sessionNums(ii) = str2double(aa(1:strfind(aa, '.mat')-1));
        end
        sessionNums
        fprintf('\n%s\n', cmap_name)
        save(cmap_name, 'cmap', 'sessionList', 'sessionNums')
    else
        fprintf('\n%s found, skipping. \n', cmap_name)
    end
end


%% Aligning sessions for Shock group
% animalNums = [6 7 8 9 12 13 15 18 30 31 32 34 35 36];% 37];
% sess_ref   = [6 7 4 5 10 10 10 13 14 6  20 18 12 16];% 15];
animalNums = [32];% 37];
sess_ref   = [20];% 15];
sessns = [];
sessns{6} = [3:9];
sessns{7} = [4:10];
sessns{8} = [1:10];
sessns{9} = [2:8];
sessns{12} = [7:13];
sessns{13} = [7:13];
sessns{15} = [6:12];
sessns{18} = [10:16];
sessns{30} = [11:17];
sessns{31} = [3:9];
% sessns{32} = [17:23];
sessns{32} = [16:24];
sessns{34} = [13:26];
sessns{35} = [9:15]; % also second shock
sessns{36} = [13:26];
% sessns{37} = [14:26];
% sesscomp = [5:12];%[9:15];
top_dir = 'F:\GB Minidata\linear_data';
save_dir = 'F:\GB Minidata\matching\manmatched_contours_shock';


for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    file_str = sprintf('%s_linear*', aname);
    sess = sessns{animalNums(i)};
    refsess = find(sess==sess_ref(i));
    manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
end

for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    cmap_name = sprintf('%s\\%s_shock_cmap.mat', save_dir, aname);
    if exist(cmap_name, 'file')~=2
        cellregfname = dir(sprintf('%s\\%s\\cellreg\\cellRegistered*', save_dir, aname));
        temp = load([cellregfname.folder '\' cellregfname.name]);
        cmap = temp.cell_registered_struct.cell_to_index_map;
        
        alignedfname = dir(sprintf('%s\\%s\\cellreg\\aligned_data*', save_dir, aname));
        temp = load([alignedfname.folder '\' alignedfname.name]);
        sessionList = temp.aligned_data_struct.sessions_list;
        sessionNums = NaN(length(sessionList),1);
        for ii = 1:length(sessionList)
            aa = sessionList{ii}(strfind(sessionList{ii}, 'linear')+6:end);
            sessionNums(ii) = str2double(aa(1:strfind(aa, '.mat')-1));
        end
        sessionNums
        fprintf('\n%s\n', cmap_name)
        save(cmap_name, 'cmap', 'sessionList', 'sessionNums')
    else
        fprintf('\n%s found, skipping. \n', cmap_name)
    end
end
%% Aligning sessions for Shock 2 group
% sesscomp = [ 1 2 7; 1 2 9; 1 2 11; 2 3 8; 4 5 26; 4 5 19; 3 4 14; 15 16 25];
animalNums = [31 35];% 37];
sess_ref   = [10 12];% 15];
sessns = [];
sessns{31} = [7:13];
sessns{35} = [15:21]; % also second shock
% sessns{37} = [14:26];
% sesscomp = [5:12];%[9:15];
top_dir = 'F:\GB Minidata\linear_data';
save_dir = 'F:\GB Minidata\matching\manmatched_contours_shock2';

% animalNums = [35];

for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    file_str = sprintf('%s_linear*', aname);
    sess = sessns{animalNums(i)};
    refsess = find(sess==sess_ref(i));
    manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
end

for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    cmap_name = sprintf('%s\\%s_shock2_cmap.mat', save_dir, aname);
    if exist(cmap_name, 'file')~=2
        cellregfname = dir(sprintf('%s\\%s\\cellreg\\cellRegistered*', save_dir, aname));
        temp = load([cellregfname.folder '\' cellregfname.name]);
        cmap = temp.cell_registered_struct.cell_to_index_map;
        
        alignedfname = dir(sprintf('%s\\%s\\cellreg\\aligned_data*', save_dir, aname));
        temp = load([alignedfname.folder '\' alignedfname.name]);
        sessionList = temp.aligned_data_struct.sessions_list;
        sessionNums = NaN(length(sessionList),1);
        for ii = 1:length(sessionList)
            aa = sessionList{ii}(strfind(sessionList{ii}, 'linear')+6:end);
            sessionNums(ii) = str2double(aa(1:strfind(aa, '.mat')-1));
        end
        sessionNums
        fprintf('\n%s\n', cmap_name)
        save(cmap_name, 'cmap', 'sessionList', 'sessionNums')
    else
        fprintf('\n%s found, skipping. \n', cmap_name)
    end
end
%% Aligning sessions for Barrier group
% sesscomp = [ 1 2 7; 1 2 9; 1 2 11; 2 3 8; 4 5 26; 4 5 19; 3 4 14; 15 16 25];
animalNums = [12 13 18 30 34 35]; % 35
sessns = [];
sessns{12} = [3:9, 10];
sessns{13} = [3:9, 10];
sessns{18} = [4:10, 13];
sessns{30} = [3:9, 14];
sessns{34} = [5:11, 18];
sessns{35} = [5:11, 12];
sess_ref = [6 6 7 6 8 8];
% sesscomp = [5:12];%[9:15];
top_dir = 'F:\GB Minidata\linear_data';
save_dir = 'F:\GB Minidata\matching\manmatched_contours_barrier';

% animalNums = [35];

for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    file_str = sprintf('%s_linear*', aname);
    sess = sessns{animalNums(i)};
    refsess = find(sess==sess_ref(i));
    manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
end


% for i = 1:length(animalNums)
%     aname = sprintf('Hipp%d', animalNums(i));
%     cellregfname = dir(sprintf('%s\\%s\\cellreg\\cellRegistered*', save_dir, aname));
%     temp = load([cellregfname.folder '\' cellregfname.name]);
%     cmap_name = sprintf('%s\\%s_barrer_cmap.mat', save_dir, aname);
%     cmap = temp.cell_registered_struct.cell_to_index_map;
%     
%     alignedfname = dir(sprintf('%s\\%s\\cellreg\\aligned_data*', save_dir, aname));
%     temp = load([alignedfname.folder '\' alignedfname.name]);
%     sessionList = temp.aligned_data_struct.sessions_list;
%     
%     fprintf('\n%s', cmap_name)
%     save(cmap_name, 'cmap', 'sessionList')
% end



%% Aligning sessions for Scoposhock group
% sesscomp = [ 1 2 7; 1 2 9; 1 2 11; 2 3 8; 4 5 26; 4 5 19; 3 4 14; 15 16 25];
% animalNums = [30 31 32 34 36]; % 35
animalNums = [12 13 15 18]; % 35
sessns = [];
sessns{12} = [20:26, 10];
sessns{13} = [20:26, 10];
sessns{15} = [20:26, 9, 10];
sessns{18} = [20:26, 13];
% sessns{30} = [7:13, 14];
% sessns{31} = [17:23, 6];
% sessns{32} = [8:14, 20];
% sessns{34} = [9:15, 18];
% sessns{36} = [9:15, 16];
% sess_ref = [14 6 20 18 16]; % 30-36
sess_ref = [23 23 23 23];
% sesscomp = [5:12];%[9:15];
top_dir = 'F:\GB Minidata\linear_data';
save_dir = 'F:\GB Minidata\matching\manmatched_contours_scoposhock';

% animalNums = [35];

for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    file_str = sprintf('%s_linear*', aname);
    sess = sessns{animalNums(i)};
    refsess = find(sess==sess_ref(i));
    manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
end
%%
for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    cellregfname = dir(sprintf('%s\\%s\\cellreg\\cellRegistered*', save_dir, aname));
    temp = load([cellregfname.folder '\' cellregfname.name]);
    cmap_name = sprintf('%s\\%s_scopshock_cmap.mat', save_dir, aname);
    cmap = temp.cell_registered_struct.cell_to_index_map;
    
    alignedfname = dir(sprintf('%s\\%s\\cellreg\\aligned_data*', save_dir, aname));
    temp = load([alignedfname.folder '\' alignedfname.name]);
    sessionList = temp.aligned_data_struct.sessions_list;
    sessionNums = NaN(length(sessionList),1);
    for ii = 1:length(sessionList)
        aa = sessionList{ii}(strfind(sessionList{ii}, 'linear')+6:end);
        sessionNums(ii) = str2double(aa(1:strfind(aa, '.mat')-1));
    end
    sessionNums
    fprintf('\n%s\n', cmap_name)
    save(cmap_name, 'cmap', 'sessionList', 'sessionNums')
end

%% Aligning sessions for Scopo group
% sesscomp = [ 1 2 7; 1 2 9; 1 2 11; 2 3 8; 4 5 26; 4 5 19; 3 4 14; 15 16 25];
animalNums = [30 31 32 36]; % 35
sessns = [];
sessns{30} = [20:26, 14];
sessns{31} = [13:19, 20];
sessns{32} = [4:10, 18];
sessns{36} = [5:11, 16];
sess_ref = [14 20 18 16];
% sesscomp = [5:12];%[9:15];
top_dir = 'F:\GB Minidata\linear_data';
save_dir = 'F:\GB Minidata\matching\manmatched_contours_scopo';

% animalNums = [35];

for i = 3%length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    file_str = sprintf('%s_linear*', aname);
    sess = sessns{animalNums(i)};
    refsess = find(sess==sess_ref(i));
    manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
end

for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    cellregfname = dir(sprintf('%s\\%s\\cellreg\\cellRegistered*', save_dir, aname));
    temp = load([cellregfname.folder '\' cellregfname.name]);
    cmap_name = sprintf('%s\\%s_scopo_cmap.mat', save_dir, aname);
    cmap = temp.cell_registered_struct.cell_to_index_map;
    
    alignedfname = dir(sprintf('%s\\%s\\cellreg\\aligned_data*', save_dir, aname));
    temp = load([alignedfname.folder '\' alignedfname.name]);
    sessionList = temp.aligned_data_struct.sessions_list;
    sessionNums = NaN(length(sessionList),1);
    for ii = 1:length(sessionList)
        aa = sessionList{ii}(strfind(sessionList{ii}, 'linear')+6:end);
        sessionNums(ii) = str2double(aa(1:strfind(aa, '.mat')-1));
    end
    sessionNums
    fprintf('\n%s\n', cmap_name)
    save(cmap_name, 'cmap', 'sessionList', 'sessionNums')
end


%% Aligning sessions for Drug reexposure group
% sesscomp = [ 1 2 7; 1 2 9; 1 2 11; 2 3 8; 4 5 26; 4 5 19; 3 4 14; 15 16 25];
animalNums = [31 32 36]; % 35
sessns = [];
sessns{31} = [14 16 18 20];
sessns{32} = [5 7 9 11];
sessns{36} = [6 8 10 12];
sess_ref = [16 7 8];
% sesscomp = [5:12];%[9:15];
top_dir = 'F:\GB Minidata\linear_data';
save_dir = 'F:\GB Minidata\matching\manmatched_contours_statedep';

% animalNums = [35];

for i =1:2%length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    file_str = sprintf('%s_linear*', aname);
    sess = sessns{animalNums(i)};
    refsess = find(sess==sess_ref(i));
    manual_contour_matching(top_dir, save_dir, aname, file_str, [], sess, refsess)
end

for i = 1:length(animalNums)
    aname = sprintf('Hipp%d', animalNums(i));
    cellregfname = dir(sprintf('%s\\%s\\cellreg\\cellRegistered*', save_dir, aname));
    temp = load([cellregfname.folder '\' cellregfname.name]);
    cmap_name = sprintf('%s\\%s_statedep_cmap.mat', save_dir, aname);
    cmap = temp.cell_registered_struct.cell_to_index_map;
    
    alignedfname = dir(sprintf('%s\\%s\\cellreg\\aligned_data*', save_dir, aname));
    temp = load([alignedfname.folder '\' alignedfname.name]);
    sessionList = temp.aligned_data_struct.sessions_list;
    sessionNums = NaN(length(sessionList),1);
    for ii = 1:length(sessionList)
        aa = sessionList{ii}(strfind(sessionList{ii}, 'linear')+6:end);
        sessionNums(ii) = str2double(aa(1:strfind(aa, '.mat')-1));
    end
    sessionNums
    fprintf('\n%s\n', cmap_name)
    save(cmap_name, 'cmap', 'sessionList', 'sessionNums')
end




















