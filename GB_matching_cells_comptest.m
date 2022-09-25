animal_name = 'TestMouse2';%animals{animalLoop};
experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/AlejandroGrau/';
dir_list_fname = '_directory_list.csv';

dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
processedDir        = sprintf('%s%s/processed_files/', experiment_folder, animal_name);
contourDir          = sprintf('%s%s/matching_contours/', experiment_folder, animal_name);

AnimalDir = setup_imaging_Sessionfiles(animal_name, dir_file, processedDir, contourDir);

datadir = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\TestMouse2\processed_files\';
manmatchDir = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\TestMouse2\matching_contours\TestMouse1\';
matchingfname = sprintf('%smatching_matrix.mat', processedDir);
load(matchingfname);
%%
s1 = 1;
s2 = 5;
fname1 = AnimalDir.processedFile{s1};
fname2 = AnimalDir.processedFile{s2};
cname1 = sprintf('%s%s_ms_placecells_data.mat.mat', manmatchDir, AnimalDir.Sessname{s1});
cname2 = sprintf('%s%s_ms_placecells_data.mat.mat', manmatchDir, AnimalDir.Sessname{s2});
c1 = load(cname1);
c2 = load(cname2);

f1 = load(fname1);
f2 = load(fname2);
%
matched = sum(cellmap>0,2)==size(cellmap,2); %cellmap(:,s1)>0 & cellmap(:,s2)>0;
spks1 = f1.ms.neuron.S_matw(cellmap(matched,s1), :);
spks2 = f2.ms.neuron.S_matw(cellmap(matched,s2), :);
spks1 = normalize_rows(spks1);
spks2 = normalize_rows(spks2);

cc1 = c1.contours_shifted(cellmap(matched,s1),:,:);
cc1(cc1<=.07) = 0;
cc2 = c2.contours_shifted(cellmap(matched,s2),:,:);
cc2(cc2<=.07) = 0;
cim1 = squeeze(sum(cc1, 1))>0;
cim2 = squeeze(sum(cc2, 1))>0;
cim = cat(3, cim1, cim1+cim2, cim2);

[pco_tau1, ~, sumspks1] = Fenton_pco(spks1, round(1/median(f1.ms.dt))*5, false, 'Kendall');
[pco_tau2, ~, sumspks2] = Fenton_pco(spks2, round(1/median(f2.ms.dt))*5, false, 'Kendall');

corcoef1 = corr(sumspks1');
corcoef2 = corr(sumspks2');
corcoef1(find(eye(size(corcoef1))==1)) = NaN;
corcoef2(find(eye(size(corcoef2))==1)) = NaN;
figure; 
subplot(131); 
image(cim); %imagesc(corcoef1); 
subplot(132); imagesc(corcoef2)
subplot(133); histogram2(corcoef1(:), corcoef2(:))