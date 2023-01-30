% animal_name = 'Acc19947';
animal_name = 'Acc20832';
experiment_folder = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water/';
dir_list_fname = '_directory_list.csv';

dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
processedDir        = sprintf('%s%s/processed_files/', experiment_folder, animal_name);
contourDir          = sprintf('%s%s/matching_contours/', experiment_folder, animal_name);

AnimalDir = setup_imaging_Sessionfiles(animal_name, dir_file, experiment_folder);%DAT_Dir, processedDir, contourDir);


datadir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\' animal_name '\processed_files\'];
manmatchDir = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\' animal_name '\matching_contours\manual_alignment\'];
matchingfname = sprintf('%smatching_matrix.mat', manmatchDir);
load(matchingfname);
%%
nsess = 12;
figure(8); clf; 
figure(9); clf;
figure(10); clf;
[overlap_mat] = pop_overlap(cellmap);
subplot(121);
imagesc(overlap_mat, [.2 .8]); axis square; colorbar
for i = 1:nsess-1
    for j = i+1:nsess
        %%
s1 = i;
s2 = j;
fname1 = AnimalDir.processedFile{s1};
fname2 = AnimalDir.processedFile{s2};
cname1 = sprintf('%s\\%s\\%s_@placecells.mat.mat', manmatchDir, animal_name, AnimalDir.Sessname{s1});
cname2 = sprintf('%s\\%s\\%s_@placecells.mat.mat', manmatchDir, animal_name, AnimalDir.Sessname{s2});
c1 = load(cname1);
c2 = load(cname2);

f1 = load(fname1);
f2 = load(fname2);
%
matched = cellmap(:,s1)>0 & cellmap(:,s2)>0; % sum(cellmap>0,2)==size(cellmap,2); %
% matched = sum(cellmap>0,2)==size(cellmap,2); %
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
figure(8);
subplot_tight(nsess, nsess, (i-1)*nsess + j, [.002 .002]); cla
image(cim); %imagesc(corcoef1); 
axis image off

fs = round(1/median(f1.ms.dt));
integration_t = 3; % sec
% [popcorr, popcorr_prob, cellcorr, cellcorr_prob, pop_sumspks, cell_sumspks] = Fenton_pco(spks, pop_window_size, cell_window_size, sliding_method)
[pco_r1, ~, cellcorr_tau1, tauprob1, sumspks1, temp] = Fenton_pco(spks1, fs*integration_t, fs*integration_t, false);
[pco_r2, ~, cellcorr_tau2, tauprob2, sumspks2] = Fenton_pco(spks2, fs*integration_t, fs*integration_t, false);
sig = tauprob1<=.05 | tauprob2<=.05;
% % corcoef1 = corr(sumspks1');
% % corcoef2 = corr(sumspks2');
% % corcoef1(find(eye(size(corcoef1))==1)) = NaN;
% % corcoef2(find(eye(size(corcoef2))==1)) = NaN;
% % 
% % figure(9); clf 
% % subplot(131); 
% % image(cim); %imagesc(corcoef1); 
% % subplot(132); 
% % imagesc(corcoef2)
% % subplot(133); 
% % histogram2(corcoef1(:), corcoef2(:))



a = cellcorr_tau1(sig);
b = cellcorr_tau2(sig);
gi = ~isnan(a.*b);
corrcorr(s1, s2) = corr(a(gi), b(gi));
corrcorr(s2, s1) = corrcorr(s1, s2);
figure(9);
subplot_tight(nsess-1, nsess-1, (i-1)*(nsess-1) + j - 1, [.002 .002]); cla
scatter(a, b, 1, 'k.')
axis(.5*[-1 1 -1 1])
axis square
drawnow
i
    end
end
figure(10);
subplot(122);
imagesc(corrcorr, [.1 .6]); axis square; colorbar

















