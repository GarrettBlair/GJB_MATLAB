%%
datadir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\OCAM imaging\Hipp18240\';
cellregfile = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\OCAM imaging\Hipp18240\matching_contours\matching_matrix.mat';
alignfile = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\OCAM imaging\Hipp18240\matching_contours\cellreg\aligned_data_struct.mat';
load(cellregfile);
load(alignfile);
conts = aligned_data_struct.spatial_footprints_corrected;
centers = aligned_data_struct.centroid_locations_corrected;
% pcell_files =  {[datadir 'processed_files\2022_09_28_H17_02_49_TR24_@placecells.mat'],...
%                 [datadir 'processed_files\2022_09_28_H17_31_07_TR25_@placecells.mat'],...
%                 [datadir 'processed_files\2022_09_29_H19_23_01_TR99_@placecells.mat'],...
%                 [datadir 'processed_files\2022_09_30_H18_18_17_TR26_@placecells.mat'],...
%                 [datadir 'processed_files\2022_09_30_H18_50_38_TR27_@placecells.mat']};
cellmap = cellmap(:,[1 3 4]); % for using TR25 99 26
conts = conts([1 3 4]); % for using TR25 99 26
centers = centers([1 3 4]); % for using TR25 99 26
pcell_files =  {[datadir 'processed_files\2022_09_28_H17_02_49_TR24_@placecells.mat'],...
                [datadir 'processed_files\2022_09_29_H19_23_01_TR99_@placecells.mat'],...
                [datadir 'processed_files\2022_09_30_H18_18_17_TR26_@placecells.mat']};

ns = length(pcell_files);
matched = sum(cellmap>0,2) == ns; % for using all 5 sessions; TR24 25 99 26 27
% matched = sum(cellmap(:,1:2)>0,2) == 2 | sum(cellmap(:,2:3)>0,2) == 2; % for using all 5 sessions; TR24 25 99 26 27
% matched = sum(cellmap(:,2:4)>0,2) == 3; % for using TR25 99 26
cmap_ID = find(matched);
sess_ID = cellmap(matched,:);
ncells = sum(matched,1);
%%
corrs       = cell(ns,1);
probs       = cell(ns,1);
spks        = cell(ns,1);
t           = cell(ns,1);
cell_contours    = cell(ns,1);
cell_center = cell(ns,1);

int_time = 1; % seconds for binning
prob_thresh = .05;

figure(1); clf
for i = 1:ns
    temp = load(pcell_files{i});
    ispks = temp.ms.neuron.S_matw(sess_ID(:,i), :);
    ispks = normalize_rows(ispks);
    i_t = temp.ms.timestamps./1000;
    i_t = i_t - i_t(1);
    [i_binspks, ~] = bin_spks_time(ispks, int_time, i_t ,'false');
    [bin_t, ~] = bin_spks_time(i_t', int_time, i_t ,'false');
%     corrs = corr(spks', 'Type', 'Kendall');
    [i_corrs, i_prob] = corr(i_binspks', 'Type', 'Kendall');
    
    i_corrs(find(eye(size(ispks,1)))) = NaN;
    thresh_corrs = i_corrs;
    thresh_corrs(i_prob>prob_thresh) = NaN;
    subplot(1,ns,i);
    imagesc(thresh_corrs, [-.2 .2])
    drawnow
    cent = centers{i}(sess_ID(:,i), :);
    c = conts{i}(sess_ID(:,i), :, :); % zeros(ncells, size(conts,2), size(conts,3));
    for j = 1:ncells
        c(j,:,:) = uint8(normalize_matrix(c(j,:,:))*255);
    end
    corrs{i}    = i_corrs;
    probs{i}    = i_prob;
    spks{i}     = i_binspks;
    t{i}        = bin_t;
    cell_contours{i} = c;
    cell_center{i} = cent;
end

if ns == 5
    r = probs{1}<=prob_thresh | probs{2}<=prob_thresh;
    g = probs{3}<=prob_thresh;
    b = probs{4}<=prob_thresh | probs{5}<=prob_thresh;
elseif ns == 3
    r = probs{1}<=prob_thresh;
    g = probs{2}<=prob_thresh;
    b = probs{3}<=prob_thresh;
end
im = cat(3, squeeze(sum(cell_contours{1},1)), squeeze(sum(cell_contours{2},1)), squeeze(sum(cell_contours{3},1)));
figure(2); image(cat(3, r, g, b))
figure(3); image(im./255)
for i = 1:ncells
    text(cell_center{1}(i,1), cell_center{1}(i,2), 1, num2str(i), 'Color', 'red')
end

vars2save = {'cellmap', 'cell_contours', 'cell_center', 'spks', 'probs', 'corrs', 't', 'pcell_files'};
save('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\OCAM imaging\Hipp18240\OCAM_matched_data.mat', vars2save{:})

























