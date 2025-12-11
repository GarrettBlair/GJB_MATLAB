h_fn = ["D:\GarrettBlair\APA\HPCACC24505\processed_files\2024_01_15_H15_20_24_HAB0_@placecells_HPC_miniscope1.mat",...
"D:\GarrettBlair\APA\HPCACC24505\processed_files\2024_01_17_H13_16_00_TR3_@placecells_HPC_miniscope1.mat",...
"D:\GarrettBlair\APA\HPCACC24505\processed_files\2024_01_26_H17_44_32_TR6_@placecells_HPC_miniscope1.mat"];

a_fn = ["D:\GarrettBlair\APA\HPCACC24505\processed_files\2024_01_17_H13_16_00_TR3_@placecells_ACC_miniscope2.mat",...
"D:\GarrettBlair\APA\HPCACC24505\processed_files\2024_01_26_H17_44_32_TR6_@placecells_ACC_miniscope2.mat",...
"D:\GarrettBlair\APA\HPCACC24505\processed_files\2024_01_15_H15_20_24_HAB0_@placecells_ACC_miniscope2.mat"];


cmap_hpc = load("D:\GarrettBlair\APA\HPCACC24505\matching_contours\manual_alignment_HPC\cellreg\cellRegistered_20240411_124138.mat");
cmap_acc = load("D:\GarrettBlair\APA\HPCACC24505\matching_contours\manual_alignment_ACC\cellreg\cellRegistered_20240411_122632.mat");
conts_hpc = load("D:\GarrettBlair\APA\HPCACC24505\matching_contours\manual_alignment_HPC\cellreg\aligned_data_struct.mat");
conts_acc = load("D:\GarrettBlair\APA\HPCACC24505\matching_contours\manual_alignment_ACC\cellreg\aligned_data_struct.mat");
%%
% [1, 4, 7] = HAB0, TR3, TR6
sessns = [1,4,7];
c1 = cmap_hpc.cell_registered_struct.cell_to_index_map(:, sessns);

c2 = cmap_acc.cell_registered_struct.cell_to_index_map(:, sessns);

bg_hpc = zeros(conts_hpc.aligned_data_struct.adjusted_y_size, conts_hpc.aligned_data_struct.adjusted_x_size, 3);
bg_acc = zeros(conts_acc.aligned_data_struct.adjusted_y_size, conts_acc.aligned_data_struct.adjusted_x_size, 3);
%%
figure(1234); clf
set(gcf, 'Color', 'w', 'Position', [186   361   885   617])
hpc_roi = [25, 70];
acc_roi = [40, 53];
cropsize = 150;
for snum = 1:3
%     htemp = load(h_fn(snum), 'ms');
%     htemp = htemp.ms.neuron.minFrame/255;
    matched1 = sum(c1>0, 2)==length(sessns);
    nomatched2 = c1(:, snum)>0 & matched1==0;
    hpc_match = conts_hpc.aligned_data_struct.spatial_footprints_corrected{sessns(snum)}(c1(matched1, snum), :, :);
    hpc_match = norm_contours(hpc_match, .7);
    hpc_nomatch = conts_hpc.aligned_data_struct.spatial_footprints_corrected{sessns(snum)}(c1(nomatched2, snum), :, :);
    hpc_nomatch = norm_contours(hpc_nomatch, .7);
    
    hpc_match = squeeze(sum(hpc_match,1));
    hpc_nomatch = squeeze(sum(hpc_nomatch,1));
    bg_hpc = cat(3, hpc_nomatch, hpc_nomatch, hpc_nomatch)/2;
    bg_hpc(bg_hpc(:)>.5)=.5;
    bg_hpc(:,:,2) = bg_hpc(:,:,1) + hpc_match;
    bg_hpc(:,:,3) = bg_hpc(:,:,2);
    %
    matched2 = sum(c2>0, 2)==length(sessns);
    nomatched2 = c2(:, snum)>0 & matched2==0;
    acc_match = conts_acc.aligned_data_struct.spatial_footprints_corrected{sessns(snum)}(c2(matched2, snum), :, :);
    acc_match = norm_contours(acc_match, .7);
    acc_nomatch = conts_acc.aligned_data_struct.spatial_footprints_corrected{sessns(snum)}(c2(nomatched2, snum), :, :);
    acc_nomatch = norm_contours(acc_nomatch, .7);
    
    acc_match = squeeze(sum(acc_match,1));
    acc_nomatch = squeeze(sum(acc_nomatch,1));
    bg_acc = cat(3, acc_nomatch, acc_nomatch, acc_nomatch)/2;
    bg_acc(bg_acc(:)>.5)=.5;
    bg_acc(:,:,1) = bg_acc(:,:,1) + acc_match;
    bg_acc(:,:,3) = bg_acc(:,:,1);
    if cropsize>0
        bg_hpc = bg_hpc(hpc_roi(1):hpc_roi(1)+cropsize, hpc_roi(2):hpc_roi(2)+cropsize, :);
        bg_acc = bg_acc(acc_roi(1):acc_roi(1)+cropsize, acc_roi(2):acc_roi(2)+cropsize, :);
    end
    subplot_tight(2,3,snum, [.01, .01])
    image(bg_hpc)
    axis image off
    subplot_tight(2,3,snum+3, [.01, .01])
    image(bg_acc)
    axis image off
    
end



function conts = norm_contours(conts, thresh)

for i = 1:size(conts,1)
    c = squeeze(conts(i,:,:));
    c = c/max(c(:));
    c(c(:)<thresh) = 0;
    conts(i,:,:) = c;
end

end