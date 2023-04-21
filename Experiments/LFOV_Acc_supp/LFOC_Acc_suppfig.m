load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\matching_contours\manual_alignment\learn\cellreg\cellRegistered_20230130_165503.mat')
load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\matching_contours\manual_alignment\learn\cellreg\aligned_data_struct.mat')
load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\matching_contours\manual_alignment\learn\cellreg\modeled_data_struct.mat')
% f1 = "C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\2023_01_08_H14_51_42_TR1_@placecells.mat";
% f2 = "C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\2023_01_09_H14_32_47_TR3_@placecells.mat";
% f3 = "C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\2023_01_11_H14_50_30_TR5_@placecells.mat";
% temp = load(f1);
% ms1 = temp.ms;
% temp = load(f2);
% ms2 = temp.ms;
% temp = load(f3);
% ms3 = temp.ms;
% 
% cmap = cell_registered_struct.cell_to_index_map(:,[3 5 7]);
% shared = find(sum(cmap>0,2)==3);
% segs1 = cmap(shared,1);
% segs2 = cmap(shared,2);
% segs3 = cmap(shared,3);
% c1 = squeeze(sum(aligned_data_struct.spatial_footprints_corrected{3}(segs1,:,:), 1));
% c2 = squeeze(sum(aligned_data_struct.spatial_footprints_corrected{5}(segs2,:,:), 1));
% c3 = squeeze(sum(aligned_data_struct.spatial_footprints_corrected{7}(segs3,:,:), 1));
% % c1 = aligned_data_struct.spatial_footprints_corrected{3}(segs1,:,:);
% % c2 = aligned_data_struct.spatial_footprints_corrected{5}(segs2,:,:);
% % c3 = aligned_data_struct.spatial_footprints_corrected{7}(segs3,:,:);
% im1 = aligned_data_struct.adjusted_footprints_projections{3};
% im2 = aligned_data_struct.adjusted_footprints_projections{5};
% im3 = aligned_data_struct.adjusted_footprints_projections{7};

f1 = "C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\2023_01_09_H13_54_33_TR2_@placecells.mat";
f2 = "C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\2023_01_11_H14_11_36_TR4_@placecells.mat";
f3 = "C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\2023_01_13_H18_48_32_TR6_@placecells.mat";
temp = load(f1);
ms1 = temp.ms;
temp = load(f2);
ms2 = temp.ms;
temp = load(f3);
ms3 = temp.ms;

cmap = cell_registered_struct.cell_to_index_map(:,[4 6 8]);
shared = find(sum(cmap>0,2)==3);
segs1 = cmap(shared,1);
segs2 = cmap(shared,2);
segs3 = cmap(shared,3);
c1 = squeeze(sum(aligned_data_struct.spatial_footprints_corrected{4}(segs1,:,:), 1));
c2 = squeeze(sum(aligned_data_struct.spatial_footprints_corrected{6}(segs2,:,:), 1));
c3 = squeeze(sum(aligned_data_struct.spatial_footprints_corrected{8}(segs3,:,:), 1));
% c1 = aligned_data_struct.spatial_footprints_corrected{3}(segs1,:,:);
% c2 = aligned_data_struct.spatial_footprints_corrected{5}(segs2,:,:);
% c3 = aligned_data_struct.spatial_footprints_corrected{7}(segs3,:,:);
im1 = aligned_data_struct.adjusted_footprints_projections{4};
im2 = aligned_data_struct.adjusted_footprints_projections{6};
im3 = aligned_data_struct.adjusted_footprints_projections{8};

%%

conts1 = gbContours(ms1.neuron.fullA, ms1.neuron.dims, [], .6);
conts2 = gbContours(ms2.neuron.fullA, ms2.neuron.dims, [], .6);
conts3 = gbContours(ms3.neuron.fullA, ms3.neuron.dims, [], .6);


shared_im = cat(3, c1, c2, c3);
all_im = cat(3, im1, im2, im3)>0;

r1 = ms1.neuron.C(segs1,:)+ms1.neuron.YrA(segs1,:);
r2 = ms2.neuron.C(segs2,:)+ms2.neuron.YrA(segs2,:);
r3 = ms3.neuron.C(segs3,:)+ms3.neuron.YrA(segs3,:);

[overlap_mat] = pop_overlap(cmap);
seps2 = [overlap_mat(2,1) overlap_mat(2,3)];
seps4 = [overlap_mat(1,3)];

countmap = histcounts2(modeled_data_struct.neighbors_x_displacements, modeled_data_struct.neighbors_y_displacements, [-9:.25:9], [-9:.25:9]);
modeled_data_struct.normalized_maximal_distance=modeled_data_struct.maximal_distance/2;
x_y_displacements=plot_x_y_displacements(modeled_data_struct.neighbors_x_displacements,modeled_data_struct.neighbors_y_displacements,2,modeled_data_struct.normalized_maximal_distance,modeled_data_struct.number_of_bins,modeled_data_struct.centers_of_bins,pwd,true);

%%
rim1 = ms1.neuron.maxFrame;
rim2 = ms2.neuron.maxFrame;
rim3 = ms3.neuron.maxFrame;

figure(99); clf
subplot(2,4,1);
imshow(rim1, [0 255]);
title(sprintf('Day0 max proj., n=%d', size(ms1.neuron.A,2)))
subplot(2,4,2);
imshow(rim2, [0 255]);
title(sprintf('Day2 max proj., n=%d', size(ms2.neuron.A,2)))
subplot(2,4,3);
imshow(rim3, [0 255]);
title(sprintf('Day4 max proj., n=%d', size(ms3.neuron.A,2)))
subplot(2,4,4);
image(all_im);
title('\color{red}Day0,   \color{green}Day1,   \color{blue}Day2')
axis image off

subplot(2,4,8); cla; hold on
bar([0 2 4], [1 mean(seps2) seps4])
plot([0 2 4], [1 mean(seps2) seps4], 'k')
scatter(0, 1, 'ko', 'MarkerFaceColor', 'k')
scatter([2 2], seps2, 'ko', 'MarkerFaceColor', 'k')
scatter(4, seps4, 'ko', 'MarkerFaceColor', 'k')
axis([-1.5 5.5 0 1])
set(gca, 'XTick', [0 2 4])
xlabel('Session separation (days)')
ylabel('Proportion of cells shared')
% subplot(2,1,2);
% stacked_traces(normalize_rows(r1))

% d=400;
% a1 = make_color_im([], normalize_matrix(c1), c1*0, c1*0);
% a2 = make_color_im(rim2./d, rim2*0, squeeze(sum(conts2,1)), rim2*0);
% a3 = make_color_im(rim3./d, rim3*0, rim3*0, squeeze(sum(conts3,1)));
% subplot(2,4,5);
% image(a1);
% subplot(2,4,6);
% image(c2);
% subplot(2,4,7);
% image(c3);


%%

ms = ms1;
% frames = 1:length(ms3.frameNum);
frames = [100:2:440, 442:11:22*188+42];%length(ms3.frameNum)];
v = ms.fileName; 
cn = normalize_matrix(ms.neuron.C);
figure(8); clf; 
set(gcf, 'Position', [384   264   457   563])

% fname = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\output_figs\acc_example.avi';
fname = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\output_figs\acc_example';
vo = VideoWriter(fname);
vo.Quality = 75;
vo.open()
% cn = double(normalize_matrix(ms.neuron.S));
for i = 1:length(frames)
    %%
    a = ms.neuron.A*cn(:,frames(i));
    c = reshape(a, ms.neuron.dims);
    clf
    set(gcf, 'Color', 'k')
    frame   = single(imread(v, frames(i)));
    dff     = frame-ms.neuron.meanFrame;
    subplot_tight(2,2,1); 
    imshow(frame, [0 255])
%     imshow(frame, [60 240])
    title('\color{white}Raw, motion corrected')
    subplot_tight(2,2,2); 
    imshow(dff, [-20 20])
    title('\color{white}DF/F')
    subplot_tight(2,2,3); cla;
    plot(ms.room.x(frames(1:i)), ms.room.y(frames(1:i)), 'Color', [.5 .3 .3], 'LineWidth', 1.75); hold on
    rectangle('Position', [-42 -42 84 84], 'Curvature', 1, 'EdgeColor', 'w', 'LineWidth',2)
    scatter(ms.room.x(frames(i)), ms.room.y(frames(i)), 60, 'or', 'MarkerFaceColor', 'r')
    t = ms.timestamps(frames(i))./1000 - ms.timestamps(frames(1))./1000;
    title('\color{white}Position')
    text(-45, -45, sprintf('%.1f sec', t), 'Color', 'w', 'FontSize', 14)
    axis([-45 45 -45 45])
    axis square off
    hold off

    subplot_tight(2,2,4)
    imshow(c, [-.025 .1]); %colormap viridis
    title('\color{white}CNMF-E Denoised')
    drawnow
    temp = getframe(gcf);
    im = temp.cdata;
%     vo.writeVideo(im);
end
vo.close()