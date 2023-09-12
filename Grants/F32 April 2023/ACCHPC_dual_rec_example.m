clear

% processedFile_acc  = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_15_H11_24_36_TR12_@placecells_ACC_miniscope2.mat';
% processedFile_hpc  = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_15_H11_24_36_TR12_@placecells_HPC_miniscope1.mat';
processedFile_acc  = "D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_20_H13_11_53_TR15_@placecells_ACC_miniscope2.mat";
processedFile_hpc  = "D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_20_H13_11_53_TR15_@placecells_HPC_miniscope1.mat";

A = load(processedFile_acc, 'ms', 'params');
H = load(processedFile_hpc, 'ms', 'params');
%%
cA = gbContours(normalize_cols( A.ms.neuron.fullA), A.ms.neuron.dims, [], .5);
cH = gbContours(normalize_cols( H.ms.neuron.fullA), H.ms.neuron.dims, [], .5);
na = size(cA,1);
nh = size(cH,1);

% asegs = [1 7 51 69 76 91 112 121];
% hsegs = [1 2 5  7  24 41 133 193];
asegs = [1:8];
hsegs = [1:8];

noisea = movmean(  A.ms.neuron.YrA, 3, 2);
noiseh = movmean(  H.ms.neuron.YrA, 3, 2);
sA = normalize_rows(A.ms.neuron.C+noisea); sA = sA(asegs, 1:10000);
sH = normalize_rows(H.ms.neuron.C+noiseh); sH = sH(hsegs, 1:10000);
for i = 1:size(sA,1)
    sA(i,:) = sA(i,:) + i;
end
for i = 1:size(sH,1)
    sH(i,:) = sH(i,:) + i;
end
%
[~, ord] = sort(rand(na,1));
cm_a = plasma(ceil(na*2.85));
cm_a = cm_a(na*1:na*2, :);
cm_a = cm_a(ord,:);
a_im = color_contours_im(cA, cm_a);
[~, ord] = sort(rand(nh,1));
cm_h = viridis(ceil(nh*3.5));
cm_h = cm_h(nh*1:nh*2, :);
cm_h = cm_h(ord,:);
h_im = color_contours_im(cH, cm_h);


figure(444); clf
set(gcf, 'Color', 'w', 'Units', 'inches', 'OuterPosition', [1 4 8 5])

axisargs = {'Color', 'none', 'Box', 'off', 'FontSize', 8, 'FontName', 'Arial', 'XColor', 'k', 'YColor', 'k'};

subplot_tight(2,3,2, [.055 .025]);
image(h_im+H.ms.neuron.minFrame./285)
x = H.ms.valid_contour_bounds;
axis image off
% axis([x.mw1 x.mw2 x.mh1 x.mh2])
axis([20 260 11 207])
set(gca, axisargs{:})
title(sprintf('Rat#1, HPC cells= %i', nh))

subplot_tight(2,3,5, [.05 .025]);
image(a_im+A.ms.neuron.minFrame./325)
x = A.ms.valid_contour_bounds;
axis image off
axis([x.mw1 x.mw2 x.mh1 x.mh2])
set(gca, axisargs{:})
title(sprintf('Rat#1, ACC cells= %i', na))


ii = round(length(H.ms.timestamps)*1);
e = time_to_indx(H.ms.timestamps, H.ms.room.entranceTimes);
e = e(e<=ii);

subplot_tight(2,3,4, [.125 .025]); hold on
x = H.ms.arena.x(1:ii); y = H.ms.arena.y(1:ii);
plot(x-45,y, 'Color', [0 .2 1]*.8);
scatter(x(e)-45,y(e), 20, 'o', 'MarkerFaceColor', [1 0 .1], 'MarkerEdgeColor', 'k');

x = H.ms.room.x(1:ii); y = H.ms.room.y(1:ii);
plot([20 45 70], [45 0 45], 'k', 'LineWidth', 2);
plot(x+45,y, 'Color', [1 .1 0]*.95);
scatter(x(e)+45,y(e), 20, 'o', 'MarkerFaceColor', [1 0 .1], 'MarkerEdgeColor', 'k');
axis([-90 90 -42 42])
axis equal
title('{\color{blue}Arena} and {\color{red}Room}')
set(gca, 'XTick', [-90 -10 10 90], 'XTickLabel', [-40 40 -40 40], 'YTick', [-40 0 40], axisargs{:})
xlabel('Position, cm')

subplot_tight(1,3,3, [.125 .075]); hold on
t = A.ms.timestamps./60000;
t = t - t(1);
for i = 1:size(sA,1)
    plot(t(1:size(sA,2)), sA(i,:), 'Color', cm_a(i,:), 'LineWidth', .5)
end
t = H.ms.timestamps./60000;
t = t - t(1);
for i = 1:size(sA,1)
    plot(t(1:size(sA,2)), size(sA,1)+sH(i,:), 'Color', cm_h(i,:), 'LineWidth', .5)
end
xlabel('Time (seconds)')
ylabel(sprintf('    ACC cells             HPC cells       '))
title('Example dual recording')
set(gca, 'YTick', [1 8 9 16], 'YTickLabel', [1 8 1 8], axisargs{:},...
'XTick', [ .5 1 1.5 2], 'XTickLabel', [ .5 1 1.5 2]*30, 'XLim', [0 2.5])

















