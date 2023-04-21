clear

% a = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\ipos_files\Acc20832_TR6_int_0.25sec.mat");
% ddir2 = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\' a.fname];
% ms_a = load(ddir2);
% ad = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\ipos_decoding\Acc20832_TR3_decoding_0.25sec.mat");
% bins = a.params.pos_bins;
% a_vmap = a.ms_room_temp.vmap;
h = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\ipos_files\Hipp18240_TR3_int_0.25sec.mat");
a = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\ipos_files\Acc19947_TR3_int_0.25sec.mat");
ddir1 = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\' h.fname];
ddir2 = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\' a.fname];
ms_h = load(ddir1);
ms_a = load(ddir2);
ad = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\ipos_decoding\Acc19947_TR3_decoding_0.25sec.mat");
hd = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\ipos_decoding\Hipp18240_TR3_decoding_0.25sec.mat");

iposmean = nanmean(a.ipos_smooth,1);
[acc_ipos_map, avmap]  = make_occupancymap_2D(a.ms_room_temp.x, a.ms_room_temp.y, iposmean, bins, bins);
figure(108); clf
imagesc(a.ipos_ra, [-1 1])
axis image
imagesc(a.ipos_ra, [-1 1])
axis image

figure(97); clf
spks = a.ipos_smooth;%ms_h.ms.neuron.C + ms_h.ms.neuron.YrA/3;
spks = a.ipos_ra;%ms_h.ms.neuron.C + ms_h.ms.neuron.YrA/3;
t = a.ms_room_temp.t;%    ms_h.ms.timestamps./1000;% + ms_h.ms.neuron.YrA/4;
spks = zscore(spks, 1, 2);
% spks = spks([1:6, 19:27],1:600*11/4);
tt = mod(t, 60);
tt = find(tt(2:end) < tt(1:end-1))+1;
% stacked_traces(spks, 3, {'c'})
imagesc(spks, [-2 2])
set(gca, 'XTick', tt, 'XTickLabel', round(t(tt)))
set(gcf,'units','centimeters','position',[8,8,12,8]);
set(gca,'FontSize',10,'FontName','Arial');


%%
clear
h = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\ipos_files\Hipp18240_TR3_int_0.25sec.mat");
a = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\ipos_files\Acc19947_TR3_int_0.25sec.mat");
ddir1 = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\' h.fname];
ddir2 = ['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\' a.fname];
ms_h = load(ddir1);
ms_a = load(ddir2);
ad = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\ipos_decoding\Acc19947_TR3_decoding_0.25sec.mat");
hd = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\ipos_decoding\Hipp18240_TR3_decoding_0.25sec.mat");

%%
bins = h.params.pos_bins;
h_vmap = h.ms_room_temp.vmap;
a_vmap = a.ms_room_temp.vmap;

[hipp_ipos_map, hvmap] = make_occupancymap_2D(h.ms_room_temp.x, h.ms_room_temp.y, nanmean(h.ipos_ra,1), bins, bins);
[acc_ipos_map, avmap]  = make_occupancymap_2D(a.ms_room_temp.x, a.ms_room_temp.y, nanmean(a.ipos_ra,1), bins, bins);
% figure(108); clf
% imagesc(hipp_ipos_map./hvmap)
% imagesc(acc_ipos_map./avmap)

figure(86); clf
im1 = ms_h.ms.neuron.maxFrame./255;
[conts, conts_outline] = contour_outlines(ms_h.ms.neuron.fullA, .35, .65, ms_h.ms.neuron.dims);
conts_outline = conts_outline/3;
conts = conts./10;
im1 = cat(3, im1+conts_outline, im1+conts+conts_outline, im1+conts+conts_outline);
im1 = im1(94:360, 123:420,:);
image(im1); axis image off
set(gcf,'units','centimeters','position',[8,12,8,8]);
set(gca,'FontSize',10,'FontName','Arial');
% export_fig 'D:\F32 FIGS\hipp_conts' -eps


figure(87); clf
spks = ms_h.ms.neuron.C + ms_h.ms.neuron.YrA/3;
t = ms_h.ms.timestamps./1000;% + ms_h.ms.neuron.YrA/4;
spks = normalize_rows(spks);
spks = spks([1:6, 19:27],1:600*11);
tt = mod(t, 60);
tt = find(tt(2:end) < tt(1:end-1))+1;
stacked_traces(spks, 3, {'c'})
set(gca, 'XTick', tt, 'XTickLabel', round(t(tt)))
set(gcf,'units','centimeters','position',[8,4,12,8]);
set(gca,'FontSize',10,'FontName','Arial');
% export_fig 'D:\F32 FIGS\hipp_traces' -eps

figure(88); clf
im2 = ms_a.ms.neuron.maxFrame./255;
[conts, conts_outline] = contour_outlines(ms_a.ms.neuron.fullA, .35, .65, ms_a.ms.neuron.dims);
conts_outline = conts_outline/3;
conts = conts./10;
im2 = cat(3, im2+conts+conts_outline, im2+conts_outline, im2+conts+conts_outline);
im2 = im2(33:283, 23:302,:);
image(im2); axis image off
set(gcf,'units','centimeters','position',[16,12,8,8]);
set(gca,'FontSize',10,'FontName','Arial');
% export_fig 'D:\F32 FIGS\acc_conts' -eps

figure(89); clf
spks = ms_a.ms.neuron.C + ms_a.ms.neuron.YrA/3;
t = ms_a.ms.timestamps./1000;% + ms_h.ms.neuron.YrA/4;
spks = normalize_rows(spks);
spks = spks([23:35, 41:42],1:600*22);
tt = mod(t, 60);
tt = find(tt(2:end) < tt(1:end-1))+1;
stacked_traces(spks, 3, {'m'})
set(gca, 'XTick', tt, 'XTickLabel', round(t(tt)))
set(gcf,'units','centimeters','position',[16,4,12,8]);
set(gca,'FontSize',10,'FontName','Arial');
% export_fig 'D:\F32 FIGS\acc_traces' -eps
%%

% behavior
numsess = 4;
ns = NaN(numsess, 3);
mta = NaN(numsess, 3);

% f = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\';
% fs = {'2022_09_08_H16_09_56_HAB_@placecells.mat',...
% '2022_09_08_H16_52_24_TR0_@placecells.mat',...
% '2022_09_12_H15_51_50_TR3_@placecells.mat',...
% '2022_09_13_H17_52_09_TR6_@placecells.mat'};
f = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\';
fs = {'2023_01_05_H12_52_30_HAB0_@placecells.mat',...
'2023_01_08_H14_16_32_TR0_@placecells.mat',...
'2023_01_09_H14_32_47_TR3_@placecells.mat',...
'2023_01_13_H18_48_32_TR6_@placecells.mat'};
ip = load([f 'ipos_files\Acc20832_TR3_int_0.25sec.mat']);
fname = sprintf('%s%s', f, fs{3});
temp = load(fname, 'ms');
% spks = temp.ms.neuron.C + temp.ms.neuron.YrA/3;
spks = ip.ms_room_temp.spks_bin;%   temp.ms.neuron.S_matw;
ipos = ip.ipos_ra;%   temp.ms.neuron.S_matw;
t = temp.ms.timestamps./1000;% + ms_h.ms.neuron.YrA/4;
spks = normalize_rows(spks);
t = t(1:60*round(22*.25));
spks = spks(:,1:60*round(22*.25));
% ipos = normalize_rows(ipos);
ipos = ipos(:,1:60*round(22*.25));
ipos_mean = normalize_rows(nanmean(ipos,1))-.5; 
ipos_mean = ipos_mean(:,1:60*round(22*.25));
ipos_mean_scale = ipos_mean.*size(ipos,1)/2 + size(ipos,1)/2; 
tt = mod(t, 5);
tt = find(tt(2:end) < tt(1:end-1))+1;
% imagesc(spks)

figure(108);  clf
im1 = (temp.ms.neuron.maxFrame-temp.ms.neuron.minFrame)./75;
[conts, conts_outline] = contour_outlines(temp.ms.neuron.fullA, .35, .65, temp.ms.neuron.dims);
conts_outline = conts_outline/3;
conts = conts./10;
im1 = cat(3, im1+conts_outline+conts, im1+conts_outline, im1+conts_outline+conts);
% im1 = cat(3, im1+conts_outline, im1, im1+conts_outline);
image(im1); axis image off
% imshow(temp.ms.neuron.maxFrame./255)
% title(fname);
set(gcf,'units','centimeters','position',[8,8,16,16]);
set(gca,'FontSize',10,'FontName','Arial');
% figname = sprintf('%s_%d', 'D:\F32 FIGS\FIG5_1_Acc20832_contours', i);
figname = sprintf('%s_%d', 'D:\F32 FIGS\FIG5_1_Acc20832_dffonly', i);
export_fig(figname, '-tiff') % -eps


figure(108);  clf
title(fname);
% stacked_traces(spks, 10, {'k'})
% subplot(1,3,1)
imagesc(1-spks, [0 1]); 
colormap gray
colorbar
set(gca, 'XTick', tt, 'XTickLabel', round(t(tt)), 'YTick', [1 size(spks,1)], 'YDir', 'normal')
set(gcf,'units','centimeters','position',[8,4,8,8]);
set(gca,'FontSize',10,'FontName','Arial');
figname = sprintf('%s_%d', 'D:\F32 FIGS\FIG5_2_Acc20832_binspk', i);
export_fig(figname, '-eps') % -eps
% subplot(1,3,2)
figure(108);  clf
imagesc(ipos, [-.5 .5]); hold on 
colormap redblue
colorbar
set(gca, 'XTick', tt, 'XTickLabel', round(t(tt)), 'YTick', [1 size(ipos,1)], 'YDir', 'normal')
yyaxis('right')
plot(ipos_mean*0, 'k-'); 
plot(ipos_mean, 'k-'); 
set(gca, 'YTick', [-.5 .5], 'YTickLabel', [], 'YDir', 'normal')
ylim([ -1 1])
set(gcf,'units','centimeters','position',[8,4,8,8]);
set(gca,'FontSize',10,'FontName','Arial');
figname = sprintf('%s_%d', 'D:\F32 FIGS\FIG5_3_Acc20832_ipos', i);
export_fig(figname, '-eps') % -eps



figure(104); clf
for i = 1:numsess
    fname = sprintf('%s%s', f, fs{i});
    temp = load(fname, 'ms');
    size(temp.ms.neuron.C)
    ns(i,1) = length(temp.ms.room.shockTimes);
    st = temp.ms.room.shockTimes;
    sind = interp1(temp.ms.timestamps, 1:length(temp.ms.timestamps), st, 'nearest') ;
    mta(i,1) = max(abs(diff(temp.ms.room.shockTimes)))./1000;
    x = temp.ms.room.x;
    y = temp.ms.room.y;
    xn = NaN*x;
    yn = NaN*y;
    [a,r] = cart2pol(x,y);
    inzone = a>(pi/2-pi/6) & a<(pi/2+pi/6);
    xn(inzone) = x(inzone);
    yn(inzone) = y(inzone);
    figure(104);clf
%     subplot(2, 1, 1)
    hold on;
    plot(x, y, 'k')
    scatter(x(sind), yn(sind), 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
%     plot(xn, yn, 'r')
    axis square
%     plot(temp.ms.room.x, temp.ms.room.y, 'k')
    text(30, 40, sprintf('%d shk', ns(i)), 'FontSize', 6)
    text(30, 30, sprintf('%3.1f mta', mta(i)), 'FontSize', 6)
    axis([-45 45 -45 45])
    set(gca, 'XTick', [-40 0 40], 'YTick', [-40 0 40], 'Color', 'none')
set(gcf,'units','centimeters','position',[8,8,4,4]);
set(gca,'FontSize',10,'FontName','Arial');
figname = sprintf('%s_%d', 'D:\F32 FIGS\path_plots_ex', i);
% export_fig(figname, '-eps') % -eps

%     subplot(2, 1, 2)
% %     clf
% %     imagesc(temp.ms.room.vmap, [0 20]); colorbar
% %     colormap viridis
% %     set(gca, 'YDir', 'normal')
% %     axis image off
% % %     subplot(2, numsess, i)
% % set(gcf,'units','centimeters','position',[8,8,6,4]);
% % set(gca,'FontSize',10,'FontName','Arial');
% % figname = sprintf('%s_%d', 'D:\F32 FIGS\path_plots_ex_occ', i);
% export_fig(figname, '-tiff') % -eps

end
%%
% colorbar
% % export_fig 'D:\F32 FIGS\path_plots_ex_colorbar' -eps

f = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\';
fs = {'2023_01_05_H14_53_39_HAB0_@placecells.mat',...
'2023_01_08_H15_29_11_TR0_@placecells.mat',...
'2023_01_09_H15_40_18_TR3_@placecells.mat',...
'2023_01_13_H14_38_46_TR6_@placecells.mat'};
for i = 1:numsess
    fname = sprintf('%s%s', f, fs{i});
    temp = load(fname, 'ms');
    ns(i,2) = length(temp.ms.room.shockTimes);
    mta(i,2) = max(abs(diff(temp.ms.room.shockTimes)))./1000;
end

f = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\';
fs = {'2022_09_08_H16_09_56_HAB_@placecells.mat',...
'2022_09_08_H16_52_24_TR0_@placecells.mat',...
'2022_09_12_H15_51_50_TR3_@placecells.mat',...
'2022_09_13_H17_52_09_TR6_@placecells.mat'};
for i = 1:numsess
    fname = sprintf('%s%s', f, fs{i});
    temp = load(fname, 'ms');
    ns(i,3) = length(temp.ms.room.shockTimes);
    mta(i,3) = max(abs(diff(temp.ms.room.shockTimes)))./1000;
end

% figure(105)
% plot(mta, 'k.-')
% figure(105)
% plot(ns, 'k.-')























