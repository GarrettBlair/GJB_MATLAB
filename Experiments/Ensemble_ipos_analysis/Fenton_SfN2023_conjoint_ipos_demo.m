clear
% load('D:\Sample Data\ensemble_prob_linear\standard\256_0.mat');
% load("D:\Sample Data\ensemble_prob_linear\noise_mod\500_0.mat");
% load("D:\Sample Data\ensemble_prob_openfield\standard\256_3.mat");
load("D:\Sample Data\ensemble_prob_linear\standard\256_1.mat")
figure(29); clf; 
p = squeeze(nansum(pf1,1));
[~, pf_ord] = max(p);
[~, pf_ord] = sort(pf_ord);
SAVE_FIGS = false; 
imagesc(bin_spks(pf_ord,200:400)); colormap viridis
hold on; 
nx1 = normalize_matrix(x1)*size(bin_spks,1);
nt = t(200:400); nt = nt-nt(1);
set(gca, 'XTick', 1:60:length(nt), 'XTickLabel', nt(1:60:end))
xlabel('Time (sec)')
ylabel('Unit #')
yyaxis('right') 
plot(-1*x1(200:400), 'r', 'LineWidth', 1)
set(gca, 'YColor', 'r', 'YTick', [-10 0 10])
ylabel('Position')
%
iso = normalize_rows(isoMap') - .5;
sig =  switch_binned;
xx1 = discretize(x1, xbins);
cm = plasma(length(xbins));
ind = find(xx1==8);
idx = (xx1==8);
indind=28;

ensem_dists = squareform(pdist(bin_spks', 'cosine'));
ensem_dists(find(eye(size(bin_spks,2))))=NaN;
% ensem_dists = squareform(pdist(iso', 'cosine'));
% ensem_dists(find(eye(size(iso,2))))=NaN;
ensem_dists_pos = ensem_dists(ind, ind);

this_dists = ensem_dists_pos(indind,:);
other_dists = ensem_dists_pos(:,indind);
all_dists  = ensem_dists(ind(indind),:);

p_x  = sum(idx)/length(idx);
ps  = sum( nanmedian(all_dists)  >= ensem_dists(:) )  / length(ensem_dists(:));
ps_x = sum( nanmedian(this_dists) >= ensem_dists_pos(:) ) / length(ensem_dists_pos(:));
p_x*ps_x*log2(2^ps_x / 2^ps);

figure(3012); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
nc=1; nr=1;
subplot_tight(nr,nc,1, [.075 .075]); hold on
hold on
% scatter3(iso(1,sig), iso(2,sig), iso(3,sig), 30, '.', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', .2)
% scatter3(iso(1,~sig), iso(2,~sig), iso(3,~sig), 30, '.', 'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', .2)
scatter(iso(1,sig), iso(2,sig), 30, '.', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', .2)
scatter(iso(1,~sig), iso(2,~sig), 30, '.', 'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', .2)
axis([-.6 .6 -.6 .6])
% set(gca, 'View', [-12, 50], 'Box', 'on', 'YTickLabel', [], 'XTickLabel', [], 'ZTickLabel', [])
set(gca, 'View', [0, 90], 'Box', 'on', 'YTickLabel', [], 'XTickLabel', [], 'ZTickLabel', [])
axis square
text(-.4, .7, 'Ctx A', 'Color', 'b')
text(.1, .7, 'Ctx B', 'Color', 'r')
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\1-2_iso_proj', 'epsc'); end

figure(3011); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
hold on
scatter(iso(1,:), iso(2,:), 30, cm(xx1,:), '.', 'MarkerFaceAlpha', .2)
set(gca, 'View', [0, 90], 'Box', 'on', 'YTickLabel', [], 'XTickLabel', [], 'ZTickLabel', [])
axis([-.6 .6 -.6 .6])
axis square
subplot_tight(nr*20,nc,1, [.001 .075]); cla; hold on
image(cat(3, cm(end:-1:1,1)', cm(end:-1:1,2)', cm(end:-1:1,3)'))
axis tight off 
text(1, .15, '-10', 'Color', 'k')
text(19, .15, '10', 'Color', 'k')
text(8, .15, 'Position', 'Color', 'k')
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\1-1_iso_proj', 'epsc');end
%
figure(3013); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
hold on
scatter(iso(1,~idx), iso(2,~idx), 30, [.7 .7 .7], '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,idx), iso(2,idx), 90, cm(ind(1),:), '.', 'MarkerFaceAlpha', .2)
% scatter3(iso(1,ind), iso(2,ind), iso(3,ind), 30, 'o', 'MarkerFaceColor', cm(find(ind==1,1),:), 'MarkerFaceAlpha', .3, 'MarkerEdgeColor', 'none')
% scatter3(iso(1,ind), iso(2,ind), iso(3,ind), 30, 'o', 'MarkerFaceColor', cm(find(ind==1,1),:), 'MarkerFaceAlpha', .3, 'MarkerEdgeColor', 'none')
set(gca, 'View', [0, 90], 'Box', 'on', 'YTickLabel', [], 'XTickLabel', [])
axis([-.6 .6 -.6 .6])
axis square
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\1-3_iso_proj', 'epsc'); end

figure(3014); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
hold on
scatter(iso(1,~idx), iso(2,~idx), 30, [.7 .7 .7], '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,idx), iso(2,idx), 90, light_colormap(cm(ind(1),:), 2), '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,ind(indind)), iso(2,ind(indind)), 240, cm(ind(1),:), '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,ind(indind)), iso(2,ind(indind)), 90, 'ro', 'MarkerFaceAlpha', 0, 'MarkerEdgeColor', 'r')
set(gca, 'View', [0, 90], 'Box', 'on', 'YTickLabel', [], 'XTickLabel', [])
axis([-.6 .6 -.6 .6])
axis square
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\1-4_iso_proj', 'epsc'); end

nc=1; nr=1;
figure(3021); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
ind_all = find(xx1>0);
for i = 1:24:length(ind_all)
    ii = ind_all(i);
    plot([ iso(1, (ind(indind))), iso(1,ii)], [ iso(2, (ind(indind))), iso(2,ii)],...
        '-', 'Color', [.9 .3 .3], 'LineWidth', .1)
end
scatter(iso(1,ind_all(1:6:end)), iso(2,ind_all(1:6:end)), 90, [.1 .1 .1], '.', 'MarkerFaceAlpha', .1)
scatter(iso(1,ind(indind)), iso(2,ind(indind)), 240, cm(ind(1),:), '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,ind(indind)), iso(2,ind(indind)), 90, 'ro', 'MarkerFaceAlpha', 0, 'MarkerEdgeColor', 'r')
set(gca, 'View', [0, 90], 'Box', 'on', 'YTickLabel', [], 'XTickLabel', [])
axis([-.6 .6 -.6 .6])
axis square
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\2-1_prob_spk', 'epsc'); end

figure(3022); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
for i = 1:15:length(ind_all)
    ii = ind_all(i);
    sublines = randperm(length(ind_all));
for j = 1:2
    jj = ind_all(sublines(j));
    plot([ iso(1, ii), iso(1, jj)],[ iso(2, ii), iso(2, jj)], '-', 'Color', [.5 .5 .5], 'LineWidth', .1)
end
end
scatter(iso(1,ind_all), iso(2,ind_all), 90, [.1 .1 .1], '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,ind(indind)), iso(2,ind(indind)), 240, cm(ind(1),:), '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,ind(indind)), iso(2,ind(indind)), 90, 'ro', 'MarkerFaceAlpha', 0, 'MarkerEdgeColor', 'r')
set(gca, 'View', [0, 90], 'Box', 'on', 'YTickLabel', [], 'XTickLabel', [])
axis([-.6 .6 -.6 .6])
axis square
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\2-2_prob_spk', 'epsc'); end

figure(3023); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
imagesc(ensem_dists); colormap viridis
rectangle('Position', [0, ind(indind)-1, length(ind_all), 2], 'EdgeColor', [.9 .3 .3], 'LineWidth', 1)
axis image off
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\2-3_prob_spk', 'epsc'); end

figure(3024); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
imagesc(ensem_dists, [.2 1]); colormap viridis
rectangle('Position', [0, 0, length(ind_all), length(ind_all)], 'EdgeColor', 'k', 'LineWidth', 2)
axis image off
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\2-4_prob_spk', 'epsc'); end

figure(3025); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   5   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.15 .15]); hold on
% histogram(ensem_dists(:), [0:.025:1], 'Normalization', 'probability', 'FaceColor', 'k'); 
[hc, bc] = histcounts(ensem_dists(:), [0:.025:1], 'Normalization', 'probability'); 
bc = bc(1:end-1) + abs(diff(bc));
bar(bc, hc, 'BarWidth', 1, 'FaceColor', [.5 .5 .5])
axis([-.1 1.1 0 .2]) 
ylabel('Probability')
yyaxis('right') 
plot([nanmedian(all_dists) nanmedian(all_dists)], [0 .3], '-', 'Color', [.9 .3 .3], 'LineWidth', 2)
set(gca, 'YColor', 'k', 'YTickLabel', [], 'YTick', [])
text(.5, .15, sprintf('p(s) = %01.2f', ps))
axis([-.1 1.1 0 .2])
xlabel('Distance (AU)')
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\2-5_prob_spk', 'epsc'); end

%
figure(3031); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
for i = 1:6:length(ind)
    ii = ind(i);
    plot([ iso(1, (ind(indind))), iso(1,ii)], [ iso(2, (ind(indind))), iso(2,ii)],...
        '-', 'Color', [.9 .3 .3], 'LineWidth', .1)
end
scatter(iso(1,ind(1:6:end)), iso(2,ind(1:6:end)), 90, light_colormap(cm(ind(1),:), 2), '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,ind(indind)), iso(2,ind(indind)), 240, cm(ind(1),:), '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,ind(indind)), iso(2,ind(indind)), 90, 'ro', 'MarkerFaceAlpha', 0, 'MarkerEdgeColor', 'r')
set(gca, 'View', [0, 90], 'Box', 'on', 'YTickLabel', [], 'XTickLabel', [], 'ZTickLabel', [])
axis([-.6 .6 -.6 .6])
axis square
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\3-1_prob_spk_position', 'epsc'); end

figure(3032); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
for i = 1:length(ind)
    ii = ind(i);
    sublines = randperm(length(ind));
for j = 1:2
    jj = ind(sublines(j));
    plot([ iso(1, ii), iso(1, jj)],[ iso(2, ii), iso(2, jj)], '-', 'Color', [.5 .5 .5], 'LineWidth', .1)
end
end
scatter(iso(1,ind), iso(2,ind), 90, light_colormap(cm(ind(1),:), 2), '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,ind(indind)), iso(2,ind(indind)), 240, cm(ind(1),:), '.', 'MarkerFaceAlpha', .2)
scatter(iso(1,ind(indind)), iso(2,ind(indind)),90, 'ro', 'MarkerFaceAlpha', 0, 'MarkerEdgeColor', 'r')
set(gca, 'View', [0, 90], 'Box', 'on', 'YTickLabel', [], 'XTickLabel', [], 'ZTickLabel', [])
axis([-.6 .6 -.6 .6])
axis square
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\3-2_prob_spk_position', 'epsc'); end

figure(3033); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
imagesc(ensem_dists_pos); colormap viridis
rectangle('Position', [0, indind-1, length(ind), 2], 'EdgeColor', 'r', 'LineWidth', 1)
axis image off
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\3-3_prob_spk_position', 'epsc'); end

figure(3034); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   3   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.075 .075]); hold on
imagesc(ensem_dists_pos); colormap viridis
rectangle('Position', [0, 0, length(ind), length(ind)], 'EdgeColor', 'k', 'LineWidth', 2)
axis image off
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\3-4_prob_spk_position', 'epsc'); end

figure(3035); clf;
set(gcf, 'Units', 'inches', 'OuterPosition', [1   1   5   4])%, 'PaperSize', [1,1])
subplot_tight(nr,nc,1, [.15 .15]); hold on
% histogram(ensem_dists_pos(:), [0:.025:1], 'Normalization', 'probability', 'FaceColor', cm(ind(1),:))
[hc, bc] = histcounts(ensem_dists_pos(:), [0:.025:1], 'Normalization', 'probability'); 
bc = bc(1:end-1) + abs(diff(bc));
bar(bc, hc, 'BarWidth', 1, 'FaceColor', cm(ind(1),:))
axis([-.1 1.1 0 .2]) 
ylabel('Probability')
yyaxis('right') 
plot([nanmedian(this_dists) nanmedian(this_dists)], [0 .3], '-', 'Color', [.9 .3 .3], 'LineWidth', 2)
text(.5, .15, sprintf('p(s|x) = %01.2f', ps_x))
axis([-.1 1.1 0 .2]) 
set(gca, 'YColor', 'k', 'YTickLabel', [], 'YTick', [])
xlabel('Distance (AU)')
if SAVE_FIGS==true; saveas(gcf, 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\3-5_prob_spk_position', 'epsc'); end


%%
clear
% load('D:\Sample Data\ensemble_prob_linear\standard\256_0.mat');
% load("D:\Sample Data\ensemble_prob_linear\noise_mod\500_0.mat");
load("D:\Sample Data\ensemble_prob_openfield\standard\256_3.mat");
% load("D:\Sample Data\ensemble_prob_linear\standard\256_1.mat")

ee = abs(ensem_av1) - abs(ensem_av2);
ii = nanmean(abs(ipos1) - abs(ipos2), 1);


p = squeeze(nansum(pf1,1));
[~, pf1_ord] = max(p);
[~, pf1_ord] = sort(pf1_ord);
p = squeeze(nansum(pf2,1));
[~, pf2_ord] = max(p);
[~, pf2_ord] = sort(pf2_ord);




s =  normalize_rows(bin_spks);
sig =  switch_binned;
ipos = abs(ipos1) - abs(ipos2);
iso = normalize_rows(isoMap') - .5;
figure(10); clf;
set(gcf, 'Color', 'w', 'Position', [62   130   843   789])
win = 90;
% fs = 1700-win:1:1850;
fs = 1900-win:1:2350;
last_idx1 = fs(1);
last_idx0 = fs(1);

% v = Fast_Tiff_Write('C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\test.tiff');
if SAVE_FIGS==true
v = VideoWriter('C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Ensemble_ipos_analysis\figs SfN2023\openfield', 'MPEG-4');
% v.Quality = 100;
v.FrameRate = 15;
v.open();
end
for i = 1:length(fs)
    %%
    idx = fs(i);
    idx_win = fs(i)-win:fs(i)+win;
    
    s1_win= s(pf1_ord, idx_win);
    s2_win= s(pf2_ord, idx_win);
    s1_win = 1-cat(3, s1_win, s1_win, s1_win*0);
    s2_win = 1-cat(3, s2_win*0, s2_win, s2_win);
    s_win = cat(1, s2_win, s1_win);
    
    figure(10); clf;
    
    %     isowin = floor(idx-win/3):ceil(idx+win/3);
    isowin = floor(idx-3):ceil(idx);
    
    
    subplot_tight(2,2,1, [.075 .075]); hold on
    if var(y1)==0
        
        if sig(idx)==0
            plot([-10 10], [0 0]+1, '-', 'Color', [.7 .7 .9])
            %     plot(x1(isowin), y1(isowin)+1, '-b')
            text(-5, 0.5, '\color{blue}CTX A')
            scatter(x1(idx), y1(idx)+1, 250, '.', 'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', .2)
            
            plot([-10 10], [0 0]-1, '-k')
            scatter(x2(last_idx0), y2(last_idx0)-1, 100, '.', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .2)
            ar = annotation('arrow', 'Position', [.29 .71 0 .05], 'Color', 'b');
            last_idx1 = idx;
        else
            plot([-10 10], [0 0]+1, '-k')
            scatter(x1(last_idx1), y1(last_idx1)+1, 100, '.', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .2)
            
            text(-5, -.5, '\color{red}CTX B')
            plot([-10 10], [0 0]-1, '-', 'Color', [.9 .7 .7])
            %     plot(x1(isowin), y1(isowin)-1, '-r')
            scatter(x2(idx), y2(idx)-1, 250, '.', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', .2)
            ar = annotation('arrow', 'Position', [.29 .76 0 -.05], 'Color', 'r', 'LineWidth', 2);
            last_idx0 = idx;
        end
        ylim([-3 3])
        axis off
    else
        if sig(idx)==0
%             plot([-10 10], [0 0]+1, '-', 'Color', [.7 .7 .9])
            rectangle('Position', [-10 2 20 20], 'EdgeColor', [.7 .7 .9])
            plot(x1(isowin), y1(isowin)+12, '-', 'Color', [.7 .7 .9])
            text(-5, 0.5, '\color{blue}CTX A')
            scatter(x1(idx), y1(idx)+12, 250, '.', 'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', .2)
            
%             plot(x2, y2-10, '-k')
            rectangle('Position', [-10 -22 20 20], 'EdgeColor', 'k')
            scatter(x2(last_idx0), y2(last_idx0)-12, 100, '.', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .2)
            ar = annotation('arrow', 'Position', [.29 .71 0 .05], 'Color', 'b');
            last_idx1 = idx;
        else
%             plot([-10 10], [0 0]+1, '-k')
%             plot(x1, y1-10, '-k')
            rectangle('Position', [-10 2 20 20], 'EdgeColor', 'k')
            scatter(x1(last_idx1), y1(last_idx1)+10, 100, '.', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .2)
            
            text(-5, -.5, '\color{red}CTX B')
            plot(x2(isowin), y2(isowin)-12, '-', 'Color', [.9 .7 .7])
            rectangle('Position', [-10 -22 20 20], 'EdgeColor', [.9 .7 .7])
%             plot(x2, y2-10, '-', 'Color', [.9 .7 .7])
            scatter(x2(idx), y2(idx)-12, 250, '.', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', .2)
            ar = annotation('arrow', 'Position', [.29 .76 0 -.05], 'Color', 'r', 'LineWidth', 2);
            last_idx0 = idx;
        end
        ylim([-25 25])
        axis off
    end
%     if var(y1)>0
%     set(gca, 'View', [-12+i, 50], 'Box', 'off')
%     else
%     set(gca, 'View', [0 75], 'Box', 'off')
%     end
    subplot_tight(2,2,2, [.075 .075]); hold on
    scatter3(iso(1,sig), iso(2,sig), iso(3,sig), 60, '.', 'MarkerEdgeColor', 'r', 'MarkerEdgeAlpha', .15)
    scatter3(iso(1,~sig), iso(2,~sig), iso(3,~sig), 60, '.', 'MarkerEdgeColor', 'b', 'MarkerEdgeAlpha', .15)
%     scatter3(iso(1,:), iso(2,:), iso(3,:), 3, '.', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .2)
    plot3(iso(1,isowin), iso(2,isowin), iso(3,isowin), 'k-', 'LineWidth', 2)
    scatter3(iso(1,idx), iso(2,idx), iso(3,idx), 80, 'ok', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1, 'MarkerFaceColor', [.2 .2 .2])
    axis([-.6 .6 -.6 .6 -.6 .6])
    set(gca, 'View', [-12+i, 50], 'Box', 'on')
    
    subplot_tight(2,2,3:4, [.1 .1]); hold on
    image(s_win);
    plot([win+1 win+1], [.5 size(s_win,1)+.5], 'k-', 'LineWidth', 1)
    ylabel('Unit #')
    set(gca, 'XTick', win+1, 'XTickLabel', sprintf('Time %3.0f sec',  round(t(idx) - t(fs(1)))),...
        'YTick', [1 128 230 270 256+128 256*2], 'YTickLabel', [0 128 256 0 128 256])
%     text(win-25, -22, 'Time (sec) - ')
    axis tight

    
    subplot_tight(15,1,15, [.01 .1]); hold on
    plot(ensem_av1(idx_win), '-', 'LineWidth', 2, 'Color', [.5 .8 1]);
    plot([win+1 win+1], [min(ensem_av1) 1.2*max(ensem_av1)], 'k-', 'LineWidth', 1)
    text(win+win*.3, max(ensem_av1)*1.2, 'Conjoint Ipos, \color{blue}CTX A')
    axis tight off
    
    drawnow
if SAVE_FIGS==true
    temp = getframe(gcf);
v.writeVideo(temp.cdata)
end
end
if SAVE_FIGS==true
v.close()
end
% end

% [models] = gb_conjoint_ensemble_model_fitting(ensem_av1, thresh_e, edges_of_bins, plotting)