clear
% load('D:\Sample Data\ensemble_prob_linear\standard\256_0.mat');
% load("D:\Sample Data\ensemble_prob_linear\noise_mod\500_0.mat");
% load("D:\Sample Data\ensemble_prob_openfield\standard\256_3.mat");
% load("D:\Sample Data\ensemble_prob_linear\standard\256_1.mat")
load("D:\Sample Data\ensemble_prob_linear\noise_mod\750_0.mat")
% load('D:\Sample Data\ensemble_prob_linear\radial_pos\radial_0.mat')

% if ~isvarname('x1') && isname('phi1')
%     x1 = phi1; x2 = phi2; y1 = rho1; y2 = rho2;
% end
ee = abs(ensem_av1) - abs(ensem_av2);
ii = nanmean(abs(ipos1) - abs(ipos2), 1);


p = squeeze(nansum(pf1,1));
[~, pf_ord] = max(p);
[~, pf_ord] = sort(pf_ord);




s =  normalize_rows(bin_spks);
sig =  switch_binned;
ipos = abs(ipos1) - abs(ipos2);
iso = normalize_rows(isoMap') - .5;
figure(10); clf;

fs = 1850:1:2150;
win = 30;
for i = win+1:length(fs)-win
    %%
    idx = fs(i);
    idx_win = fs(i-win:i+win);
    
    s_win= s(pf_ord, idx_win);
    
    figure(10); clf;
    
%     isowin = floor(idx-win/3):ceil(idx+win/3);
    isowin = floor(idx-2):ceil(idx+2);
    
    
    subplot_tight(2,2,1, [.075 .075]); hold on
    if sig(idx)==0
    scatter3(x1, y1, x1*0+1, 30, '.', 'MarkerEdgeColor', [.7 .7 .9], 'MarkerFaceAlpha', .2)
    plot3(x1(isowin), y1(isowin), x1(isowin)*0+1, 'b')
    scatter3(x1(idx), y1(idx), x1(idx)*0+1, 100, '.', 'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', .2)

    scatter3(x2, y2, x2*0+3, 30, '.', 'MarkerEdgeColor', [.7 .7 .7], 'MarkerFaceAlpha', .2)
    plot3(x2(isowin), y2(isowin), x2(isowin)*0+3, 'k')
    scatter3(x2(idx), y2(idx), x2(idx)*0+3, 100, '.', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .2)
    else
    scatter3(x1, y1, x1*0+1, 30, '.', 'MarkerEdgeColor', [.7 .7 .7], 'MarkerFaceAlpha', .2)
    plot3(x1(isowin), y1(isowin), x1(isowin)*0+1, 'k')
    scatter3(x1(idx), y1(idx), x1(idx)*0+1, 100, '.', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .2)

    scatter3(x2, y2, x2*0+3, 30, '.', 'MarkerEdgeColor', [.9 .7 .7], 'MarkerFaceAlpha', .2)
    plot3(x2(isowin), y2(isowin), x2(isowin)*0+3, 'r')
    scatter3(x2(idx), y2(idx), x2(idx)*0+3, 100, '.', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', .2)
    end
    axis off
    if var(y1)>0
    set(gca, 'View', [-12+i, 50], 'Box', 'off')
    else
    set(gca, 'View', [0 75], 'Box', 'off')
    end
    subplot_tight(2,2,2, [.075 .075]); hold on
    scatter3(iso(1,sig), iso(2,sig), iso(3,sig), 30, '.', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', .2)
    scatter3(iso(1,~sig), iso(2,~sig), iso(3,~sig), 30, '.', 'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', .2)
%     scatter3(iso(1,:), iso(2,:), iso(3,:), 3, '.', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .2)
    plot3(iso(1,isowin), iso(2,isowin), iso(3,isowin), 'r-')
    scatter3(iso(1,idx), iso(2,idx), iso(3,idx), 'ro', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', .7, 'MarkerFaceColor', [.8 .2 .2])
    axis([-.6 .6 -.6 .6 -.6 .6])
    set(gca, 'View', [-12+i, 50], 'Box', 'on')
    
    subplot_tight(2,2,3:4, [.1 .1]); hold on
    imagesc(s_win, [0 1]);
    plot([win+1 win+1], [.5 size(s,1)+.5], 'r-', 'LineWidth', 2)
    axis tight
    colormap plasma
    
    drawnow
    pause(1/30)
    
end