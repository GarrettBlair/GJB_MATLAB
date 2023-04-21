%% Plotting for cell pair correlation analysis
animal_num_shk = cellmat_Anum(sesslabel>0);
animal_num_scop = cellmat_Anum(sesslabel==0);

shock_pre_corr = train_pre_dist(sesslabel>0);
shock_post_corr = train_post_dist(sesslabel>0);
scop_pre_corr = train_pre_dist(sesslabel==0);
scop_post_corr = train_post_dist(sesslabel==0);
% 0==scopo+shock, 1==shock, 2==2nd shock, 3==barrier, 4==scopo alone
figure(213); clf; hold on; 
w = .15;
b1 = bar(-5, 0, 'FaceColor', [1 0 0], 'EdgeColor', 'k', 'BarWidth', .2, 'LineWidth', 2);
b2 = bar(-5, 0, 'FaceColor', [1 0 1], 'EdgeColor', 'k', 'BarWidth', .2, 'LineWidth', 2);

gb_scatterbar([shock_pre_corr shock_post_corr],  1.5+[-w +w], [1 0 0])
gb_scatterbar([scop_pre_corr scop_post_corr],  2.5+[-w +w], [1 0 1])
% gb_scatterbar(scoppre_corr, 2.4+[-w +w], [1 0 1])
% gb_scatterbar(scopshock_corr, 2.7+[-w +w], [1 0 1])
maxx = .04;%max([train_pre_dist; train_post_dist])*1.2;
minn = -.01;%min([train_pre_dist; train_post_dist])*1.2;
axis([.75 3.25 minn maxx])
set(gca, 'XTick', [1.3 1.7 2.3 2.7], 'XTickLabel', {'PRE-TRAIN' 'PRE-POST' 'PRE-TRAIN' 'PRE-POST' },...
    'XTickLabelRotation', 70, 'YTick', [-.1:.1:.6]);
title('Corr of cell-pair corr matrices')
ylabel('Corr. of cellpair cofiring')
xlabel('Sessions compared')
hl = legend([b1 b2], {'Shk' 'Scopo+Shk'});
% 0==scopo+shock, 1==shock, 2==2nd shock, 3==barrier, 4==scopo alone

% [h, p] = ttest2(scop_pre_corr - scop_post_corr, shock_pre_corr - shock_post_corr)

vars2save = {'cellmat_shockcols' 'sesslabel' 'is_shock' 'cellmat_Anum' 'deltaT' 'corr_pre' 'corr_post'...
    'animal_num_shk' 'animal_num_scop' 'shock_pre_corr' 'shock_post_corr' 'scop_pre_corr' 'scop_post_corr'};

% save([topdir '\cellpair_corrs.mat'], vars2save{:});
% savefig(gcf, [topdir '\corr_of_cellpair_corrs.fig']);



%%
shock_pre_overlap = train_in_pre(sesslabel>0);
shock_post_overlap = post_in_pre(sesslabel>0);
scop_pre_overlap = train_in_pre(sesslabel==0);
scop_post_overlap = post_in_pre(sesslabel==0);
% 0==scopo+shock, 1==shock, 2==2nd shock, 3==barrier, 4==scopo alone
figure(214); clf; hold on; 
w = .15;
b1 = bar(-5, 0, 'FaceColor', [1 0 0], 'EdgeColor', 'k', 'BarWidth', .2, 'LineWidth', 2);
b2 = bar(-5, 0, 'FaceColor', [1 0 1], 'EdgeColor', 'k', 'BarWidth', .2, 'LineWidth', 2);

% gb_scatterbar([shock_pre_overlap shock_post_overlap],  1.5+[-w +w], [1 0 0])
% gb_scatterbar([scop_pre_overlap scop_post_overlap],  2.5+[-w +w], [1 0 1])
gb_scatterbar([shock_pre_overlap-shock_post_overlap],  1.5+[-w +w], [1 0 0])
gb_scatterbar([scop_pre_overlap-scop_post_overlap],  2.5+[-w +w], [1 0 1])
maxx = 1;%max([train_pre_dist; train_post_dist])*1.2;
minn = 0;%min([train_pre_dist; train_post_dist])*1.2;

axis([.75 3.25 minn maxx])
set(gca, 'XTick', [1.3 1.7 2.3 2.7], 'XTickLabel', {'PRE-TRAIN' 'PRE-POST' 'PRE-TRAIN' 'PRE-POST' },...
    'XTickLabelRotation', 70, 'YTick', [-.1:.1:.6]);
title('Proportion of cell pairs significantly corr in each session')
xlabel('Sessions compared')
ylabel('Proportion of cellpairs')
hl = legend([b1 b2], {'Shk' 'Scopo+Shk'});
% 0==scopo+shock, 1==shock, 2==2nd shock, 3==barrier, 4==scopo alone
% [h, p] = ttest2(scop_pre_overlap - scop_post_overlap, shock_pre_overlap - shock_post_overlap)
%%%%%% Hipp12 nice example
% vars2save = {'cellmat_shockcols' 'sesslabel' 'is_shock' 'cellmat_Anum' 'deltaT' 'overlap_pre' 'overlap_post'...
%     'animal_num_shk' 'animal_num_scop' 'shock_pre_overlap' 'shock_post_overlap' 'scop_pre_overlap' 'scop_post_overlap'};
% 
% save([topdir '\intersection_sig_cell_pairs.mat'], vars2save{:});
% savefig(gcf, [topdir '\intersection_sig_cell_pairs.fig']);
%%
%%
shock_pre_overlap_e = pred_pretr_eigs(sesslabel>0);
shock_post_overlap_e = pred_prepost_eigs(sesslabel>0);
scop_pre_overlap_e = pred_pretr_eigs(sesslabel==0);
scop_post_overlap_e = pred_prepost_eigs(sesslabel==0);
shock_pre_overlap = pred_pretr_spks(sesslabel>0);
shock_post_overlap = pred_prepost_spks(sesslabel>0);
scop_pre_overlap = pred_pretr_spks(sesslabel==0);
scop_post_overlap = pred_prepost_spks(sesslabel==0);
% 0==scopo+shock, 1==shock, 2==2nd shock, 3==barrier, 4==scopo alone
figure(214); clf; hold on; 
w = .15;
b1 = bar(-5, 0, 'FaceColor', [1 0 0].*.8, 'EdgeColor', 'k', 'BarWidth', .2, 'LineWidth', 2);
b2 = bar(-5, 0, 'FaceColor', [1 0 1].*.8, 'EdgeColor', 'k', 'BarWidth', .2, 'LineWidth', 2);
b1 = bar(-5, 0, 'FaceColor', [1 0 0], 'EdgeColor', 'k', 'BarWidth', .2, 'LineWidth', 2);
b2 = bar(-5, 0, 'FaceColor', [1 0 1], 'EdgeColor', 'k', 'BarWidth', .2, 'LineWidth', 2);

% gb_scatterbar([shock_pre_overlap shock_post_overlap],  1.5+[-w +w], [1 0 0])
% gb_scatterbar([scop_pre_overlap scop_post_overlap],  2.5+[-w +w], [1 0 1])
gb_scatterbar([shock_pre_overlap, shock_post_overlap],  .75+[-w +w], [1 0 0].*.8)
gb_scatterbar([scop_pre_overlap, scop_post_overlap],  2.5+[-w +w], [1 0 1].*.8)
gb_scatterbar([shock_pre_overlap_e, shock_post_overlap_e],  1.5+[-w +w], [1 0 0])
gb_scatterbar([scop_pre_overlap_e, scop_post_overlap_e],  3.25+[-w +w], [1 0 1])
plot([-1 6], [.5 .5], 'k:')
maxx = 1.2;%max([train_pre_dist; train_post_dist])*1.2;
minn = 0;%min([train_pre_dist; train_post_dist])*1.2;

axis([0 4 minn maxx])
set(gca, 'XTick', [1.3 1.7 2.3 2.7], 'XTickLabel', {'PRE-TRAIN' 'PRE-POST' 'PRE-TRAIN' 'PRE-POST' },...
    'XTickLabelRotation', 70, 'YTick', [.2:.2:1]);
title('Proportion of cell pairs significantly corr in each session')
xlabel('Sessions compared')
ylabel('Decoding of session')
hl = legend([b1 b2], {'Shk' 'Scopo+Shk'});



%%
function gb_scatterbar(d, x, colors)
% d = data matrix, subject by condition
% x = x value to plot on
%%
dm = nanmean(d,1);
[ns, nd] = size(d);
barcolor = shift_colormap(colors, 4);
linecolor = shift_colormap(colors, -4);

cond_jitter = NaN(ns,nd);
for i = 1:nd
    cond_jitter(:,i) = gb_rand_jitter(d(:,i), 24);
end

for i = 1:nd
    bar(x(i), dm(i), 'FaceColor', barcolor, 'EdgeColor', 'k', 'BarWidth', .2, 'LineWidth', 2)
end

for i = 1:ns
    plot(x + cond_jitter(i,:), d(i,:), 'Color', linecolor)
end

for i = 1:nd
    xs = x(i) + zeros(ns,1) + cond_jitter(:,i);
    scatter(xs, d(:,i), 'MarkerFaceColor', barcolor, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .8);
end
end

function cmap = shift_colormap(cmap, kappa)
if kappa>0
cmap = cmap + (1-cmap)/kappa;
elseif kappa<0
cmap = cmap - (cmap)/abs(kappa);
cmap = cmap - min(cmap);
end
end
