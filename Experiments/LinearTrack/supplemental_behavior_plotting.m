%%
presess = shock_sess_diff==-1;
postsess = shock_sess_diff==1;
A_names = Short_rate.animal_id; %%%%%%%%%%%
is_female = Short_rate.is_female; %%%%%%%%%%%
scopo_first = Short_rate.scopo_shock_first; %%%%%%%%%%%
A_num = length(A_names);

all_sess_done = ~any(isnan([Short_rate.SAL_pre Short_rate.SAL_post Short_rate.SCOP_pre Short_rate.SCOP_post]), 2);
sal_done = ~any(isnan([Short_rate.SAL_pre Short_rate.SAL_post]), 2);
scop_done = ~any(isnan([Short_rate.SCOP_pre Short_rate.SCOP_post]), 2);

e_sr_sal = [Short_rate.SAL_pre Short_rate.SAL_post]; % shortrate_saline(:, presess | postsess);
e_sr_sal = e_sr_sal(sal_done,:);
[p_sr_sal] = ranksum(e_sr_sal(:,1), e_sr_sal(:,2));

e_sp_sal = [Short_proportion.SAL_pre Short_proportion.SAL_post]; % e_sp_sal = shortprop_saline(:, presess | postsess);
e_sp_sal = e_sp_sal(sal_done,:);
[p_sp_sal] = ranksum(e_sp_sal(:,1), e_sp_sal(:,2));

e_sr_scop = [Short_rate.SCOP_pre   Short_rate.SCOP_post]; % e_sr_scop = shortrate_scopo(:, presess | postsess); 
e_sr_scop = e_sr_scop(scop_done,:);
[p_sr_scop] = ranksum(e_sr_scop(:,1), e_sr_scop(:,2));

e_sp_scop = [Short_proportion.SCOP_pre   Short_proportion.SCOP_post]; % e_sp_scop = shortprop_scopo(:, presess | postsess); 
e_sp_scop = e_sp_scop(scop_done,:);
[p_sp_scop] = ranksum(e_sp_scop(:,1), e_sp_scop(:,2));


saline_sr_mean = nanmean(shortrate_saline,1);
saline_sp_mean = nanmean(shortprop_saline,1);
scopo_sr_mean  = nanmean(shortrate_scopo,1);
scopo_sp_mean  = nanmean(shortprop_scopo,1);
figure(963); clf; 
d1 = shock_sess_diff(1)-.5;
d2 = shock_sess_diff(end)+.5;
scopocolor = [1 .5 1];
scopomeancolor = [1 0 1];
controlcolor = [.85 .85 .85];
controlmeancolor = [.45 .45 .45];

n_saline = size(e_sr_sal,1); 
n_scopo = size(e_sr_scop,1); 
subplot(221); hold on
rectangle('Position', [shock_sess_diff(presess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
rectangle('Position', [shock_sess_diff(postsess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
title(sprintf('Saline- ranksum(n=%d) p=%0.3f', n_saline, p_sr_sal))
for j = 1:A_num
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    plot(shock_sess_diff, shortrate_saline(j,:), 'Color', controlcolor, 'LineStyle', linetype, 'LineWidth', 1.5); 
end
plot(shock_sess_diff, saline_sr_mean, 'Color', controlmeancolor, 'LineWidth', 2); 
ylabel('short laps / min')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 6])

subplot(222); hold on
rectangle('Position', [shock_sess_diff(presess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
rectangle('Position', [shock_sess_diff(postsess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
title(sprintf('Saline- ranksum(n=%d) p=%0.3f', n_saline, p_sp_sal))
for j = 1:A_num
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    plot(shock_sess_diff, shortprop_saline(j,:), 'Color', controlcolor, 'LineStyle', linetype, 'LineWidth', 1.5); 
end
plot(shock_sess_diff, saline_sp_mean, 'Color', controlmeancolor, 'LineWidth', 2); 
ylabel('short laps / all laps')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 1])
set(gca, 'YTick', [0:.2:1])

subplot(223); hold on
rectangle('Position', [shock_sess_diff(presess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
rectangle('Position', [shock_sess_diff(postsess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
title(sprintf('Scopo- ranksum(n=%d) p=%0.3f', n_scopo, p_sr_scop))
for j = 1:A_num
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    plot(shock_sess_diff, shortrate_scopo(j,:), 'Color', scopocolor, 'LineStyle', linetype, 'LineWidth', 1.5); 
end
plot(shock_sess_diff, scopo_sr_mean, 'Color', scopomeancolor, 'LineWidth', 2); 
ylabel('short laps / min')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 6])

subplot(224); hold on
rectangle('Position', [shock_sess_diff(presess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
rectangle('Position', [shock_sess_diff(postsess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
title(sprintf('Scopo- ranksum(n=%d) p=%0.3f', n_scopo, p_sp_scop))
for j = 1:A_num
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    plot(shock_sess_diff, shortprop_scopo(j,:), 'Color', scopocolor, 'LineStyle', linetype, 'LineWidth', 1.5); 
end
plot(shock_sess_diff, scopo_sp_mean, 'Color', scopomeancolor, 'LineWidth', 2); 
ylabel('short laps / all laps')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 1])
set(gca, 'YTick', [0:.2:1])
set(gcf, 'Color', 'w')

fname = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\fig_scopo_shock_paired_supp.fig';
savefig(gcf, fname)
fname = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\fig_scopo_shock_paired_supp.tiff';
temp = getframe(gcf);
imwrite(temp.cdata, fname)
%%
figure(964); clf; 
plotx_prepost = shock_sess_diff(presess|postsess);
d1 = plotx_prepost(1)-.5;
d2 = plotx_prepost(end)+.5;
scopocolor = [1 .2 1];
controlcolor = [.6 .6 .6];
plot_line_diff = .1;
marker_size = 50;
subplot(121); hold on
t1 = sprintf('Saline- ranksum(n=%d) p=%0.3f', n_saline, p_sr_sal);
t2 = sprintf('Scopo- ranksum(n=%d) p=%0.3f', n_scopo, p_sr_scop);
title(sprintf('%s\n%s', t1, t2))

for j = 1:size(e_sr_sal,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost-plot_line_diff, e_sr_sal(j,:), marker_size/4, 'o', 'MarkerFaceColor', controlcolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost-plot_line_diff, e_sr_sal(j,:), 'Color', controlcolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sr_sal,1), marker_size, 'o', 'MarkerFaceColor', controlcolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sr_sal,1), 'Color', controlcolor/1.3, 'LineWidth', 2); 
ylabel('short laps / min')
xlabel('sessions to shock')
axis([d1 d2 0 5])

subplot(122); hold on
t1 = sprintf('Saline- ranksum(n=%d) p=%0.3f', n_saline, p_sp_sal);
t2 = sprintf('Scopo- ranksum(n=%d) p=%0.3f', n_scopo, p_sp_scop);
title(sprintf('%s\n%s', t1, t2))
for j = 1:size(e_sp_sal,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost-plot_line_diff, e_sp_sal(j,:), marker_size/4, 'o', 'MarkerFaceColor', controlcolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost-plot_line_diff, e_sp_sal(j,:), 'Color', controlcolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sp_sal,1), marker_size, 'o', 'MarkerFaceColor', controlcolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sp_sal,1), 'Color', controlcolor/1.3, 'LineWidth', 2); 
ylabel('short laps / all laps')
xlabel('sessions to shock')
axis([d1 d2 0 1])
set(gca, 'YTick', [0:.2:1])

subplot(121); hold on
for j = 1:size(e_sr_scop,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost+plot_line_diff, e_sr_scop(j,:), marker_size/4, 'd', 'MarkerFaceColor', scopocolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost+plot_line_diff, e_sr_scop(j,:), 'Color', scopocolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sr_scop,1), marker_size, 'd', 'MarkerFaceColor', scopocolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sr_scop,1), 'Color', scopocolor/1.3, 'LineWidth', 2); 
ylabel('short laps / min')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 5])

subplot(122); hold on
for j = 1:size(e_sp_scop,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost+plot_line_diff, e_sp_scop(j,:), marker_size/4, 'd', 'MarkerFaceColor', scopocolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost+plot_line_diff, e_sp_scop(j,:), 'Color', scopocolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sp_scop,1), marker_size, 'd', 'MarkerFaceColor', scopocolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sp_scop,1), 'Color', scopocolor/1.3, 'LineWidth', 2); 
ylabel('short laps / all laps')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 1])
set(gca, 'YTick', [0:.2:1])

set(gcf, 'Color', 'w')

fname = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\fig_pre_post_paired_supp.fig';
savefig(gcf, fname)
fname = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\fig_pre_post_paired_supp.tiff';
temp = getframe(gcf);
imwrite(temp.cdata, fname)
