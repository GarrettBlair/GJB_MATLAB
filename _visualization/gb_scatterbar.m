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