function gb_scatterbar(d, x, colors, mean_style)

% d = data matrix, subject by condition
% x = x value to plot on
%%
dm = nanmean(d,1);
ds = nanstd(d,[],1);
[ns, nd] = size(d);
barcolor = shift_colormap(colors, 4);
linecolor = shift_colormap(colors, -4);

cond_jitter = NaN(ns,nd);
for i = 1:nd
    cond_jitter(:,i) = gb_rand_jitter(d(:,i), 24);
end
for i = 1:nd
    if strcmp(mean_style, 'bar')
        bar(x(i), dm(i), 'FaceColor', barcolor, 'EdgeColor', 'k', 'BarWidth', .2, 'LineWidth', 2)
    elseif strcmp(mean_style, 'boxplot')
        boxplot(d(:,i), 'Positions', x(i))
    elseif strcmp(mean_style, 'meanline')
        plot([x(i) x(i)], [dm(i)-ds(i), dm(i)+ds(i)], 'Color', barcolor, 'LineWidth', 10)
    end
end

for i = 1:ns
    plot(x + cond_jitter(i,:), d(i,:), 'Color', linecolor)
end

for i = 1:nd
    xs = x(i) + zeros(ns,1) + cond_jitter(:,i);
    scatter(xs, d(:,i), 'MarkerFaceColor', barcolor, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .8);
end