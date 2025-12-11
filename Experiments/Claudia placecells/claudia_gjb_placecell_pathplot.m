clear
% close all
ddir = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Claudia placecells\';
load([ddir 'H1_1206_gb.mat'])
FR_matrix = FR_matrix(116:215,:);
indeces = indeces-115;
time_sec = (0:length(X)-1)*(1/30);
diff_t = ones(length(X),1)*(1/30);
%%
% figure(1336); clf;
[Th, R] = cart2pol(X,Y);
T = rad2deg(Th);
T = mod(T+180-270, 360);
% dt = (diff(T));
% dt = speed_calc_gb(X,Y);
% figure;
% plot(dt); hold on;
% kern = gausswin(61); kern = kern./sum(kern(:));
% dt = conv(dt, kern, 'same');
% plot(dt);

% plot(X, Y, '-', 'Color', [.7 .7 .7])
X1 = T<30;
X2 = T>330;
% X1 = X<0 & Y>.3;
% X2 = X>0 & Y>.3;
x1_s = find(X1(1:end-1)==1 & X1(2:end)==0);
x2_s = find(X2(1:end-1)==1 & X2(2:end)==0);
x1_e = find(X1(2:end)==1   & X1(1:end-1)==0);
x2_e = find(X2(2:end)==1   & X2(1:end-1)==0);

x1_start_stop = [];
valid = zeros(size(FR_matrix,2), 1);
for i = 1:length(x1_s)
    idx = x1_s(i);
    idx1 = x1_s( find(x1_s>idx, 1) );
    idx2 = x2_e( find(x2_e>idx, 1) );
    if ~isempty(idx2) & idx2<idx1
        x1_start_stop = cat(1, x1_start_stop, [idx, idx2]);
        valid(idx:idx2) = 1;
    end
end
x2_start_stop = [];
for i = 1:length(x2_s)
    idx = x2_s(i);
    idx1 = x2_s( find(x2_s>idx, 1) );
    idx2 = x1_e( find(x1_e>idx, 1) );
    if ~isempty(idx2) & idx2<idx1
        x2_start_stop = cat(1, x2_start_stop, [idx, idx2]);
        valid(idx:idx2) = 2;
    end
end
% hold on
% plot(X(valid==1), Y(valid==1), '-', 'Color', [1 .7 .7])
% 
% plot(X(valid==2), Y(valid==2), '-', 'Color', [.7 .7 1])
%%
figure(1337); clf;
set(gcf, 'Color', 'w')
dim1 = false;
FR_matrix2 = FR_matrix(:, valid>0);
dt2 = diff_t(valid>0);
x = X(valid>0);
y = Y(valid>0);
th = Th(valid>0);
t = T(valid>0);
r = R(valid>0);
kern = gausswin(11); kern = kern./sum(kern(:));
% r = conv(r, kern, 'same');
r = .44+(rand(1,length(r))-.5)/10;
% r = r+(rand(1,length(r))-.5)/10;
[x,y] = pol2cart(th,r);

%
% inds = [56 79 70 57 36 53 93 80];
inds = [56 70 57 53 93];


indeces = 1:1:size(FR_matrix,1);
indeces = inds;
ni = length(indeces);
ni2 = ceil(sqrt(ni));

clrs = viridis(ni);
% clrs = jet(ni);
% clrs(1:find(clrs(:,1)==0.5),1) = 0.5;
clrs = clrs(randperm(ni),:);
% clrs = [clrs(1:2:end,:); clrs(2:2:end,:)];
clrs = light_colormap(clrs, 5);
subplot(1,2,1); hold on
hold on
if dim1
    plot(1:length(th), th, '-', 'Color', [.7 .7 .7])
else
%     plot(x, y, '-', 'Color', [.7 .7 .7])
end
plot(X, Y, '-', 'Color', [.7 .7 .7])
for i = 0:ni
    %%
    if i > 0 
    spk = FR_matrix2(indeces(i),:);
    spk_i = spk>0;
    t2 = dt2;
    t2(spk==0) = 0;
    subplot(1,2,1); hold on
    if dim1
        scatter(find(spk_i), th(spk_i), 30, 'o', 'MarkerFaceColor', clrs(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .8)
    else
        %     scatter(x(spk_i), y(spk_i), 30, 'o', 'MarkerFaceColor', clrs(i,:), 'MarkerEdgeColor', clrs(i,:)*0, 'MarkerFaceAlpha', .8)
        scatter(x(spk_i), y(spk_i), 30, 'o', 'MarkerFaceColor', clrs(i,:), 'MarkerEdgeColor', clrs(i,:)*0, 'MarkerFaceAlpha', 1)
    end
    axis square
    
    subplot(1,2,2); hold on
    bins = [0:15:360];
    s  = make_Ca_map_1D(t, bins, spk);
%     tt = make_Ca_map_1D(t, bins, t2);
    s = s/sum(s); % ./tt;
    
%     s(tt==0) = 0;
    patch([ 0 bins(1:end-1)+2.5 360], [-1 s -1], [-1 (s*0)-1 -1], 'FaceColor', clrs(i,:), 'EdgeColor', clrs(i,:)*.7, 'FaceAlpha', .5, 'LineWidth', 2)
    end
    subplot(1,2,1); hold on
    axis([-1 1 -1 1]*.6)
    set(gca, 'XTick', [], 'YTick', [])
    xlabel('X')
    ylabel('Y')
    subplot(1,2,2); hold on
    axis([0 360 0 .7])
    set(gca, 'XTick', [0:90:360], 'YTick', [0:.2:.8])
    ylabel('Area norm rate')
    xlabel('Angular position (deg)')
    drawnow()
    fname = sprintf('%spathplot_cell%d.png', ddir, i);
    saveas(gcf, fname)
end

