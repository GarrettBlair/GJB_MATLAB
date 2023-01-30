function gb_ziv_topoclustering_plotting(eig_struct, topo_out, spks, position_active_frames)

%%
NumOfBins = max(position_active_frames);
cmap_position = jet(NumOfBins);

TopoClustering = topo_out;
TopoClusterVer=1;
Indexing_Ver=3;
ChosenType=3; %cluster of running
v2 = eig_struct{3};

if ~isempty(topo_out)
NumOfClusters=max(TopoClustering.global.ind)-1;
%cmap_clusters=lines(NumOfClusters+1);
cmap_clusters=jet(NumOfClusters+1);

switch Indexing_Ver
    case 0
        cluster_ind=TopoClustering.global.ind;
    case 1
        cluster_ind=TopoClustering.global.ind_augmented;
    case 2
        cluster_ind=TopoClustering.global.ind_augmented2;
    case 3
        cluster_ind=TopoClustering.global.ind_augmented3;
end
end


%% plot eigs
figure,
plot3(v2(:,2),v2(:,3),v2(:,4),'k.','markersize',5)
hold on
view([20,20])
xlabel('Comp. 1')
ylabel('Comp. 2')
zlabel('Comp. 3')
% Panel1A_lim=0.025;
% xlim([-Panel1A_lim Panel1A_lim])
% ylim([-Panel1A_lim Panel1A_lim])
% zlim([-Panel1A_lim Panel1A_lim])
set(gca,'XTick',[])
set(gca,'XTickLabel',[]);
set(gca,'YTick',[])
set(gca,'YTickLabel',[]);
set(gca,'ZTick',[])
set(gca,'ZTickLabel',[]);
axis square
%%
bg_color = [.2 .2 .2];
write_vid = false;

cm = cmap_position(position_active_frames,:);
% cm = cmap_clusters(cluster_ind,:);
cm_light = (cm/3)+2/3;
eig_lag = 30;
im_lag = 100;
markerfade = round(linspace(3,100, eig_lag+1));
if write_vid
    filename = 'F:\ScaryData2\Hipp12to18\Hipp8_lin7_nonsharedsegs.tiff';
    f = getframe(gcf);
    %     imwrite(f.cdata, filename, 'Compression', 'none', 'WriteMode', 'overwrite');
end
figure(678); clf; set(gcf, 'Position', [460 551 775 504], 'Color', bg_color);
colormap viridis
hold on
Panel1A_lim= max(max([v2(:,2),v2(:,3),v2(:,4)], [], 1));
for ii = im_lag+1:3:length(v2(:,2))-im_lag
    %%
    subplot_tight(2,2,1); cla; hold on
    scatter(1:length(position_active_frames), position_active_frames, 25, cm_light,'fill')
    plot(position_active_frames, 'k-')
    scatter(ii, position_active_frames(ii), 100, cm(ii,:),'fill')
    
    axis tight off
    subplot_tight(2,2,3); cla; hold on
    imagesc(spks(:, ii-im_lag:ii+im_lag));
    hold on; plot([im_lag+1 im_lag+1], [0 size(spks,1)+1], 'r-', 'LineWidth', 3);
    axis tight off
    
    subplot_tight(2,2,[2 4], [.05, .1]); cla; hold on
    scatter3(v2(:,2),v2(:,3),v2(:,4),3,cm_light, 'MarkerFaceAlpha', .05)
%     scatter3(v2(ii-eig_lag:ii,2),v2(ii-eig_lag:ii,3),v2(ii-eig_lag:ii, 4), markerfade, cmap_position(position_active_frames(ii-eig_lag:ii),:),'fill')
    scatter3(v2(ii-eig_lag:ii,2),v2(ii-eig_lag:ii,3),v2(ii-eig_lag:ii, 4), markerfade, cm(ii-eig_lag:ii,:),'fill')
    % plot3(v2(ii-30:ii,2),v2(ii-30:ii,3),v2(ii-30:ii,4),'r-','markersize',5)
    view([mod(round(ii/4),360) , 20])
    xlabel('Comp. 1', 'Color', 'w')
    ylabel('Comp. 2', 'Color', 'w')
    zlabel('Comp. 3', 'Color', 'w')
    %     Panel1A_lim=0.05;
    xlim([-Panel1A_lim Panel1A_lim])
    ylim([-Panel1A_lim Panel1A_lim])
    zlim([-Panel1A_lim Panel1A_lim])
    set(gca,'XTick',[])
    set(gca,'XTickLabel',[]);
    set(gca,'YTick',[])
    set(gca,'YTickLabel',[]);
    set(gca,'ZTick',[])
    set(gca,'ZTickLabel',[]);
    axis square
    set(gca, 'Color', bg_color)
    drawnow
    if write_vid
        f = getframe(gcf);
        imwrite(f.cdata, filename, 'Compression', 'none', 'WriteMode', 'append');
    end
end
%%
figure; hold on
cm = jet(10);
x1 = 1:30;
x2 = 60:-1:31;
dist_vect = NaN(length(x1), 10);
for i = 1:length(x1)
    inds1 = position_active_frames == x1(i);
    inds2 = position_active_frames == x2(i);
    for j = 2:10
%         dists = pdist2(v2(inds1,j), v2(inds2,j), 'mahalanobis');
        dists = abs(median(v2(inds1,j)) - median(v2(inds2,j)));
        dist_vect(i, j) = nanmedian(dists(:));
%         scatter(i*ones(sum(inds1),1), v2(inds1,j), 5, [.4 .4 .4], 'fill')
%         scatter(i, nanmean(v2(inds1,j)), 25, cm(j,:), 'fill')
%         scatter(i, nanmean(v2(inds2,j)), 25, cm(j,:), 'fill')
    end
end

plot(dist_vect)
plot(nanmean(dist_vect,2), 'k', 'LineWidth', 3)
xlabel('Bin Location    start-->finish')
ylabel('Median abs difference between bins')
title('Eigen vector difference along track, 9 eigs')
% x1 = position_active_frames(ismember(position_active_frames, 1:30));
% x2 = position_active_frames(ismember(position_active_frames, 31:60));


%%
%
figure, scatter3(v2(:,2),v2(:,3),v2(:,4),5,cmap_position(position_active_frames,:),'fill')
title('Eigs with position label colored')
view([20,20])
xlabel('Comp. 1')
ylabel('Comp. 2')
zlabel('Comp. 3')
% Panel1A_lim=0.025;
% xlim([-Panel1A_lim Panel1A_lim])
% ylim([-Panel1A_lim Panel1A_lim])
% zlim([-Panel1A_lim Panel1A_lim])
set(gca,'XTick',[])
set(gca,'XTickLabel',[]);
set(gca,'YTick',[])
set(gca,'YTickLabel',[]);
set(gca,'ZTick',[])
set(gca,'ZTickLabel',[]);
axis square
drawnow

%%
switch Indexing_Ver
    case 0
        cluster_ind=TopoClustering.global.ind;
    case 1
        cluster_ind=TopoClustering.global.ind_augmented;
    case 2
        cluster_ind=TopoClustering.global.ind_augmented2;
    case 3
        cluster_ind=TopoClustering.global.ind_augmented3;
end


figure
scatter3(v2(:,2),v2(:,3),v2(:,4),5,cmap_clusters(cluster_ind,:),'fill')
view([20,20])
xlabel('Comp. 1')
ylabel('Comp. 2')
zlabel('Comp. 3')
% Panel1A_lim=0.025;
% xlim([-Panel1A_lim Panel1A_lim])
% ylim([-Panel1A_lim Panel1A_lim])
% zlim([-Panel1A_lim Panel1A_lim])
set(gca,'XTick',[])
set(gca,'XTickLabel',[]);
set(gca,'YTick',[])
set(gca,'YTickLabel',[]);
set(gca,'ZTick',[])
set(gca,'ZTickLabel',[]);
axis square
drawnow

%%
figure
NumOfRows=ceil(sqrt(NumOfClusters));
NumOfColumns=ceil(NumOfClusters/NumOfRows);
for run_cluster=1:NumOfClusters
    subplot(NumOfRows,NumOfColumns,run_cluster)
    scatter3(v2(cluster_ind==run_cluster,2),v2(cluster_ind==run_cluster,3),v2(cluster_ind==run_cluster,4),...
        5,cmap_position(position_active_frames(cluster_ind==run_cluster),:),'fill')
    hold on
    scatter3(v2(cluster_ind~=run_cluster,2),v2(cluster_ind~=run_cluster,3),v2(cluster_ind~=run_cluster,4),...
        1,'fill')
    view([20,20])
    xlabel('Comp. 1')
    ylabel('Comp. 2')
    zlabel('Comp. 3')
    Panel1A_lim=0.025;
    xlim([-Panel1A_lim Panel1A_lim])
    ylim([-Panel1A_lim Panel1A_lim])
    zlim([-Panel1A_lim Panel1A_lim])
    set(gca,'XTick',[])
    set(gca,'XTickLabel',[]);
    set(gca,'YTick',[])
    set(gca,'YTickLabel',[]);
    set(gca,'ZTick',[])
    set(gca,'ZTickLabel',[]);
    axis square
end
drawnow

