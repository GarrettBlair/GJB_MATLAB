%get manifold embedding
% Computes Isomap embedding using the algorithm of Tenenbaum, de Silva, and Langford (2000). 
% clear; close all; 
clc;
tic
animateplot = 0; 
compDist = 0;
record = 0;
% sizpost = size(frxPost,1);
% idx = randperm(sizpost); % create index to randomize order of units during light
%     shuffrxPost = frxPost;
%     shuffrxPost(idx,:) = frxPost(:,:);
global X %globally defined data
X = frx;

options.dims = 3; %1:3 %range of dimensions to calculate typ 2:20 for dim calculations, 3 for plot
options.landmarks = 1:size(X,2); %recommended value
nNeighbours = 7;

[Y, R, E] = IsomapII('dfun', 'k', nNeighbours, options); 

%% Compute variances

%isomap
vexpISM = 1-R;
% figure
% plot(1:length(variance),variance,'-o')
%PCA
transX = X';
[coeff,~,latent,~,explained,mu] = pca(transX);
vexpPCA = cumsum(explained');
%% Compute means
coordinates = Y.coords{1};
RemovedSamples = [length(X)-length(coordinates)]
indexk = find(Y.index==k)-1; indexk2 = find(Y.index==k2)-1;% in case there are missing values
Xmean = mean(coordinates(1,:)); Ymean = mean(coordinates(2,:)); Zmean = mean(coordinates(3,:)); %mean coordinates
PreXmean = mean(coordinates(1,1:indexk)); PreYmean = mean(coordinates(2,1:indexk)); PreZmean = mean(coordinates(3,1:indexk)); %mean coordinates
PreCentroid = [PreXmean;PreYmean;PreZmean];% just so I can check maually
LiXmean = mean(coordinates(1,indexk:indexk2)); LiYmean = mean(coordinates(2,indexk:indexk2)); LiZmean = mean(coordinates(3,indexk:indexk2));
LiCentroid = [LiXmean;LiYmean;LiZmean];
PostXmean = mean(coordinates(1,indexk2:end)); PostYmean = mean(coordinates(2,indexk:end)); PostZmean = mean(coordinates(3,indexk:end)); %mean coordinates
PostCentroid = [PostXmean;PostYmean;PostZmean];
%% Build reference distance to centroid
movedlight = coordinates;
movedlight(:,indexk:indexk2) = coordinates(:,indexk:indexk2)-(LiCentroid-PreCentroid);
if compDist == 1
  % total
    % distances = coordinates; %get for every point in spaces its distance to mean coordinate
    % distances(1,:) = (coordinates(1,:) - Xmean).^2; distances(2,:) = (coordinates(2,:) - Ymean).^2;  distances(3,:) = (coordinates(3,:) - Zmean).^2;
    % undostres = [mean(distances(1,:)),mean(distances(2,:)),mean(distances(3,:))];%also just to check
    % ans1 = sum(distances); ans2 = sqrt(ans1); meandist = mean(ans2);stdist = std(ans2);

  %pre: distance from all points in pre to pre centroid
    distpre = coordinates(:,1:k-1);
    distpre(1,:) = (coordinates(1,1:k-1)-PreXmean).^2;
    distpre(2,:) = (coordinates(2,1:k-1)-PreYmean).^2;
    distpre(3,:) = (coordinates(3,1:k-1)-PreZmean).^2;
    distpre2 = sqrt(sum(distpre));
    figure
    histogram(distpre2,'BinWidth',4,'Normalization','Probability','FaceColor','black');
    hold on

  %light: distance from all points in light to pre centroid
    distli2pre = coordinates(:,k:k2-1);
    distli2pre(1,:) = (coordinates(1,k:k2-1)-PreXmean).^2;
    distli2pre(2,:) = (coordinates(2,k:k2-1)-PreYmean).^2;
    distli2pre(3,:) = (coordinates(3,k:k2-1)-PreZmean).^2;
    distli2pres2 = sqrt(sum(distli2pre)); 
    histogram(distli2pres2,'BinWidth',4,'Normalization','Probability','FaceColor','blue');
    meanLitoPre = mean(distli2pres2);

  %post: distance from all points in post to pre centroid
    distpo2pre = coordinates(:,k2:end);
    distpo2pre(1,:) = (coordinates(1,k2:end)-PreXmean).^2;
    distpo2pre(2,:) = (coordinates(2,k2:end)-PreYmean).^2;
    distpo2pre(3,:) = (coordinates(3,k2:end)-PreZmean).^2;
    distpo2pre2 = sqrt(sum(distpo2pre));
    histogram(distpo2pre2,'BinWidth',4,'Normalization','Probability','FaceColor',[0.5 0.5 0.5]);

% %% A) Compute distance Light/Post centroid to Pre centroid distance in
% stdev, before displacement
sddistpre = std(distpre2);
sddistpost = std(distpo2pre2);
MeanPretoPre = mean(distpre2);
lictoprec = sqrt(sum((LiCentroid-PreCentroid).^2));
poctoprec = sqrt(sum((PostCentroid-PreCentroid).^2));
ltpsd = (lictoprec-MeanPretoPre)/sddistpre;
ptpsd = (poctoprec-MeanPretoPre)/sddistpre;
summarydist = [lictoprec,ltpsd,poctoprec,ptpsd,mean(distpre2),sddistpre,sddistpost];
%% translate light and find overlap
movedlight = coordinates;
movedlight(:,indexk:indexk2) = coordinates(:,indexk:indexk2)-(LiCentroid-PreCentroid);
MovedLiCentroid = [mean(movedlight(1,indexk:indexk2)),mean(movedlight(2,indexk:indexk2)),mean(movedlight(3,indexk:indexk2))]';
rmovedlight = int16(movedlight);
trmovedlight = rmovedlight';
if length(frxPre)>length(frxLight)
 [commcoordz, ia, ib] = intersect(movedlight(:,k-(k2-k):k)',movedlight(:,k:k2)','rows');
 comcor = sum(ismember(trmovedlight(k-(k2-k):k,:),trmovedlight(k:k2,:),'rows'))
else
 [commcoordz, ia, ib] = intersect(movedlight(:,1:k)',movedlight(:,k:k+k)','rows');
 comcor = sum(ismember(trmovedlight(1:k,:),trmovedlight(k:k+k,:),'rows'))
end

%% B) Distance Light/Post centroid to Pre centroid distance in stdev
sddistpre = std(distpre2);
sddistpost = std(distpo2pre2);
MeanPretoPre = mean(distpre2);
molictoprec = sqrt(sum((MovedLiCentroid-PreCentroid).^2));
moltpsd = (molictoprec-MeanPretoPre)/sddistpre;
summarydist2 = [molictoprec,moltpsd];
%% Calculate mean distance to Pre (from Light, for every point, thing above is distance between centroids, which makes no sense)
distmoli2pre = movedlight(:,k:k2-1);
distmoli2pre(1,:) = (movedlight(1,k:k2-1)-PreXmean).^2;
distmoli2pre(2,:) = (movedlight(2,k:k2-1)-PreYmean).^2;
distmoli2pre(3,:) = (movedlight(3,k:k2-1)-PreZmean).^2;
distmoli2pres2 = sqrt(sum(distmoli2pre));
MeaDismoLitoCentr = mean(distmoli2pres2);
MDMLtcenSD = (MeaDismoLitoCentr-MeanPretoPre)/sddistpre;
newsummary = [MeaDismoLitoCentr,MDMLtcenSD]
figure
histogram(distpre2,'BinWidth',4,'Normalization','Probability','FaceColor','black');hold on; histogram(distmoli2pres2,'BinWidth',4,'Normalization','Probability','FaceColor','blue');
else
end
%% plot isomap either statically (animateplot = 0) or dynamically (animateplot = 1)
figure; hold on;
% subplot(2,2,2)
dimID = 1; %index of dim to plot from options.dims
% arestar = [5;5;5];
% tcoordinates = coordinates(:,1:k)';
% F = donly;
% T=Mesh2Tetra(tcoordinates,F);
if animateplot == 0
   plot3(coordinates(1,1:indexk),coordinates(2,1:indexk),coordinates(3,1:indexk),'.black');hold on
   plot3(coordinates(1,indexk:indexk2),coordinates(2,indexk:indexk2)',coordinates(3,indexk:indexk2),'.b')
   plot3(coordinates(1,indexk2:end),coordinates(2,indexk2:end)',coordinates(3,indexk2:end),'.','Color',[0.8 0.8 0.8])
   figure
   
%%Bounding box,removing outliers, after displacing light
rmcoordinatespre = rmoutliers(movedlight(:,1:indexk-1)');
rmcoordinatesli = rmoutliers(movedlight(:,indexk:indexk2)');
[rotmat1,cornerpoints1,volume1,surface1,edgelength1] = minboundbox(rmcoordinatespre(:,1),rmcoordinatespre(:,2),rmcoordinatespre(:,3));
[rotmat2,cornerpoints2,volume2,surface2,edgelength2] = minboundbox(rmcoordinatesli(:,1),rmcoordinatesli(:,2),rmcoordinatesli(:,3));
plotminbox(cornerpoints1,'black'); hold on
plotminbox(cornerpoints2,'b');

plot3(movedlight(1,1:indexk),movedlight(2,1:indexk),movedlight(3,1:indexk),'.black');hold on
plot3(movedlight(1,indexk:indexk2),movedlight(2,indexk:indexk2)',movedlight(3,indexk:indexk2),'.b')
plotminbox(cornerpoints1,'black'); hold on
plotminbox(cornerpoints2,'b');
diffvolume = (volume1-volume2);
summvol = [volume1,volume2,diffvolume];

else
   
h = animatedline;
h2 = animatedline;
h3 = animatedline;
c1 = h.Color;
c2 = h2.Color;
c3 = h3.Color;
h.Color = 'black';
h2.Color = 'b';
h3.Color = [0.6510    0.6510    0.6510];


for i= 1:indexk
    addpoints(h,coordinates(1,i),coordinates(2,i),coordinates(3,i));
    drawnow
    view(45+(i/5),30)
    pause(0.0001)
end
hold on
for i= 1:(indexk2-indexk)
    addpoints(h2,coordinates(1,indexk+i),coordinates(2,indexk+i),coordinates(3,indexk+i));
    drawnow
    view(45+(indexk/5)+(i/5),30)
    pause(0.0001)
end
for i= 1:length(coordinates)-indexk2
    addpoints(h3,coordinates(1,indexk2+i),coordinates(2,indexk2+i),coordinates(3,indexk2+i));
    drawnow
    view(45+(indexk2/5)+(i/5),30)
    pause(0.0001)
end


% for i= 1:indexk
%     addpoints(h,movedlight(1,i),movedlight(2,i),movedlight(3,i));
%     drawnow
%     view(45+(i/4),30)
%     pause(0.002)
% end
% hold on
% for i= 1:(indexk2-indexk)
%     addpoints(h2,movedlight(1,indexk+i),movedlight(2,indexk+i),movedlight(3,indexk+i));
%     drawnow
%     view(45+(indexk/4)+(i/4),30)
%     pause(0.002)
% end
% % for i= 1:length(movedlight)-indexk2
% %     addpoints(h3,movedlight(1,indexk2+i),movedlight(2,indexk2+i),movedlight(3,indexk2+i));
% %     drawnow
% %     view(45+(indexk2/4)+(i/4),30)
% %     pause(0.002)
% % end


% Bounding box with original recording
% [rotmat1,cornerpoints1,volume1,surface1,edgelength1] = minboundbox(coordinates(1,1:indexk)',coordinates(2,1:indexk)',coordinates(3,1:indexk)')
% [rotmat2,cornerpoints2,volume2,surface2,edgelength2] = minboundbox(coordinates(1,indexk:indexk2)',coordinates(2,indexk:indexk2)',coordinates(3,indexk:indexk2)');
% plotminbox(cornerpoints1,'black'); hold on
% plotminbox(cornerpoints2,'r');
% figure
% plot3(coordinates(1,1:indexk),coordinates(2,1:indexk),coordinates(3,1:indexk),'.black');hold on
% plot3(coordinates(1,indexk:indexk2),coordinates(2,indexk:indexk2),coordinates(3,indexk:indexk2),'.r'),
% plotminbox(cornerpoints1,'black'); hold on
% plotminbox(cornerpoints2,'r');

% %% Bounding box removing outliers
% rmcoordinatespre = rmoutliers(coordinates(:,1:indexk-1)','quartiles');
% rmcoordinatesli = rmoutliers(coordinates(:,indexk:indexk2)','quartiles');
% [rotmat1,cornerpoints1,volume1,surface1,edgelength1] = minboundbox(rmcoordinatespre(:,1),rmcoordinatespre(:,2),rmcoordinatespre(:,3));
% [rotmat2,cornerpoints2,volume2,surface2,edgelength2] = minboundbox(rmcoordinatesli(:,1),rmcoordinatesli(:,2),rmcoordinatesli(:,3));
% plotminbox(cornerpoints1,'black'); hold on
% plotminbox(cornerpoints2,'b');
% figure
% plot3(coordinates(1,1:indexk),coordinates(2,1:indexk),coordinates(3,1:indexk),'.black');hold on
% plot3(coordinates(1,indexk:indexk2),coordinates(2,indexk:indexk2)',coordinates(3,indexk:indexk2),'.b')
% plotminbox(cornerpoints1,'black'); hold on
% plotminbox(cornerpoints2,'b');

%%Bounding box,removing outliers, after displacing light
rmcoordinatespre = rmoutliers(movedlight(:,1:indexk-1)');
rmcoordinatesli = rmoutliers(movedlight(:,indexk:indexk2)');
[rotmat1,cornerpoints1,volume1,surface1,edgelength1] = minboundbox(rmcoordinatespre(:,1),rmcoordinatespre(:,2),rmcoordinatespre(:,3));
[rotmat2,cornerpoints2,volume2,surface2,edgelength2] = minboundbox(rmcoordinatesli(:,1),rmcoordinatesli(:,2),rmcoordinatesli(:,3));
plotminbox(cornerpoints1,'black'); hold on
plotminbox(cornerpoints2,'b');
figure
plot3(movedlight(1,1:indexk),movedlight(2,1:indexk),movedlight(3,1:indexk),'.black');hold on
plot3(movedlight(1,indexk:indexk2),movedlight(2,indexk:indexk2)',movedlight(3,indexk:indexk2),'.b')
plotminbox(cornerpoints1,'black'); hold on
plotminbox(cornerpoints2,'b');

%%
diffvolume = (volume1-volume2);
summvol = [volume1,volume2,diffvolume];
toc

end

%% if record ==1 then copy, modify and run Code Simon sent you.




