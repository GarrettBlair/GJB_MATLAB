clear;  clc
close all;
addpath('IsoMap/');

pathAll = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/CircleSquare/analysisFiles/';

noHMC = 0;

animalList = {'M19', 'M20', 'M29','M35','M34','M39'};

col = [0,0,0;0,1,0;1,0.5,0];
col2 = [0,0,0;0.2,1,0;1,0.5,0.2];

% create a smoothing gaussian
lenGauss = 50;
sG = 3; % careful with this
convGauss = gaussmf(0:2:lenGauss,[sG lenGauss/2]);

envIdx = [-1,0,1];

iAni = 3

%load data
load([pathAll, '1000ms_CircleSquare_', animalList{iAni}, '.mat'])

iDay = 8
%     b = cat(1,a{i,1,:});
%     b2 = cat(1,a{i,2,:});
%     c = cat(2,b,b2);

%get data
s = spike{iDay}';
sessIndic = sess{iDay};

if noHMC
    s = s(:,sess{iDay}>=0);
    sessIndic = sessIndic(:,sess{iDay}>=0);
end

nCells = size(s,1);

% create smoothed version of spiking
sConv = nan(size(s));
for iCell = 1:nCells
    sConv(iCell,:) = (conv(s(iCell,:), convGauss, 'same')-mean(s(iCell,:)))/std(s(iCell,:));
    %         sConv(iCell,:) = conv(s(iCell,:), convGauss, 'same')/max(s(iCell,:));
    %             sConv(iCell,:) = conv(s(iCell,:), convGauss, 'same');
end

% take square root of data

%% run and plot PCA
[coeff,score,latent,tsquared,explained,mu] = pca(sConv');

for iEnvIdx = 1:3
    iEnv = envIdx(iEnvIdx);
    if iDay == 8
        % 3d plot first 3 components in different colors for each env
        figure(20+iAni); hold on
        subplot(121); hold on
        plot3(score(sessIndic==iEnv,1),score(sessIndic==iEnv,2),score(sessIndic==iEnv,3),'color',col(iEnvIdx,:))
    end
    %     % compute mean/std of PCA components
    %     meanPCA(iEnv+2,iDay,:) =  mean(score(sessIndic==iEnv,1:3))./std(score(sessIndic==iEnv,1:3));
end

%% RUN IsoMap
clear X
global X %globally defined data

X = sConv;
%remove neurons with no spikes
%s = sum(X,2);
%X = X(s>0,:);

options.dims = 1:25; %1:3 %range of dimensions to calculate typ 2:20 for dim calculations, 3 for plot
options.landmarks = 1:size(X,2); %recommended value
options.display = 0
nNeighbours = 5;
[Y, R, E] = IsomapII('dfun', 'k', nNeighbours, options);


%3D plot
% correct for missing dt after isomap transformation
sessIndicIM = sessIndic(Y.index);
dimID = 3;
for iEnvIdx = 1:3
    iEnv = envIdx(iEnvIdx);
    if iDay == 8
        % 3d plot first 3 components in different colors for each env
        figure(20+iAni)
        subplot(122); hold on
        plot3(Y.coords{dimID}(1,sessIndicIM==iEnv),Y.coords{dimID}(2,sessIndicIM==iEnv),Y.coords{dimID}(3,sessIndicIM==iEnv),'color',col(iEnvIdx,:))
    end
    %     % compute mean/std of Isomap components
    %     meanIM(iEnv+2,iDay,:) = mean(Y.coords{dimID}(1:3,sessIndicIM==iEnv),2)./std(Y.coords{dimID}(1:3,sessIndicIM==iEnv),[],2);
end

%% Compute spatial average of PCA or Isomap
% nBins = 8;
% 
% envList = [-1, 0, 1];
% figure
% for iIdxEnv = 1:length(envList)
%     
%     Pos = pos{iDay}(sessIndicIM==envList(iIdxEnv),:);
%     r = Y.coords{10}(:,sessIndicIM==envList(iIdxEnv));
%     %     r = score(sessIndicIM==iEnv,1:8)';
%     
%     boundaries(:,1) = min(Pos,[],1)';
%     boundaries(:,2) = max(Pos,[],1)+0.1;
%     
%     [nCells, T] = size(r);
%     steps = diff(boundaries,1,2)/nBins;
%     
%     Occupation = nan(nBins,nBins);
%     RateMaps = nan(nBins,nBins,nCells);
%     Location = nan(T,2);
%     
%     Xbin = boundaries(1,1):steps(1):boundaries(1,2);
%     Ybin = boundaries(2,1):steps(2):boundaries(2,2);
%     
%     for x = 1:nBins
%         for y = 1:nBins
%             loc_tmp = Pos(:,1)>=Xbin(x) & Pos(:,1)<Xbin(x+1) ...
%                 & Pos(:,2)>=Ybin(y) & Pos(:,2)<Ybin(y+1);
%             Occupation(y,x) = sum(loc_tmp);
%             RateMaps(y,x,:) =  nanmean(r(:,loc_tmp),2);
%             
%             Location(loc_tmp,:) = repmat([x,y],sum(loc_tmp),1);
%         end
%     end
%     
%     for iC = 1:nCells
%         subplot(length(envList),nCells,iC+((iIdxEnv-1)*nCells))
%         imagesc(RateMaps(:,:,iC))
%     end
% end

%%
sConvAll = [];
nDim = 3
for iDim =1:nDim
    %find where a component is most active
    f = Y.coords{nDim}(iDim,:)>quantile(Y.coords{nDim}(iDim,:),0.90);
    
    % sort component by activity
    [~,sort_idx] = sort(Y.coords{nDim}(iDim,:),'descend');
    
    % extract cell activity when component is most active
    sConvDim = sConv(:,f);
    sConvAll = [sConvAll,sConvDim]; %append
    
    % sort cells by their zscore during these times when the component is highly active
    [~, cell_sort]=sort((mean(sConvDim,2)-mean(sConv,2))./std(sConv,[],2), 'descend');
    if iDim == 1 %save sort for Dim1
        cell_sort1 = cell_sort;
    end
    
    %plot times most active for each dimension, cells ranked as per dim 1
    figure(4)
    subplot(1,nDim,iDim)
    imagesc(sConvDim(cell_sort1,:))
    
    % compute correlations within times most active for each component
    figure(6)
    subplot(1,nDim,iDim)
    imagesc(corr(sConvDim))
    % imagesc(corr(sConv(cell_sort,sort_idx)))
    caxis([-.5 .5])
    
    % plot entire recording, cells sorted by activity, components sorted by
    % activity
    figure(5)
    subplot(nDim,1,iDim)
    imagesc(sConv(cell_sort1,sort_idx))
    
    
    figure(15)
    subplot(nDim,1,iDim); hold on
    % plot the mean & mean+std of each neurons to compare w/ next plots
    plot(mean(sConv,2)); plot(mean(sConv,2)+std(sConv,[],2))
    
    nCells = 10;
    for iIdxCell=1:nCells % but will use
        
        % Pick most active cell (A) in moments when dimension highly active
        iCell = cell_sort(iIdxCell);
        
        % find time when cell A is very active
        iCellHigh = sConv(iCell,:)>quantile(sConv(iCell,:),0.98);
        
        %         % sort these times
        %         [~, iCellSort] = sort(sConv(iCell,iCellHigh), 'descend');
        
        % find cells whose activity is high when cell A active
        cellPairList = find(mean(sConv(:,iCellHigh),2)>mean(sConv,2)+std(sConv,[],2));
        cellPairList(cellPairList == iCell) = []; % remove cell A from list of cells
        
        % plot times when cell A is high, cells ordered by activity in dim
        figure(30+iDim)
        subplot(1,nCells,iIdxCell)
        imagesc(sConv(cell_sort,iCellHigh))
        
        figure(35+iDim)
        subplot(1,nCells,iIdxCell)
        [~, iCellSortDim] = sort(sConvDim(iCell,:), 'descend');
        
        imagesc(sConvDim(cell_sort1,iCellSortDim(sConvDim(iCell,iCellSortDim)>0)))
        
        figure(15)
        subplot(nDim,1,iDim)
        %         % plot mean activity of each cell when cell A is high
        %         plot(mean(sConv(:,iCellHigh),2), '.-')
        %         % put an red * on cell A
        %         plot(iCell, mean(sConv(iCell,iCellHigh),2), 'r*')
        
        % or ordered
        plot(mean(sConv(cell_sort,iCellHigh),2), '.-')
        plot(iIdxCell, mean(sConv(iCell,iCellHigh),2), 'r*')
        
        figure(25)
        subplot(nDim,nCells,iIdxCell+(iDim-1)*nCells); hold on
        for iCellPair = 1:length(cellPairList)
            plot(-20:20,xcorr(s(iCell,:),s(cellPairList(iCellPair),:),20))
            title(num2str(iCell))
        end
        
        figure;
        subplot(1,length(cellPairList)+1,1)
        imagesc(sConv(cell_sort,iCellHigh))
        
        for iCellPair = 1:length(cellPairList)
            [~,xx] = sort(sConv(cellPairList(iCellPair),iCellHigh),'descend');
            f =find(iCellHigh);
            subplot(1,length(cellPairList)+1,iCellPair+1)
            imagesc(sConv(cell_sort,f(xx(1:10))))
        end
        
    end
end

%%
iCell = cell_sort(1);
iCellHigh = sConv(iCell,:)>quantile(sConv(iCell,:),0.98);

cellPairList = find(mean(sConv(:,iCellHigh),2)>mean(sConv,2)+2*std(sConv,[],2));
cellPairList(cellPairList == cell_sort(iCell)) = [];

figure;
subplot(1,length(cellPairList)+1,1)
imagesc(sConv(cell_sort,iCellHigh))

for iCellPair = 1:length(cellPairList)
    [~,xx] = sort(sConv(cellPairList(iCellPair),iCellHigh),'descend')
    f =find(iCellHigh);
    subplot(1,length(cellPairList)+1,iCellPair+1)
    imagesc(sConv(cell_sort,f(xx)))
end
%%
c = corr(sConvAll);
figure;
imagesc(c)
figure;hold on
plot(mean(c(1:180,:),1))
plot(mean(c(181:360,:),1))
plot(mean(c(361:end,:),1),'k')

%% correlation
[c, p] = corrcoef(sConv');
%% Plot correlation ribon
figure(1)
nDim = 3;
nCells = 100; %size(sConvDim,1);
cellCount = zeros(nCells,2,nDim);
for iDim =1:nDim
    %find where a component is most active
    f = Y.coords{nDim}(iDim,:)>quantile(Y.coords{nDim}(iDim,:),0.90);
    
    % sort component by activity
    [~,sort_idx] = sort(Y.coords{nDim}(iDim,:),'descend');
    
    % extract cell activity when component is most active
    sConvDim = s(:,f);
    
    [c, p] = corrcoef(sConvDim');
    cDim(:,:,iDim) = c;
    pDim(:,:,iDim) = p;
    
    [~, cell_sort]=sort((mean(sConvDim,2)-mean(s,2))./std(s,[],2), 'descend');
    
    figure(3)
    subplot(1,3,iDim); hold on
    
    for iIdx=size(sConvDim,1):-1:1
        iCell = cell_sort(iIdx);
        
        for jIdx = 1:size(sConvDim,1)
            jCell = cell_sort(jIdx);
            if p(iCell,jCell) < 0.05/nCells
%                 plot([0,1],[iIdx,jIdx], 'k--')
                if c(iCell,jCell)>0.8
                    cellCount(iIdx,1,iDim) = cellCount(iIdx) + 1;
                    if iIdx<50
                        plot([0,1],[iCell,jCell], 'b', 'LineWidth', 2)
                    else
                        plot([0,1],[iCell,jCell], 'k', 'LineWidth', 0.5)
                    end
                elseif c(iCell,jCell)<-0.8
                    cellCount(iIdx,2,iDim) = cellCount(iIdx) + 1;
                    plot([0,1],[iCell,jCell], 'r', 'LineWidth', 0.01)
                end
            end
        end
    end
    
    figure(2)
    subplot(2,3,iDim)
    hist(cellCount(:,1,iDim),1:2:max(cellCount(:,1,iDim)))
    
    subplot(2,3,iDim+3)
    if sum(cellCount(:,2,iDim))>0
        hist(cellCount(:,2,iDim),1:2:max(cellCount(:,2,iDim)))
    end
end

%%
figure
for i=1:3
subplot(2,3,i)
imagesc(cDim(:,:,i))
end

cc1 = cDim(:,:,1);
cc2 = cDim(:,:,2);
cc3 = cDim(:,:,3);

subplot(2,3,4)
plot(cc1(:),cc2(:),'.'); hold on
subplot(2,3,5)
plot(cc1(:),cc3(:),'.'); hold on
subplot(2,3,6)
plot(cc2(:),cc3(:),'.'); hold on

%%
figure
plot3(1-cc1(:),1-cc2(:),1-cc3(:), '.')









