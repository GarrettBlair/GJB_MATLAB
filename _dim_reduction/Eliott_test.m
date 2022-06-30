clear;  clc
close all;
addpath('IsoMap/');

pathAll = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/CircleSquare/analysisFiles/';

noHMC = 0;

animalList = {'M19','M20'}%,'M29', 'M35', 'M34','M39'};

col = [0,0,0;0,1,0;1,0.5,0];
col2 = [0,0,0;0.2,1,0;1,0.5,0.2];

lenGauss = 100;
sG = 1;
convGauss = gaussmf(0:2:lenGauss,[sG lenGauss/2]);

envIdx = [-1,0,1];

nAni = length(animalList);
nDay = 9;
dayPlot = 8;
dimCompute = 25;

dayMeanIM = nan(nAni,9);
dayMeanPCA = nan(nAni,9);
dayDiffIM = nan(nAni,9);
dayDiffPCA = nan(nAni,9);

for iAni = 1:nAni
    
    disp(animalList{iAni})
    
    RAll = nan(nDay, dimCompute);
    ExpVar = nan(nDay, dimCompute);
    meanPCA = nan(3,nDay,3);
    meanIM = nan(3,nDay,3);
    
    load([pathAll, '1000ms_CircleSquare_', animalList{iAni}, '.mat'])
    dayList = days; clear days; % to avoid confusion btw day and days

    for day=1:nDay
        %     b = cat(1,a{i,1,:});
        %     b2 = cat(1,a{i,2,:});
        %     c = cat(2,b,b2);
        
        dayIdx = find(dayList==day);
        disp(day)
        
        if isempty(dayIdx)
            continue
        end
        
        s = spike{dayIdx}';
        sessIndic = sess{dayIdx};
        
        if noHMC
            s = s(:,sess{dayIdx}>=0); %remove HMC
            s = s(sum(s,2)>0,:); %remove cells active only in HMC
            sessIndic = sessIndic(:,sess{dayIdx}>=0);
        end
        
        
        nCells = size(s,1);
        sConv = nan(size(s));
        for iCell = 1:nCells
            sConv(iCell,:) = sqrt((conv(s(iCell,:), convGauss, 'same')));%-mean(s(iCell,:)))/std(s(iCell,:));
            %         sConv(iCell,:) = conv(s(iCell,:), convGauss, 'same')/max(s(iCell,:));
            %             sConv(iCell,:) = conv(s(iCell,:), convGauss, 'same');
        end
        
        [coeff,score,latent,tsquared,explained,mu] = pca(sConv');
        
        figure(10+iAni)
        subplot(2,9,day+9); hold on;
        for iEnvIdx = 1:3
            iEnv = envIdx(iEnvIdx);
            plot3(score(sessIndic==iEnv,1),score(sessIndic==iEnv,2),score(sessIndic==iEnv,3),'color',col(iEnvIdx,:))
            title('PCA')
            if day == dayPlot
                figure(20+iAni); hold on
                subplot(121); hold on
                plot3(score(sessIndic==iEnv,1),score(sessIndic==iEnv,2),score(sessIndic==iEnv,3),'color',col(iEnvIdx,:))
                title('PCA')
            end
            meanPCA(iEnv+2,day,:) = [mean(score(sessIndic==iEnv,1))/std(score(sessIndic==iEnv,1)),...
                mean(score(sessIndic==iEnv,2))/std(score(sessIndic==iEnv,2)),...
                mean(score(sessIndic==iEnv,3))/std(score(sessIndic==iEnv,3))];
        end
        
        
        clear X
        global X %globally defined data
        
        X = sConv;
        %remove neurons with no spikes
        %s = sum(X,2);
        %X = X(s>0,:);
        
        options.dims = 1:dimCompute; %1:3 %range of dimensions to calculate typ 2:20 for dim calculations, 3 for plot
        options.landmarks = 1:size(X,2); %recommended value
        options.display = 0;
        nNeighbours = 5;
        [Y, R, E] = IsomapII('dfun', 'k', nNeighbours, options);
        
        %3D plot
        figure(10+iAni)
        dimID = 3; %index of dim to plot from options.dims
        subplot(2,9,day); hold on
        sessIndicIM = sessIndic(Y.index);
        for iEnvIdx = 1:3
            iEnv = envIdx(iEnvIdx);
            plot3(Y.coords{dimID}(1,sessIndicIM==iEnv),Y.coords{dimID}(2,sessIndicIM==iEnv),Y.coords{dimID}(3,sessIndicIM==iEnv),'color',col(iEnvIdx,:))
            title('Isomap')
            if day == dayPlot
                figure(20+iAni)
                subplot(122); hold on
                plot3(Y.coords{dimID}(1,sessIndicIM==iEnv),Y.coords{dimID}(2,sessIndicIM==iEnv),Y.coords{dimID}(3,sessIndicIM==iEnv),'color',col(iEnvIdx,:))
                title('Isomap')
            end
            meanIM(iEnv+2,day,:) = [mean(Y.coords{dimID}(1,sessIndicIM==iEnv))/std(Y.coords{dimID}(1,sessIndicIM==iEnv)),...
                mean(Y.coords{dimID}(2,sessIndicIM==iEnv))/std(Y.coords{dimID}(2,sessIndicIM==iEnv)),...
                mean(Y.coords{dimID}(3,sessIndicIM==iEnv))/std(Y.coords{dimID}(3,sessIndicIM==iEnv))];
        end
        
        
        RAll(day,:) = R;
        ExpVar(day,:) = cumsum(explained(1:dimCompute)/sum(explained));
        figure(2)
        subplot(1,nAni,iAni)
        plot(1-R, 'Color',[0,1-day/9,1]); hold on
        plot(cumsum(explained/sum(explained)),'Color',[1,0,1-day/9])
        xlim([0 dimCompute])
        %     plot3(Y.coords{dimID}(1,s1+1:end),Y.coords{dimID}(2,s1+1:end),Y.coords{dimID}(3,s1+1:end))
    end
    
    figure(4); hold on
    errorbar(1:dimCompute,mean(ExpVar),std(ExpVar),'r--')
    errorbar(1:dimCompute,1-mean(RAll),std(RAll),'b--')
    
    % figure
    % for iEnv=1:3
    % subplot(3,2,iEnv)
    % plot(squeeze(meanIM(iEnv,:,:)))
    %
    % subplot(3,2,iEnv+3)
    % plot(squeeze(meanPCA(iEnv,:,:)))
    % end
    
    figure(5)
    subplot(2,1,1); hold on
    % plot(sum(meanIM(:,:,:).^2,3)','b')
    % plot(sum(meanPCA(:,:,:).^2,3)','r')
    dayMeanIM(iAni,1:nDay) = mean(sum(meanIM(:,:,2:3).^2,3),1);
    dayMeanPCA(iAni,1:nDay) = mean(sum(meanPCA(:,:,2:3).^2,3),1);
    errorbar(1:nDay, mean(sum(meanIM(:,:,1:2).^2,3),1),std(sum(meanIM(:,:,1:2).^2,3),1),'b--')
    errorbar(1:nDay,mean(sum(meanPCA(:,:,1:2).^2,3),1),std(sum(meanPCA(:,:,1:2).^2,3),1),'r--')
    
    subplot(2,1,2); hold on
    dayDiffIM(iAni,1:nDay) = sum((meanIM(3,:,1:2)-meanIM(2,:,1:2)).^2,3);
    dayDiffPCA(iAni,1:nDay) = sum((meanPCA(3,:,1:2)-meanPCA(2,:,1:2)).^2,3);
    plot(sum((meanIM(3,:,1:2)-meanIM(2,:,1:2)).^2,3),'b')
    plot(sum((meanPCA(3,:,1:2)-meanPCA(2,:,1:2)).^2,3),'r')
    %%
    figure(9)
%     figure(20+iAni)
    for iDay = 1:nDay
        subplot(2,nAni,iAni); hold on
        h = sqrt(sum(meanIM(1,iDay,1:2).^2));
        r = [meanIM(1,iDay,1)/h, meanIM(1,iDay,2)/h;...
            -meanIM(1,iDay,2)/h, meanIM(1,iDay,1)/h];
        
        for iEnv =[2,3,1]
            v = squeeze(meanIM(iEnv,iDay,1:2));
            v = r*v;
            if iEnv == 2
                s = sign(v(2));
            end
            plot([0,v(1)],[0,v(2)*s],'color',col(iEnv,:))
%             plot([0,v(1)*50],[0,v(2)*50],'color',col2(iEnv,:),'LineWidth',4)
        end
        xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
        title('Isomap')
        
        subplot(2,nAni,iAni+nAni); hold on
        h = sqrt(sum(meanPCA(1,iDay,1:2).^2));
        r = [meanPCA(1,iDay,1)/h, meanPCA(1,iDay,2)/h;...
            -meanPCA(1,iDay,2)/h, meanPCA(1,iDay,1)/h];
        for iEnv =[2,3,1]
            v = squeeze(meanPCA(iEnv,iDay,1:2));
            v = r*v;
            if iEnv == 2
                s = sign(v(2));
            end
            plot([0,v(1)],[0,v(2)*s],'color',col(iEnv,:))
%             plot([0,v(1)*50],[0,v(2)*50],'color',col2(iEnv,:),'LineWidth',4)
        end
        xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
        title('PCA')
    end
    
end

%%
figure
subplot(121)
bar(1:3,[squeeze(mean(nanmean(reshape(dayMeanPCA, [6,3,3]),2))),squeeze(mean(nanmean(reshape(dayMeanIM, [6,3,3]),2)))])
hold on
plot([0.9,1,1.1, 1.9,2,2.1, 2.9,3,3.1]+0.145,dayMeanIM','k.')
errorbar([1,2,3]+0.145,squeeze(mean(nanmean(reshape(dayMeanIM, [6,3,3]),2),1)),squeeze(std(nanmean(reshape(dayMeanIM, [6,3,3]),2),1))/sqrt(6), 'ro')
errorbar([1,2,3]-0.145,squeeze(mean(nanmean(reshape(dayMeanPCA, [6,3,3]),2),1)),squeeze(std(nanmean(reshape(dayMeanPCA, [6,3,3]),2),1))/sqrt(6), 'bo')
plot([0.9,1,1.1, 1.9,2,2.1, 2.9,3,3.1]-0.145,dayMeanPCA','k.')

subplot(122)
bar(1:3,[squeeze(mean(nanmean(reshape(dayDiffPCA, [6,3,3]),2))),squeeze(mean(nanmean(reshape(dayDiffIM, [6,3,3]),2)))])
hold on
plot([0.9,1,1.1, 1.9,2,2.1, 2.9,3,3.1]+0.145,dayDiffIM','k.')
errorbar([1,2,3]+0.145,squeeze(mean(nanmean(reshape(dayDiffIM, [6,3,3]),2),1)),squeeze(std(nanmean(reshape(dayDiffIM, [6,3,3]),2),1))/sqrt(6), 'ro')
errorbar([1,2,3]-0.145,squeeze(mean(nanmean(reshape(dayDiffPCA, [6,3,3]),2),1)),squeeze(std(nanmean(reshape(dayDiffPCA, [6,3,3]),2),1))/sqrt(6), 'bo')
plot([0.9,1,1.1, 1.9,2,2.1, 2.9,3,3.1]-0.145,dayDiffPCA','k.')