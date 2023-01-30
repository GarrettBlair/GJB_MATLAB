%%
file_dirs = {'E:/Miniscope Data/mouse_cw1/2021_04_09/09_23_46/',... % 5044 frames, running, none
    'E:/Miniscope Data/mouse_cw1/2021_04_09/09_31_26/',... % 4104 frames, sounds & scents 61726-92188 campfire,121621-200205 despacito, 221302-249368? (& -265688 & -270345) windex, 302208-331588 eoh70, 347739 claps
    'E:/Miniscope Data/mouse_cw2/2021_04_09/08_48_06/',... % 4187 frames, running, none
    'E:/Miniscope Data/mouse_cw2/2021_04_09/09_00_11/',... % 2724 frames, running, none
    'E:/Miniscope Data/mouse_cw2/2021_04_09/09_04_40/',... % 2983 frames, sounds 62440-91133 natural sounds, 129778-241354 despacito
    'E:/Miniscope Data/mouse_cw2/2021_04_09/09_09_59/'}; % 1237 frames, scents 22060-42298 windex, 61606-78970 eoh70, 91918 clap

file_title = {'mouse_cw1_rec1',... % 5044 frames, running, none
    'mouse_cw1_rec2',... % 4104 frames, sounds & scents 61726-92188 campfire,121621-200205 despacito, 221302-249368? (& -265688 & -270345) windex, 302208-331588 eoh70, 347739 claps
    'mouse_cw2_rec1',... % 4187 frames, running, none
    'mouse_cw2_rec2',... % 2724 frames, running, none
    'mouse_cw2_rec3',... % 2983 frames, sounds 62440-91133 natural sounds, 129778-241354 despacito
    'mouse_cw2_rec4'}; % 1237 frames, scents 22060-42298 windex, 61606-78970 eoh70, 91918 clap

mask_fname = {'F:\ScaryData2\Guo Blair Methods data\crainial window\vessels\cw1 mask.tif',... % 5044 frames, running, none
    'F:\ScaryData2\Guo Blair Methods data\crainial window\vessels\cw1 mask.tif',... % 4104 frames, sounds & scents 61726-92188 campfire,121621-200205 despacito, 221302-249368? (& -265688 & -270345) windex, 302208-331588 eoh70, 347739 claps
    'F:\ScaryData2\Guo Blair Methods data\crainial window\vessels\cw2 mask.tif',... % 4187 frames, running, none
    'F:\ScaryData2\Guo Blair Methods data\crainial window\vessels\cw2 mask.tif',... % 2724 frames, running, none
    'F:\ScaryData2\Guo Blair Methods data\crainial window\vessels\cw2 mask.tif',... % 2983 frames, sounds 62440-91133 natural sounds, 129778-241354 despacito
    'F:\ScaryData2\Guo Blair Methods data\crainial window\vessels\cw2 mask.tif'}; % 1237 frames, scents 22060-42298 windex, 61606-78970 eoh70, 91918 clap

cmn_name = 'caiman_cnmfe_out.mat';
saveDir = 'F:\ScaryData2\Guo Blair Methods data\crainial window\data\';
temporalDS = 2;
mm_per_pix = (200-80)/10;
plotting = true;
for sessLoop = 1%[1 3]:length(file_dirs)
    savefile = [saveDir file_title{sessLoop} '.mat'];
    if exist(savefile, 'file')~=2
        fprintf('\nRunning behavior on: \n%s\n\t\t', [file_dirs{sessLoop} 'MiniCam0/'])
        [behav] = extract_wheel_speed([file_dirs{sessLoop} 'MiniCam0/'], temporalDS, mm_per_pix, plotting);
        save( [file_dirs{sessLoop} 'MiniCam0/behav.mat'], 'behav')
        fprintf('Done!\n')
        
        fprintf('\nExtracting caiman data from: \n%s\n\t\t', [file_dirs{sessLoop}, 'Miniscope/' cmn_name])
        neuron = load([file_dirs{sessLoop}, 'Miniscope/' cmn_name]);
        numA = floor(size(neuron.A,2)/double(neuron.segs_per_fullA));
        A = []; for i = 1:numA+1; eval(sprintf('A=cat(2, A, neuron.fullA%d);', (i-1)*neuron.segs_per_fullA)); end
        neuron.fullA = A;
        contours = gbContours(neuron.fullA, neuron.dims, [], .5);
        
        mask_im = imread(mask_fname{sessLoop});
        [circularity, area, centroids] = gbContours_circularity(contours, [], maks_im);
        
        %     [mask_out] = gb_blood_vessel_detect(neuron.minFrame./255, 7, contours);
        
        ts_file = [file_dirs{sessLoop}, 'Miniscope\timestamps.csv'];
        tfile = xlsread(ts_file);
        frameNum = tfile(:,1)+1;
        time = tfile(1:temporalDS:end,2);
        dt = (( time(end)-time(1) )/1000 ) /length(time);
        spd = interp1(behav.time, behav.speed, time, 'linear');
        spd(isnan(spd)) = interp1(find(~isnan(spd)), spd(find(~isnan(spd))), find(isnan(spd)), 'nearest', 'extrap');
        dist = cumsum(spd./(1/dt));
        
        craw = double(neuron.C+neuron.YrA);
%         ispks = Tad_spk_inference(craw, 2.5, false)>0;
        ispks = neuron.spks;

        fprintf('Done!\n')
        %     good = area <= quantile(area, .85) & circularity >= quantile(circularity, .05);
        %     goodim = squeeze(sum(contours(good,:,:), 1));
        %     badim = squeeze(sum(contours(~good,:,:), 1));
        %     figure; image(cat(3, badim>0, goodim>0, goodim*0))
        
        fprintf('\nRuning speed corr analysis... ')
        nrand = 100;%500;
        spdCorr = NaN(size(craw,1), 1);
        spdCorr_rand = NaN(size(craw,1), nrand);
        spdP = NaN(size(craw,1), 1);
        spdChance = NaN(size(craw,1), 1);
        spdP_rand = NaN(size(craw,1), nrand);
        tic
        for segLoop = 1:size(craw,1)
            s = (neuron.spks(segLoop,:)>0)';
            randds = rand(nrand, length(s));
            [spdCorr(segLoop), spdP(segLoop)] = nancorr(spd, s);
            for randLoop = 1:nrand
                [~, rand_ord] = sort(randds(randLoop, :));
                [spdCorr_rand(segLoop, randLoop), spdP_rand(segLoop, randLoop)] = nancorr(spd, s(rand_ord));
            end
            spdChance(segLoop) = sum(spdP(segLoop) > spdP_rand(segLoop,:))/nrand;
        end
        toc
        fprintf('\n\t\tDone!\n')
        
        [~, spk_ord] = sort(ms.spdChance, 'ascend');
        %     [~, spk_ord] = sort(spdCorr, 'descend');
        spks_ord = neuron.spks(spk_ord, :)>0;
        c_colors = viridis(size(neuron.fullA,2));
        [~, ord] = sort(rand(size(c_colors,1), 1));
        cim = gb_make_color_contours(contours, [], c_colors(ord,:), true);
        figure;
        set(gcf, 'Name', file_title{sessLoop}, 'Position', [300 109 725 1199])
        subplot_tight(3,2,1, [.025 .1]);
        imshow(neuron.meanFrame./255)
        subplot_tight(3,2,2, [.025 .1]);
        image(cim); axis image off
        subplot_tight(3,2, 3:4, [.025 .1]);
        yyaxis('left')
        plot(ms.time./1000, ms.spd); axis tight
        ylabel('Wheel speed (cm/sec)');
        xlabel('Time (sec)')
        yyaxis('right')
        plot(ms.time./1000, ms.dist); axis tight
        ylabel('Cumulative distance (cm)');
        %     xlim([0 60])
        subplot_tight(3,2,5:6, [.025 .1]);
        imagesc(spks_ord)
        
        ms.datadir = [file_dirs{sessLoop}, 'Miniscope/'];
        ms.frameNum = frameNum;
        ms.time = time;
        ms.spd = spd;
        ms.dist = dist;
        ms.dt = dt;
        ms.spdCorr = spdCorr;
        ms.spdP = spdP;
        ms.spdChance = spdChance;
        ms.spd_rand_samples = nrand;
%         neuron.spks = ispks;
        neuron.craw = craw;
        neuron.circularity = circularity;
        neuron.area = area;
        neuron.centroids = centroids;
        
        savefile = [file_dirs{sessLoop}, 'Miniscope/' file_title{sessLoop} '.mat'];
        save(savefile, 'ms', 'neuron', 'behav', 'contours', '-v7.3')
        
        savefile = [saveDir file_title{sessLoop} '.mat'];
        save(savefile, 'ms', 'neuron', 'behav', 'contours', '-v7.3')
        figsave = [saveDir 'fig_' file_title{sessLoop} '.fig'];
        saveas(gcf, figsave)
        
    else
        %%
        fprintf('Loading existing file: %s\n', savefile);
        savefile = [saveDir file_title{sessLoop} '.mat'];
        
        load(savefile)
        if ~isfield('newC', neuron)
            [neuron.newC2, neuron.spks2, ~] = gb_cnmfe_deconv(neuron.C+neuron.YrA, 'ar1', false);
            save(savefile, 'neuron', '-append')
            ms = rmfield(ms,  'spdRate_prob');
        end
        %%
        if ~isfield(ms, 'spdRate_prob')
        mask_im = imread(mask_fname{sessLoop});
%         mask_im = uint8(squeeze(nansum(mask_im, 3)));
        mask_im = uint8(squeeze(mask_im(:,:,1)))<=30;
%         mask_im = imresize(mask_im, 2);
        [circularity, area, centroids, prop_mask, mask_aligned] = gbContours_circularity(contours, [], mask_im, neuron.minFrame);
        
        neuron.prop_mask = prop_mask;
        neuron.mask_aligned = mask_aligned;

        
%         save(savefile, 'neuron', '-append')
        
        good_seg = neuron.prop_mask<=.01;
        goodim = squeeze(sum(contours(good_seg,:,:), 1));
        badim = squeeze(sum(contours(~good_seg,:,:), 1));
        m = double(neuron.minFrame);
        m = double(neuron.minFrame./255);
        r = .95*m + 4*double(badim);
        g = .95*m + 4*double(goodim);
        b = .95*m + .5*double(neuron.mask_aligned);
        im = cat(3, r, g, b);
        figure; subplot_tight(1,1,1, [0 0]); image(1*im); axis image off
        
%         [~, spk_ord] = sort(ms.spdChance(good), 'ascend');

        nrand = 500;
        spdRate_prob = NaN(size(neuron.spks,1), 1);
        rspdRate = NaN(size(neuron.spks,1), nrand);
        rstillRate = NaN(size(neuron.spks,1), nrand);
        tic
        params.num_partitions = 3;
        params.max_spd_thresh = .5;
        params.min_spd_thresh = .1;
        params.min_samples = 5;
        spd_epochs = get_speed_epochs(ms.spd, params);
        spdRate = (sum(neuron.spks(:, spd_epochs>0)>0, 2)./sum(spd_epochs>0)).*ms.dt;
        stillRate = (sum(neuron.spks(:, spd_epochs==0)>0, 2)./sum(spd_epochs==0)).*ms.dt;
        randds = rand(nrand, size(neuron.spks,2));
        for randLoop = 1:nrand
            [~, rand_ord] = sort(randds(randLoop, :));
            s = neuron.spks(:, rand_ord);
            rspdRate(:, randLoop) = (sum(s(:, spd_epochs>0)>0, 2)./sum(spd_epochs>0)).*ms.dt;
            rstillRate(:, randLoop) = (sum(s(:, spd_epochs==0)>0, 2)./sum(spd_epochs==0)).*ms.dt;
%             a = (sum(s(:, spd_epochs>0)>0, 2)./sum(spd_epochs>0)).*ms.dt;
%             b = (sum(s(:, spd_epochs==0)>0, 2)./sum(spd_epochs==0)).*ms.dt;
%             rspdRate(:, randLoop) = a - b;
        end
        for segLoop = 1:size(neuron.spks,1)
            a = spdRate(segLoop)-stillRate(segLoop);
            b = rspdRate(segLoop,:)-rstillRate(segLoop,:);
            spdRate_prob(segLoop) = 1 - ( sum(a>b)/nrand );
        end
        toc
        ms.spdRate=spdRate;
        ms.stillRate=stillRate;
        ms.spdRate_prob=spdRate_prob;
        ms.rspdRate=rspdRate;
        ms.rstillRate = rstillRate;
        ms.speed_epochs = spd_epochs;
        ms.speed_params = params;
        end
        %%
        figure;
        %
%         goodS = ms.spdChance<=.05 & ms.spdCorr>0 & ms.spdP<=.05 & good_seg;
        good_seg = neuron.prop_mask<=.01;
        goodS = ms.spdRate_prob<.05 & good_seg;
        badS = ms.spdRate_prob>=.05 & good_seg;
%         [~, spk_ord] = sort(ms.spdP(good_seg), 'ascend');
        [~, spk_ord] = sort(ms.spdRate_prob(goodS), 'ascend');
        spks = neuron.spks(goodS, :)>0;
        spks_ord = spks(spk_ord, :);
%         c_colors = plasma(2*sum(goodS));
        c_colors = (1-ms.spdRate_prob(goodS))*[0 1 1];
        c_colors(:,[ 2 3] ) = normalize_matrix(c_colors(:,[2, 3]), 1);
        b_colors = ones(sum(badS), 1)*[1 0 0];
        m = double(neuron.minFrame./255);
%         cim = gb_make_color_contours(contours(goodS,:,:), m*.5, c_colors(spk_ord,:), true);
        gim = gb_make_color_contours(contours(goodS,:,:), m*.5, c_colors, true);
        badim = gb_make_color_contours(contours(badS,:,:), m*.5, b_colors, true);
        
        cim = (gim+badim./2);
        %
        xl = [0 360];%[120 240];
        conv_kern = [zeros(1,4), ones(1,5)]; conv_kern = conv_kern./sum(conv_kern);
        s = nansum(normalize_matrix(neuron.spks(goodS,:)>0, 1))./sum(goodS);
        s = 100*s;%conv(s, ones(ssize,1)./ssize, 'same');
        set(gcf, 'Name', file_title{sessLoop}, 'Position', [300 109 1057 1199])
        clf
        subplot_tight(6,2,[1, 3], [.025 .025]); cla
%         image(im); axis image off
        imshow(double(neuron.maxFrame./255)); axis image off
        subplot_tight(6,2, [2 4], [.025 .025]); cla
        image(cim); axis image off
        subplot_tight(6,2, 5:6, [.0025 .05]); cla
        yyaxis('left')
        plot(ms.time./1000, ms.spd, 'k-', 'LineWidth', 2); axis tight
        ylabel('Wheel speed (cm/sec)');
        ylim([-.1 4]);
        set(gca, 'YColor', 'k')
        yyaxis('right')
        plot(ms.time./1000, s, '-', 'Color', [0 .6 .6], 'LineWidth', 1); axis tight
        ylabel('Speed cells active (%)');
        
        x = [1:length(ms.time)];
        t = ms.time./1000 - ms.time(find(ms.time./1000 >= xl(1), 1))./1000;
        tchange = diff(mod(t,60))<0;
        xchange = x(tchange);
        tlabel = round(t(tchange));
        set(gca, 'XTick', ms.time(tchange)./1000, 'XTickLabel', tlabel, 'YColor', [0 .6 .6]./2)
        xlim(xl)

        subplot_tight(6,2,7:12, [.05 .05]); cla
        spk = conv2(1, conv_kern, spks_ord, 'same')>0;
        bspk = conv2(1, conv_kern, neuron.spks(badS, :), 'same')>0;
        
        spk_im = ones(size([spk; bspk],1), size([spk; bspk],2), 3);
        spk_im(1:size(spk,1),:,1) = spk_im(1:size(spk,1),:,1) - spk;
        spk_im(1:size(spk,1),:,2) = spk_im(1:size(spk,1),:,2) - spk*.4;
        spk_im(1:size(spk,1),:,3) = spk_im(1:size(spk,1),:,3) - spk*.4;
        spk_im(size(spk,1)+1:end,:,2) = spk_im(size(spk,1)+1:end,:,2) - bspk/2;
        spk_im(size(spk,1)+1:end,:,3) = spk_im(size(spk,1)+1:end,:,3) - bspk/2;
%         imagesc([spk; bspk]); colormap(1-bone)
        image(spk_im)
        set(gca, 'XTick', xchange, 'XTickLabel', tlabel)
        xl2 = [find(ms.time./1000 >= xl(1), 1) find(ms.time./1000 >= xl(2), 1)];
        xlim(xl2)
        xlabel('Time (sec)')
        %%
        
%         savefile = [file_dirs{sessLoop}, 'Miniscope/' file_title{sessLoop} '.mat'];
%         save(savefile, 'neuron', 'ms', '-append')
        
        savefile = [saveDir file_title{sessLoop} '.mat'];
        save(savefile, 'neuron', 'ms', '-append')

%         figsave = [saveDir 'fig_' file_title{sessLoop} '.png'];
        figsave = [saveDir 'fig_' file_title{sessLoop} '.fig'];
        saveas(gcf, figsave)
    end
end
%         [~, spk_ord] = sort(ms.spdP(goodS), 'ascend');
% LFOV_behav_vid(fname, ms, behav, neuron, 290, 340, 4, goodS, spk_ord)
% m = neuron.minFrame./255;
% hBig = ones(20)./(20^2); %strel('disk', 100);
% mg = filter2(hBig, m,'same');
% ms = m-mg;
% ms = ms-min(ms(:));
% ms = 255.*ms./max(ms(:));
% 
% 
% neuron = load('E:\Miniscope Data\mouse_cw1\2021_04_09\09_23_46\Miniscope\caiman_cnmfe_out.mat');
% 
% contours = gbContours(neuron.fullA, neuron.dims, [], .5);
% 
% 
% contours = gbContours(A, neuron.dims, [], .25);
% [circularity, area, centroids] = gbContours_circularity(contours, []);
% contours = gbContours(A, neuron.dims, [], .7);
% c_colors = viridis(size(A,1));
% [~, ord] = sort(rand(size(c_colors,1), 1));
% cim = gb_make_color_contours(contours, [], c_colors(ord,:), true);

%%
bg = neuron.meanFrame./255;
[inside_inds, outside_inds, contours, bounds] = Draw_contour_bounding(neuron.fullA, neuron.dims, bg, []);
neuron.good_segs = outside_inds;%intersect(good_inds, ms.neuron.idx_components+1);
neuron.bad_segs = inside_inds;%union(bad_inds, ms.neuron.idx_components_bad+1);
contours = contours(neuron.good_segs,:, :);


% neuron.ispks = Tad_spk_inference(neuron.C+neuron.YrA, 2.5, false);
% neuron.ispks = 1*(neuron.ispks>0);

s = neuron.spks>0;
conv_s = zeros(size(s));
nsegs = size(s,1);
for i = 1:nsegs
    conv_s(i,:) = conv(s(i,:), ones(3,1)./3, 'same');
end
%%

figure;
hold on
% cr = normalize_matrix(neuron.C(neuron.good_segs,:)+neuron.YrA(neuron.good_segs,:));
cr = normalize_matrix(neuron.C(neuron.good_segs,:));
act_im = cr + 2*s;

for i = 1:size(s,1)
    plot(cr(i,:)*.9 + i);
    scatter(find(s(i,:)), cr(i,s(i,:)>0).*.9 + i, 'r.');
    
end

save('F:\ScaryData2\Guo Blair Methods data\crainial window\cranial_win_1.mat', 'neuron', '-v7.3');


%%
fs = 1:4:1000;
fname = 'F:\ScaryData2\Guo Blair Methods data\crainial window\msCam_MC_1.tiff';
height = double(neuron.dims(1)); width = double(neuron.dims(2));


%%
Yim = cell(length(fs),1);
bim = zeros(length(fs), height, width);
colormap = .5.*rand(length(neuron.good_segs), 3) + .5;
for i = 1:length(fs)
    %%
   bim(i,:,:) = double(imread(fname, fs(i)))./255;
   act = zeros(height, width, 3);
   for j = 1:3
       a = neuron.fullA(:,neuron.good_segs)*(cr(:,fs(i)).*colormap(:,j));
       act(:,:,j) = reshape(a, [height width]);
   end
   Yim{i} = act;%cat(3, act{1}, act{2}, act{3});
end
%
figure(88); clf
set(gcf, 'Color', 'k')
set(gcf,'Position', [1000         399         width/2         height/2])
v = VideoWriter('F:\ScaryData2\Guo Blair Methods data\crainial window\test.avi');
v.Quality = 80;
v.FrameRate = 30;
v.open;
for i = 1:length(fs)
   image(Yim{i}*5+ squeeze(bim(i,:,:))*.5);%,[30 220]);
   axis image off
   drawnow
   temp = getframe(gcf);
   v.writeVideo(temp.cdata);
end
v.close



%% plot rate and speed %%%%%%%%%%

t1 = load('F:\ScaryData2\Guo Blair Methods data\crainial window\data\mouse_cw1_rec1.mat');
t2 = load('F:\ScaryData2\Guo Blair Methods data\crainial window\data\mouse_cw2_rec1.mat');
%%
s = [t1.ms.spd; t2.ms.spd];
% r = [t1.ms.spdRate; t2.ms.spdRate];


params.num_partitions = 3;
params.max_spd_thresh = .5;
params.min_spd_thresh = .1;
params.min_samples = 5;
% spd_epochs = get_speed_epochs(t1.ms.spd, params);


%
step = .25;
sub_speed = .5:step:3.5;%4;
figure(109); clf;

subplot(1,3,1)
hold on
histogram(t1.ms.spd, sub_speed, 'Normalization', 'probability', 'FaceColor', [0 0 .8], 'FaceAlpha', .3)
histogram(t2.ms.spd, sub_speed, 'Normalization', 'probability', 'FaceColor', [0 .8 0], 'FaceAlpha', .1)
set(gca,'YTick', [0:.05:.1], 'XTick', round(10*linspace(sub_speed(1), sub_speed(end), 4))/10,...
    'XTickLabel', round(10*linspace(sub_speed(1), sub_speed(end), 4))/10)
axis([0 4 0 .1])
title('Speed dist. (b=mouse1, g=mouse2)')
ylabel('Probability')
xlabel('Speed')

spd_epochs = get_speed_epochs(t1.ms.spd, params);
goodc1 = t1.ms.spdRate_prob<=.05 & t1.neuron.prop_mask<=.01;
notgoodc1 = t1.ms.spdRate_prob>.05 & t1.neuron.prop_mask<=.01;
for i = 2:length(sub_speed)
    inds = spd_epochs>0 & t1.ms.spd<=sub_speed(i) & t1.ms.spd>sub_speed(i-1);
    spdRate = (sum(t1.neuron.spks(:, inds)>0, 2)./sum(inds)).*t1.ms.dt;
    
%     subplot(1,3,2); hold on
%     scatter(spdRate(goodc & spdRate>0)*0 + sub_speed(i) - .0125/2,  spdRate(goodc & spdRate>0), 'o',...
%         'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0 1], 'MarkerFaceAlpha', .1)
%     subplot(1,3,3); hold on
%     scatter(spdRate(~goodc & spdRate>0)*0 + sub_speed(i) - .0125/2,  spdRate(~goodc & spdRate>0), 'x',...
%         'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.6 0 .6], 'MarkerFaceAlpha', .1)
    subplot(1,3,2); hold on
    x = sub_speed(i) - step/2 - step/8;
    inds = goodc1 & ~isnan(spdRate);
    y = spdRate(inds);% & spdRate>0);
    bar(x, mean(y), 'BarWidth', step/4, 'FaceColor', 'b')
    errorbar(x, mean(y), [], std(y), 'b')
    subplot(1,3,3); hold on
    inds = notgoodc1 & ~isnan(spdRate);
    y = spdRate(inds);
    bar(x, mean(y), 'BarWidth', step/4, 'FaceColor', 'b')
    errorbar(x, mean(y), [], std(y), 'b')

%     scatter(spdRate(~goodc & spdRate>0)*0 + sub_speed(i) - .0125/2,  spdRate(~goodc & spdRate>0), 'o',...
%         'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.6 0 .6], 'MarkerFaceAlpha', .1)
end
drawnow

spd_epochs = get_speed_epochs(t2.ms.spd, params);
goodc2 = t2.ms.spdRate_prob<=.05 & t2.neuron.prop_mask<=.01;
notgoodc2 = t2.ms.spdRate_prob>.05 & t2.neuron.prop_mask<=.01;
for i = 2:length(sub_speed)
    inds = spd_epochs>0 & t2.ms.spd<=sub_speed(i) & t2.ms.spd>sub_speed(i-1);
    spdRate = (sum(t2.neuron.spks(:, inds)>0, 2)./sum(inds)).*t2.ms.dt;
    
%     subplot(1,3,2); hold on
%     scatter(spdRate(goodc & spdRate>0)*0 + sub_speed(i) - .0125/2,  spdRate(goodc & spdRate>0), 'o',...
%         'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 1 1], 'MarkerFaceAlpha', .1)
%     subplot(1,3,3); hold on
%     scatter(spdRate(~goodc & spdRate>0)*0 + sub_speed(i) - .0125/2,  spdRate(~goodc & spdRate>0), 'o',...
%         'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 .6 .6], 'MarkerFaceAlpha', .1)
    subplot(1,3,2); hold on
    x = sub_speed(i) - step/2 + step/8;
    inds = goodc2 & ~isnan(spdRate);
    y = spdRate(inds);% & spdRate>0);
    bar(x, mean(y), 'BarWidth', step/4, 'FaceColor', 'g')
    errorbar(x, mean(y), [], std(y), 'g')
    subplot(1,3,3); hold on
    inds = notgoodc2 & ~isnan(spdRate);
    y = spdRate(inds);
    bar(x, mean(y), 'BarWidth', step/4, 'FaceColor', 'g')
    errorbar(x, mean(y), [], std(y), 'g')
end
%
subplot(1,3,2); hold on
% axis([0 4 0 .008])
axis([0 4 0 .03])
title(sprintf('Speed Modulated cells (p<=.05)\n (b) n=%d , (g) n=%d', sum(goodc1), sum(goodc2)))
set(gca,'YTick', [0:.01:.03], 'XTick', round(10*linspace(sub_speed(1), sub_speed(end), 4))/10,...
    'XTickLabel', round(10*linspace(sub_speed(1), sub_speed(end), 4))/10)
xlabel('Speed')
ylabel('Ca^{2+} rate')
subplot(1,3,3); hold on
title(sprintf('Others cells (p>.05)\n (b) n=%d , (g) n=%d', sum(notgoodc1), sum(notgoodc2)))
set(gca,'YTick', [0:.01:.03], 'XTick', round(10*linspace(sub_speed(1), sub_speed(end), 4))/10,...
    'XTickLabel', round(10*linspace(sub_speed(1), sub_speed(end), 4))/10)
% axis([0  4 0 .008])
axis([0  4 0 .03])
xlabel('Speed')
ylabel('Ca^{2+} rate')


drawnow

saveas(gcf, 'F:\ScaryData2\Guo Blair Methods data\crainial window\data\CW_speed_rate_relationship.fig')






















