% Y = load_tiffstack("E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_12\16_55_19_RET8\HPC_miniscope1\msCam_MC_cawaves.tiff");
% mc_name = "E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_12\16_55_19_RET8\HPC_miniscope1\msCam_MC.tiff";
clear
%%%%%%%% first retrieval 8
dirrr = {'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_07_25\17_44_28_TR7\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_12\16_55_19_RET8\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_13\16_28_05_HC9\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_13\17_35_20_RET10\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_13\17_51_12_CON11\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_15\14_26_24_HC15\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_15\15_16_02_CON16\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_15\16_28_51_HC17\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_16\12_19_41_CON18\HPC_miniscope1\',...%%%%
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_07_25\17_11_23_TR7\HPC_miniscope1\',...%%%%
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_08_12\17_39_31_RET8\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_08_13\16_57_02_HC9\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_08_13\18_35_09_RET10\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_08_13\18_46_02_CON11\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_08_15\14_48_54_HC15\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_08_15\15_49_46_CON16\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_08_15\16_58_04_HC17\HPC_miniscope1\',...
    'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_08_16\12_46_12_CON18\HPC_miniscope1\'};

ndir = length(dirrr);
%%
videoPlotting = false;
%%%%%% evaluates ca waves
if false
for dirLoop = 1:ndir
    dirrr{dirLoop}
    Y = load_tiffstack([dirrr{dirLoop} 'msCam_MC_cawaves.tiff']);
    data_name = [dirrr{dirLoop} 'ca_waves.mat'];
    if videoPlotting == true
        outname = [dirrr{dirLoop} 'ca_processing.tiff'];
        Yraw = load_tiffstack([dirrr{dirLoop} 'msCam_MC.tiff']);
        minbg = min(Yraw,[],3);
    end
    ys = Y(1:2:end,1:2:end,:);
    [h, w, nf] = size(ys);
    bg = mean(ys,3);
    bgs = std(single(ys), [], 3);
    mask = bgs<1 | bg<1;
    bgs(mask) = Inf;
    bgtime = reshape(bg, [h*w, 1])*ones(1, nf);
    bgstime = reshape(bgs, [h*w, 1])*ones(1, nf);
    td = 30*5; % smoothing in time
    yss = reshape(ys, [h*w, nf]);
    yss = conv2(1, ones(td,1)./td, yss, 'same');
    yss = (yss - bgtime )./bgstime;
    nf_ds = size(yss,2);
    
    yim = reshape(yss, [h, w, nf_ds]);
    %%
    z_score_thresh = 2;
    minarea = 10;
    maxarea = ceil(h*w/4);
    diff_thresh = 10;
    figure(2); clf;
    set(gcf, 'Position', [180   733   719   221], 'Color', [.8 .8 .8])
    colormap viridis
    
    % fs = 3000:2:(3000+(10*700));%nf;%30000;%850:4:2200;
    fs = 1:nf;%nf;%30000;%850:4:2200;
    % fs = 7400:8800;
    
    spd = single(NaN(length(fs), 10));
    cx = single(NaN(length(fs), 10));
    cy = single(NaN(length(fs), 10));
    area = single(NaN(length(fs), 10));
    pm = single(NaN(length(fs), 10));
    drawinds = -1; % fs(30:30:end);%1:30:nf;
    %
    area_clrs = light_colormap(jet(16), 3);
    for ind = 1:length(fs) % 1000:2000%nf_ds
        %%
        i = fs(ind);
        a = squeeze(yim(:,:,i)) >= z_score_thresh;
        if ismember(i, drawinds)
            %%
            clf
            im1 = squeeze(Yraw(:,:,i)-minbg);
            %         im1 = imresize(im1, [h, w], 'bilinear');
            subplot_tight(1,3,1);
            imagesc(im1, [0 40]);%, [0 1]);
            title(sprintf('DF/F, frame %6.0f', i))
            colorbar
            axis image
            im2 = squeeze(yim(:,:,i));
            
            subplot_tight(1,3,2);
            imagesc(im2, [-2 8]);%, [0 1]);
            hold on
            title(sprintf('Morph. filtered,\n z-score'))
            colorbar
            axis image
            
            subplot_tight(1,3,3);
            imagesc(a, [0 1]);%, [0 1]);
            title('Ca ROIs ')
            colorbar
            axis image
        end
        temp = regionprops(a, 'Centroid', 'ConvexArea', 'Perimeter', 'ConvexHull');
        num_regions = 0;
        clearvars valid_region
        for tt = 1:length(temp)
            if temp(tt).ConvexArea >= minarea && temp(tt).ConvexArea <= maxarea
                num_regions  = num_regions +1;
                valid_region(num_regions) = temp(tt);
            end
        end
        if num_regions > 0
            if num_regions>1
                bestmatch_idx = NaN(num_regions ,1);
                dx = []; dy = []; da = []; dp = [];
                try
                    prev_regions = sum(~isnan(cx(ind-1, :)));
                catch
                    prev_regions = 0;
                end
                for ttt = 1:num_regions
                    for rrr = 1:prev_regions
                        dx(ttt, rrr) = valid_region(ttt).Centroid(1) - cx(ind-1, rrr);
                        dy(ttt, rrr) = valid_region(ttt).Centroid(2) - cy(ind-1, rrr);
                        da(ttt, rrr) = sqrt(abs(valid_region(ttt).ConvexArea - area(ind-1, rrr)));
                        dp(ttt, rrr) = valid_region(ttt).Perimeter - pm(ind-1, rrr);
                    end
                end
                dd = sqrt(dx.^2 + dy.^2 + da.^2 + dp.^2);
                labels = 1:num_regions;
                if prev_regions >= num_regions
                    label_loop = num_regions;
                else
                    label_loop = prev_regions;
                end
                for rrr = 1:label_loop
                    if nanmin(dd(:,rrr)) < diff_thresh
                        best = find(dd(:,rrr)==nanmin(dd(:,rrr)));
                        bestmatch_idx(best) = rrr;
                    end
                end
                bestmatch_idx(isnan(bestmatch_idx))  = setdiff(labels, bestmatch_idx);
            else
                bestmatch_idx = 1;
            end
            for tt = 1:num_regions
                this_ind = bestmatch_idx(tt);
                if ind>1 && ~isnan(cx(ind-1))
                    dx = valid_region(tt).Centroid(1) - cx(ind-1, this_ind);
                    dy = valid_region(tt).Centroid(2) - cy(ind-1, this_ind);
                    dd = sqrt(dx.^2 + dy.^2);
                    
                    if dd<diff_thresh || isnan(dd)
                        spd(ind, this_ind) = dd;
                        cx(ind, this_ind) = valid_region(tt).Centroid(1);
                        cy(ind, this_ind) = valid_region(tt).Centroid(2);
                        area(ind, this_ind) = valid_region(tt).ConvexArea;
                        pm(ind, this_ind) = valid_region(tt).Perimeter;
                    end
                else
                    %                 idx = idx+1;
                    spd(ind, this_ind) = 0;
                    cx(ind, this_ind) = valid_region(tt).Centroid(1);
                    cy(ind, this_ind) = valid_region(tt).Centroid(2);
                    area(ind, this_ind) = valid_region(tt).ConvexArea;
                    pm(ind, this_ind) = valid_region(tt).Perimeter;
                end
                if ismember(i, drawinds)
                    hold on
                    scatter(valid_region(tt).Centroid(1), valid_region(tt).Centroid(2), 'MarkerFaceColor', area_clrs(this_ind,:))
                    plot(valid_region(tt).ConvexHull(:,1), valid_region(tt).ConvexHull(:,2), '-', 'Color', area_clrs(this_ind,:))
                    
                    a = find(isnan(cx(1:ind, this_ind)), 1, 'last');
                    if ind>a
                        plot(cx(a:ind, this_ind), cy(a:ind, this_ind), '-', 'Color', area_clrs(this_ind,:))
                    end
                end
            end
            last_region = valid_region;
        end
        drawnow;
        if ismember(i, drawinds) && videoPlotting
            temp = getframe(gcf);
            imwrite(temp.cdata, outname,'WriteMode','append')
        end
        %     pause(.01)
    end
    active_roi = nanmax(area,[],2)>0;
    roi_perc(dirLoop) = nanmean(active_roi);
    figure; plot(fs, area); xlim([1 nf]); title(sprintf('%% ROI = %3.3f',100*nanmean(nanmax(area,[],2)>0))) 
    drawnow
    save(data_name, 'cx', 'cy', 'spd', 'area', 'pm', 'fs', 'dirrr')
    
end %%%%%%%%%%%%%%%%%%%%%%%
end
%% LOAD IT
roi_perc = NaN(ndir,1);
roi1 = NaN(ndir,1);
roi2 = NaN(ndir,1);
for dirLoop = 1:ndir
    %%
    data_name = [dirrr{dirLoop} 'ca_waves.mat'];
    disp(data_name)
    load(data_name, 'cx', 'cy', 'spd', 'area', 'pm', 'fs', 't');
    nf = length(fs);
    active_roi = nanmax(area,[],2)>0;
    roi_perc(dirLoop) = nanmean(active_roi);
    split = find(t >= t(end)/2, 1);
    roi1(dirLoop) = nanmean(active_roi(1:split));
    roi2(dirLoop) = nanmean(active_roi(split+1:end));
%     figure; plot(t, active_roi); title(sprintf('%% ROI = %3.3f',100*nanmean(nanmax(area,[],2)>0))); 
%     title(data_name(44:77), 'Interpreter', 'none')
    axis tight
    drawnow    
end %%%%%%%%%%%%%%%%%%%%%%%
roi_perc = [roi_perc(1:9), roi_perc(10:18)];
roi1 = [roi1(1:9), roi1(10:18)];
roi2 = [roi2(1:9), roi2(10:18)];
%%
% temp = regionprops3(yim, 'Image', 'Centroid')
names = {'TR' 'RET' 'HC' 'RET' 'CON' 'HC' 'CON' 'HC' 'CON'};
figure(38); clf; hold on
% plot([0 9], [0 0], 'k', 'LineWidth', 1)
% plot(2:4, roi_perc(10:12), '-k', 'LineWidth', 2)
% plot(5:7, roi_perc(13:15), '-k', 'LineWidth', 2)
% plot(1:8, roi_perc(9:16), 'o--k')
% plot(2:4, roi_perc(2:4), '-b', 'LineWidth', 2)
% plot(5:7, roi_perc(5:7), '-b', 'LineWidth', 2)
% plot(1:8, roi_perc(1:8), 'o--b')
% plot([0 10], [0 0], 'k', 'LineWidth', 1)
plot(1:9, roi_perc(:,1), '-bo', 'LineWidth', 2)
plot(1:9, roi_perc(:,2), '-ko', 'LineWidth', 2)
ylabel('Prop. frames with ROI detected')
xlabel('Session (connected are same day)')
title('\color{black}mHPC23459, \color{blue}mHPC23454')

axis([ .5 9.5 -.1 .6])
set(gca, 'XTick', 1:9, 'XTickLabels', names)


% d1 = ( roi1(:,1)-roi2(:,1) ) ./ (roi1(:,1)+roi2(:,1));
% d2 = ( roi1(:,2)-roi2(:,2) ) ./ (roi1(:,2)+roi2(:,2));
% figure(39); clf; hold on
% plot([0 10], [0 0], 'k', 'LineWidth', 1)
% plot(1:9, d1, '-bo', 'LineWidth', 2)
% plot(1:9, d2, '-ko', 'LineWidth', 2)
% ylabel('Prop. frames with ROI detected')
% xlabel('Session (connected are same day)')
% axis([ .5 9.5 -1.2 1.2])
% title('\color{black}mHPC23459, \color{blue}mHPC23454')
% set(gca, 'XTick', 1:9, 'XTickLabels', names)


%%
% Need to refine foci identification to remove jumps and keep identity
% consistent
[nt, nd] = size(cx);
reach = 2;

for d = 1:nd
    cxdiff = abs(diff(cx(:,d)));
    cydiff = abs(diff(cy(:,d)));
    adiff  = abs(diff(area(:,d)));
    pdiff  = abs(diff(pm(:,d)));
    
    
end

for ind = 2:nt
    num_regions = sum(~isnan(cx(ind+1, :)));
    prev_regions = sum(~isnan(cx(ind-reach, :)));
    if num_regions < prev_regions && (num_regions>1 || prev_regions>1) % sum(  any(~isnan(cx([ind-1:ind],:)), 2) ) ==2
    dx = []; dy = []; da = []; dp = []; f = [];
    for ttt = 1:prev_regions
        for rrr = 1:prev_regions
%             f(ttt, rrr) = area(ttt,rrr);
            dx(ttt, rrr) = cx(ind-reach, ttt)   - cx(ind+1, rrr);
            dy(ttt, rrr) = cy(ind-reach, ttt)   - cy(ind+1, rrr);
            da(ttt, rrr) = area(ind-reach, ttt) - area(ind+1, rrr);
            dp(ttt, rrr) = pm(ind-reach, ttt)   - pm(ind+1, rrr);
        end
    end
    dd = sqrt(dx.^2 + dy.^2 + da.^2 + dp.^2);
    labels = 1:prev_regions;
    for rrr = 1:num_regions
            best = find(dd(:,rrr)==nanmin(dd(:,rrr)));
            if best~=rrr
                
                figure;
                imagesc(area(ind-10:ind+10,:))
            end
%             bestmatch_idx(best) = rrr;
    end
    end
end

%%
pqwerty
nc = 10;
[p.coeff, p.score, p.latent, p.tsquared, p.explained, p.mu] = pca(single(yss)', ...
    'Algorithm', 'als', 'NumComponents', nc, 'Centered', 'off');

s = p.score';
fp   = zeros(h, w, nc);
% lat  = reshape(p.latent, [h w]);
% exp  = reshape(p.explained, [h w]);
mu   = reshape(p.mu, [h w]);

for i = 1:nc
    fp(:,:,i)  = reshape(p.coeff(:,i), [h w]);
end

figure(nc); clf
for i = 1:nc
    subplot_tight(ceil(sqrt(nc)),ceil(sqrt(nc)), i)
    imagesc(squeeze(fp(:,:,i))); 
    title(i)
    axis image off
    drawnow; 
end

figure; 
stacked_traces(normalize_rows(s))
%%
cmap = jet(nc);
cmap = cmap(randperm(nc),:);
nf2 = size(s,2);
y2 = p.coeff*p.score';
aa = min(min(min(y2)));
bb = max(max(max(y2)));


yr = p.coeff*( p.score.*( ones(nf2 ,1)*cmap(:,1)' ) )';
yg = p.coeff*( p.score.*( ones(nf2 ,1)*cmap(:,2)' ) )';
yb = p.coeff*( p.score.*( ones(nf2 ,1)*cmap(:,3)' ) )';

yimr = reshape(yr, [h, w, size(yr,2)]);
yimg = reshape(yg, [h, w, size(yr,2)]);
yimb = reshape(yg, [h, w, size(yr,2)]);

figure(2); clf; 
for i = 1:4:nf2
    subplot(1,2,1)
    im = squeeze(yim(:,:,i));
    imagesc(im, [20 50]); 
    subplot(1,2,2)
    im = cat(3, yimr(:,:,i), yimg(:,:,i), yimb(:,:,i));
    im = im - aa / (bb-aa);
    image(im); 
    title(i)
    drawnow; 
end

%% Examples showing ca wave spreading depression for grant

tiffname1 = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_12\16_55_19_RET8\HPC_miniscope1\msCam_MC.tiff';
tiffname2 = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_12\16_55_19_RET8\HPC_miniscope1\msCam_MC_cawaves.tiff';
temp = readtable('E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_12\16_55_19_RET8\HPC_miniscope1\timestamps.csv');
ca = load('E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_12\16_55_19_RET8\HPC_miniscope1\ca_waves.mat');
t = temp.TimeStamp_ms_./1000;
%%
startf = 22130; % 22142;
t0 = t(startf);
tstep = 2;
tind = [startf, find(t>=t0+1*tstep,1), find(t>=t0+2*tstep,1), find(t>=t0+3*tstep,1)];
tall = [t(tind)]-t0;
x = ca.cx(tind);
y = ca.cy(tind);
%%
zooms = [90, 90+120, 96, 96+120];
figure(1); clf; colormap bone
set(gcf, 'Color', 'w', 'Position', [210         609        1050         369])
for i = 1:4
    subplot_tight(2,5,i+1)
    im = imread(tiffname1, tind(i));
%     im = im(zooms(1):zooms(2), zooms(3):zooms(4));
    imagesc(im, [0 255]); hold on
%     plot(x(1:i)*4, y(1:i)*4, 'r.-')
    title(sprintf('%d sec', round(tall(i))));
    axis image off
    axis([zooms]);
    
    subplot_tight(2,5,i+6)
    im = imread(tiffname2, tind(i));
%     im = imresize(im,2);
%     im = im(zooms(1):zooms(2), zooms(3):zooms(4));
    imagesc(im, [0 40]); hold on
%     plot(x(1:i)*2, y(1:i)*2, 'r.-')
    title(sprintf('%d sec', round(tall(i))));
    axis image off
    axis([round(zooms/2)]);
end
%     subplot_tight(2,100,200)
% colorbar

subplot_tight(2,5,1)
im = imread(tiffname1, tind(1));
imagesc(im, [0 255]); hold on
% rectangle('Position', [zooms(1) zooms(3), zooms(2)-zooms(1), zooms(4)-zooms(3)], 'EdgeColor', 'r')
title(sprintf('Full FOV, %d sec', round(tall(1))));
axis image off

subplot_tight(2,5,6)
im = imread(tiffname2, tind(1));
imagesc(im, [0 40]); hold on
% rectangle('Position', round([zooms(1) zooms(3), zooms(2)-zooms(1), zooms(4)-zooms(3)]/2), 'EdgeColor', 'r')
title(sprintf('SD filter, %d sec', round(tall(1))));
axis image off











