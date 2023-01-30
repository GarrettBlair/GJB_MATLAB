figure(5); clf
all_animal_name = {'Hipp18240'};%animals{animalLoop};

pop_sec_res = 60;
cell_sec_res = 3;
pall = cell(2,1);

rng('default')
for aLoop = 1:length(all_animal_name)%:2
animal_name = all_animal_name{aLoop};%animals{animalLoop};
experiment_folder = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\';
dir_list_fname = '_directory_list.csv';

dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
% processedDir        = sprintf('%s%s/processed_files/', experiment_folder, animal_name);
contourDir          = sprintf('%s%s/matching_contours/', experiment_folder, animal_name);
matchingfname       = sprintf('%smatching_matrix.mat', contourDir);
DAT_Dir             = sprintf('%sDAT_files/', experiment_folder);

AnimalDir = setup_imaging_Sessionfiles(animal_name, dir_file, experiment_folder);%, DAT_Dir, processedDir, contourDir);
nsess = AnimalDir.numSess;
load(matchingfname);

sess_corr1 = NaN(nsess);
% subplot(1,3,1)
colors = viridis(nsess*2);
colors = colors(floor(nsess/2):floor(nsess/2)+nsess-1, :);
light_colors = colors + (1-colors)/3;

numcells2use = [];
sess_corr = NaN(nsess);
pall{aLoop} = NaN(nsess,30);
for s1 = 1:nsess
    for s2 = s1+1:nsess
        %%
        if s1~=s2
            %%
            fname1 = AnimalDir.processedFile{s1};
            cname1 = AnimalDir.contourFile{s1};
            fname2 = AnimalDir.processedFile{s2};
            
            f1 = load(fname1);
            c1 = load(cname1);
            f2 = load(fname2);
            im = zeros(size(c1.contours,2), size(c1.contours,3), 3);
            for i = 1:size(c1.contours,1)
                r = squeeze(c1.contours(i,:,:)).*rand(1);
                g = squeeze(c1.contours(i,:,:)).*rand(1);
                b = squeeze(c1.contours(i,:,:)).*rand(1);
                
                im = im + cat(3, r,g,b);
            end
            
            matched = cellmap(:,s1)>0 & cellmap(:,s2)>0;
            cells_s1 = cellmap(matched,s1);
            cells_s2 = cellmap(matched,s2);
            spks1 = f1.ms.neuron.S_matw(cells_s1, :);
            spks2 = f2.ms.neuron.S_matw(cells_s2, :);
            
            spks1 = normalize_rows(spks1);
            spks2 = normalize_rows(spks2);
            nsamples = min(size(spks1,2), size(spks2,2));
            spks1 = spks1(:, 1:nsamples);
            spks2 = spks2(:, 1:nsamples);
            
            if ~isempty(numcells2use)
                [~, randord] = sort(rand(sum(matched), 1));
                spks1 = spks1(randord(1:numcells2use), :);
                spks2 = spks2(randord(1:numcells2use), :);
            end
            pop_dt1 = round(1/median(f1.ms.dt))*pop_sec_res;
            pop_dt2 = round(1/median(f2.ms.dt))*pop_sec_res;
            pcell_dt1 = round(1/median(f1.ms.dt))*cell_sec_res;
            pcell_dt2 = round(1/median(f2.ms.dt))*cell_sec_res;
            
            [timecorr1, ~, ~] = Fenton_pop_stability(spks1, pop_sec_res, f1.ms.timestamps(1:nsamples)./1000, false);
            [timecorr2, ~, ~] = Fenton_pop_stability(spks2, pop_sec_res, f2.ms.timestamps(1:nsamples)./1000, false);
            %             [popcorr3, ~, cellcorr3, ~, ~, ~] = Fenton_pop_stability(cat(2, spks1, spks2), pop_dt1, false);
            
            [cellcorr1, ~, ~] = Fenton_cell_corr(spks1, cell_sec_res, f1.ms.timestamps(1:nsamples)./1000, false);
            [cellcorr2, ~, ~] = Fenton_cell_corr(spks2, cell_sec_res, f2.ms.timestamps(1:nsamples)./1000, false);
            
            c1 = triu(cellcorr1,1);
            c2 = triu(cellcorr2,1);
            lowerinds = c1==0 & c2==0;
            sess_corr(s1,s2) = corr(c1(~lowerinds), c2(~lowerinds));
            
            %% place fields
            p1 = f1.ms.room.pfields_smooth(cells_s1, :, :);
            p2 = f2.ms.room.pfields_smooth(cells_s2, :, :);
            room_betweencorr = NaN(sum(matched), 1);
            for ploop = 1:sum(matched)
                p1sub = squeeze(p1(ploop,:,:));
                p2sub = squeeze(p2(ploop,:,:));
                valid_inds = ~isnan(p1sub.*p2sub);
                room_betweencorr(ploop) = corr(p1sub(valid_inds), p2sub(valid_inds));
            end
            p1 = f1.ms.arena.pfields_smooth(cells_s1, :, :);
            p2 = f2.ms.arena.pfields_smooth(cells_s2, :, :);
            arena_betweencorr = NaN(sum(matched), 1);
            for ploop = 1:sum(matched)
                p1sub = squeeze(p1(ploop,:,:));
                p2sub = squeeze(p2(ploop,:,:));
                valid_inds = ~isnan(p1sub.*p2sub);
                arena_betweencorr(ploop) = corr(p1sub(valid_inds), p2sub(valid_inds));
            end
            
        else
            fname1 = AnimalDir.processedFile{s1};
            f1 = load(fname1);
            
            matched = cellmap(:,s1)>0;
            spks1 = f1.ms.neuron.S_matw(cellmap(matched,s1), :);
            
            spks1 = normalize_rows(spks1);
            nsamples = size(spks1,2);
            spks1 = spks1(:, 1:nsamples);
            
            if ~isempty(numcells2use)
                [~, randord] = sort(rand(sum(matched), 1));
                spks1 = spks1(randord(1:numcells2use), :);
            end
            pop_dt1 = round(1/median(f1.ms.dt))*pop_sec_res;
            pcell_dt1 = round(1/median(f1.ms.dt))*cell_sec_res;
            
            [popcorr1, ~, ~] = Fenton_pop_stability(spks1, pop_sec_res, f1.ms.timestamps./1000, false);
            popcorr1_within = popcorr1;
            figure; 
            set(gcf, 'Name', sprintf('Animal - %s, s1=%d', animal_name, s1))
            imagesc(popcorr1, [-.2 1]);
            colorbar
        end
        p1 = nanmean(popcorr1_within,1);
        pall{aLoop}(s1,:) = p1;
        %             p2 = nanmean(popcorr2,1);
        %             p3 = nanmean(popcorr3,1);
        nt = length(p1);
        %             subplot(5,5, 5*(s1-1)+s2)
        figure(5);
        subplot(1,2,1)
        hold on
        plot(1:nt, p1, 'Color', colors(s1,:))
        %             plot(nt+1:2*nt, p2, 'Color', light_colors(s2,:))
        drawnow
    end
end
%             imagesc(popcorr3, [-.1, 1])
%             imagesc(cellcorr3, [-.1, 1])
if strcmp(animal_name, 'TestMouse1')
sess_corr1 = sess_corr;
else
sess_corr2 = sess_corr;
end
% plot(1:nt, nanmean(pall{aLoop},1), 'Color', mean(colors,1), 'LineWidth', 3)
end
subplot(1,2,1)
axis([-1, nt+1, -.1, 1.1])
xv = [0:15:30]; 
set(gca, 'XTick', xv, 'XTickLabel', xv)
xlabel('Session Time (min)')
ylabel('Population stability (corr)')

%%
% for i = 1:2

a = sess_corr1(~isnan(sess_corr1)); b = sess_corr2(~isnan(sess_corr2));
subplot(1,2,2); cla
hold on

aa = mean(a); bb = mean(b);

colors = viridis(20);
colors = colors(6:10,:);
c1 = mean(colors,1);
colors = plasma(20);
colors = colors(13:17,:);
c2 = mean(colors,1);

bar(0, aa, 'FaceColor', c1/1.5);
scatter(zeros(length(a),1), a, 'MarkerFaceColor', c1, 'MarkerEdgeColor', 'k');

bar(1, bb, 'FaceColor', c2/1.5);
scatter( ones(length(b),1), b, 'MarkerFaceColor', c2, 'MarkerEdgeColor', 'k')
set(gca, 'XTickLabel', {'Mouse1', 'Mouse2'}, 'XTick', [0 1])
axis([-1.5 2.5 0 1.1])
axis square
xlabel('')
ylabel('PCO (corr)')
%%



%%
matchingfile = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\matching_contours\matching_matrix.mat';
ddir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\';
filenames = ...
   {'2022_09_16_H17_03_07_TR9';...
    '2022_09_16_H17_41_09_WTR10';...
    '2022_09_16_H18_19_23_TR11'};
% filenames = ...
%    {'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\2022_09_14_H17_07_45_TR7_@placecells.mat';...
%     'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\2022_09_14_H17_52_22_WTR8_@placecells.mat';...
%     'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\2022_09_16_H17_03_07_TR9_@placecells.mat';...
%     'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\2022_09_16_H17_41_09_WTR10_@placecells.mat';...
%     'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\2022_09_16_H18_19_23_TR11_@placecells.mat';...
%     'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\2022_09_24_H16_11_01_WTR19_@placecells.mat';...
%     'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\2022_09_24_H16_46_19_TR20_@placecells.mat';...
%     'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\2022_09_27_H17_23_57_WTR21_@placecells.mat';...
%     'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\2022_09_27_H17_58_42_TR22_@placecells.mat';...
%     'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\2022_09_27_H18_26_28_TR23_@placecells.mat'};
pop_sec_res = 30;
cell_sec_res = 5;
spd_thresh = 5; %cm/sec
numRand = 10;
figdir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\WTR_manip\processed_files\temp_output\';
temp = load(matchingfile);
% cellmap = temp.cellmap;
cellmap = temp.cellmap(:,3:5);
nsess = size(cellmap,2);
sess_corr = NaN(nsess);
sess_absdiff = NaN(nsess);
time_corr_mean = NaN(nsess,1);
time_corr_std = NaN(nsess,1);

num_shared = [cellmap'>0] * [cellmap>0];
min_shared = min(num_shared(:));
numcells2use = [];%
if isempty(numcells2use)
    numRand = 1;
end
figure(95); 
clf
room_pfieldcorr  = cell(nsess);
arena_pfieldcorr = cell(nsess);
for s1 = 1:nsess
    for s2 = s1:nsess
        %%
        fname1 = sprintf('%sprocessed_files\\%s_@placecells.mat', ddir, filenames{s1});
        cname1 = sprintf('%smatching_contours\\%s_contours.mat', ddir, filenames{s1});
        temp = load(fname1, 'ms');
        ms1 = temp.ms;
        temp = load(cname1, 'contours');
        conts1 = temp.contours;
        
        if s1==s2
%             single_sess_anaylsis
            spks1 = ms1.neuron.S_matw;
            spks1 = normalize_rows(spks1);
            [timecorr1, ~, ~] = Fenton_pop_stability(spks1, pop_sec_res, ms1.timestamps./1000, false);
            time_corr_mean(s1) = nanmean(timecorr1(:));
            time_corr_std(s1) = nanstd(timecorr1(:))*2;
        else
%             cross_sess_anaylsis
        fname2 = sprintf('%sprocessed_files\\%s_@placecells.mat', ddir, filenames{s2});
        cname2 = sprintf('%smatching_contours\\%s_contours.mat', ddir, filenames{s2});
        temp = load(fname2, 'ms');
        ms2 = temp.ms;
        temp = load(cname2, 'contours');
        conts2 = temp.contours;
        
        
        matched = cellmap(:,s1)>0 & cellmap(:,s2)>0;
        matched = sum(cellmap>0,2)==nsess; % use cell matched across all sessions
        cells_s1 = cellmap(matched,s1);
        cells_s2 = cellmap(matched,s2);
        conts1 = conts1(cells_s1,:,:);
        conts2 = conts2(cells_s2,:,:);
        spks1 = ms1.neuron.S_matw(cells_s1, :);
        spks2 = ms2.neuron.S_matw(cells_s2, :);
        
        pfields1_room = ms1.room.pfields_smooth(cells_s1, :, :);
        pfields2_room = ms2.room.pfields_smooth(cells_s2, :, :);
        pfields1_arena = ms1.arena.pfields_smooth(cells_s1, :, :);
        pfields2_arena = ms2.arena.pfields_smooth(cells_s2, :, :);
        %%
        for pfieldLoop = 1:sum(matched)
            %%
%             room_alpha  = ms1.room.pfield_alpha.*ms2.room.pfield_alpha;
%             arena_alpha = ms1.arena.pfield_alpha.*ms2.arena.pfield_alpha;
            p1 = squeeze(pfields1_room(pfieldLoop,:,:));
            p2 = squeeze(pfields2_room(pfieldLoop,:,:));
            validinds = ~isnan(p1.*p2);
            cc = corr(p1(validinds), p2(validinds));
            figure([s1*100+s2]); clf;
            subplot(2,2,1)
            imagesc(p1, 'AlphaData', ms1.room.pfield_alpha); axis image off
            title(sprintf('Sess%d Room, %1.3f', s1, cc))
            subplot(2,2,2)
            imagesc(p2, 'AlphaData', ms2.room.pfield_alpha); axis image off
            title(sprintf('Sess%d Room, %1.3f', s2, cc))
            room_pfieldcorr{s1, s2} = cat(1, room_pfieldcorr{s1, s2}, cc);

            p1 = squeeze(pfields1_arena(pfieldLoop,:,:));
            p2 = squeeze(pfields2_arena(pfieldLoop,:,:));
            validinds = ~isnan(p1.*p2);
            cc = corr(p1(validinds), p2(validinds));
            subplot(2,2,3)
            imagesc(p1, 'AlphaData', ms1.arena.pfield_alpha); axis image off
            title(sprintf('Sess%d Arena, %1.3f', s1, cc))
            subplot(2,2,4)
            imagesc(p2, 'AlphaData', ms2.arena.pfield_alpha); axis image off
            title(sprintf('Sess%d Arena, %1.3f', s2, cc))
            arena_pfieldcorr{s1, s2} = cat(1, arena_pfieldcorr{s1, s2}, cc);
            temp = getframe(gcf); 
            imwrite(temp.cdata, sprintf('%s\\sess%d-sess%d_seg%d_Arena.png', figdir, s1, s2, pfieldLoop))
        end
        figure(95); 
        subplot(nsess, nsess, sub2ind([nsess, nsess], s1, s2)); hold on
        violinplot([room_pfieldcorr{s1,s2}, arena_pfieldcorr{s1,s2}])

        spks1 = normalize_rows(spks1);
        spks2 = normalize_rows(spks2);
        nsamples = min(size(spks1,2), size(spks2,2)); % minimum samles
        spks1 = spks1(:, 1:nsamples);
        spks2 = spks2(:, 1:nsamples);
        
        pop_dt1 = round(1/median(ms1.dt))*pop_sec_res;
        pop_dt2 = round(1/median(ms2.dt))*pop_sec_res;
        pcell_dt1 = round(1/median(ms1.dt))*cell_sec_res;
        pcell_dt2 = round(1/median(ms2.dt))*cell_sec_res;
        
        spd1 = ms1.arena.speed_smooth;
        spd2 = ms2.arena.speed_smooth;
        c = NaN(numRand,1);
        d = NaN(numRand,1);
        for randloop = 1:numRand
            if ~isempty(numcells2use)
                [~, randord] = sort(rand(sum(matched), 1));
                spks1_sub = spks1(randord(1:numcells2use), :);
                spks2_sub = spks2(randord(1:numcells2use), :);
            elseif numcells2use<0 % negavtive for random
                [~, randord] = sort(rand(sum(matched), 1));
                spks1_sub = spks1(randord(1:numcells2use), :);
                [~, randord] = sort(rand(sum(matched), 1));
                spks2_sub = spks2(randord(1:numcells2use), :);
            else
                spks1_sub = spks1;
                spks2_sub = spks2;
            end
            spks1_sub(:,spd1<=spd_thresh) = 0;
            spks2_sub(:,spd2<=spd_thresh) = 0;
            spks1_sub = spks1_sub>0;
            spks2_sub = spks2_sub>0;
%             [timecorr1, ~, ~] = Fenton_pop_stability(spks1_sub, pop_sec_res, ms1.timestamps(1:nsamples)./1000, false);
%             [timecorr2, ~, ~] = Fenton_pop_stability(spks2_sub, pop_sec_res, ms2.timestamps(1:nsamples)./1000, false);
            %         joint_t = ms1.timestamps(1:nsamples)./1000;
            %         joint_t = [joint_t; joint_t(end) + cell_sec_res + ms2.timestamps(1:nsamples)./1000];
            %         [timecorr3, ~, ~] = Fenton_pop_stability([spks1, spks2], pop_sec_res, joint_t, false);
            %             [popcorr3, ~, cellcorr3, ~, ~, ~] = Fenton_pop_stability(cat(2, spks1, spks2), pop_dt1, false);
            
            [cellcorr1, ~, ~] = Fenton_cell_corr(spks1_sub, cell_sec_res, ms1.timestamps(1:nsamples)./1000, false);
            [cellcorr2, ~, ~] = Fenton_cell_corr(spks2_sub, cell_sec_res, ms2.timestamps(1:nsamples)./1000, false);
            
            c1 = triu(cellcorr1,1);
            c2 = triu(cellcorr2,1);
            lowerinds = c1==0 & c2==0;
            corrs1 = c1(~lowerinds);
            corrs2 = c2(~lowerinds);
            [corrs1_sort, ord] = sort(corrs1, 'descend');
            corrs2_sort = corrs2(ord);
            tops = floor(length(ord)/4);
            corrs2_sort = corrs2_sort(1:tops);
            corrs1_sort = corrs1_sort(1:tops);
            
            pw_diff = corrs1_sort - corrs2_sort;
            c(randloop) =  corr(corrs1_sort, corrs2_sort);
            d(randloop) =  median((pw_diff));
        end
        sess_corr(s1,s2) = nanmean(c);
        sess_corr(s2,s1) = sess_corr(s1,s2);
        
        sess_absdiff(s1,s2) = nanmean(d);
        sess_absdiff(s2,s1) = sess_absdiff(s1,s2);
        end
    end
end
figure(92); clf; subplot(121); imagesc(sess_corr, [0 1]); subplot(122); imagesc(sess_absdiff);
figure; subplot(131); imagesc(timecorr1, [-.2 .6]); subplot(132); imagesc(timecorr2, [-.2 .6]); subplot(133); imagesc(timecorr3, [-.2 .6]);
figure(94); clf; subplot(131); imagesc(cellcorr1, [-.3 .3]); subplot(132); imagesc(cellcorr2, [-.3 .3]); subplot(133); imagesc(cellcorr1-cellcorr2, [-.3 .3]);
% figure; hold on; plot(corrs1_sort); plot(corrs2_sort)

%%
    figure(10); clf;
    set(gcf, 'Name', 'ARENA FRAME / ROOM FRAME')
%     colormap(plasma)
%     figure(11); clf
%     set(gcf, 'Name', 'ROOM FRAME')
    colormap(viridis)
    figure(10);
    for i = 1:ds
        subplot_tight(2, ds+1, i, [.001 .001])
        p = squeeze(p_a(i,:,:));
        %     p = squeeze(ms.arena.pfields(i,:,:));
        p = p./max(p(:));
        ii = imagesc(p, [0 .5]);
        ii.AlphaData = ms.arena.pfield_alpha;
        axis square off
    end
%     figure(11);
    for i = 1:ds
        subplot_tight(2, ds+1, ds+i+1, [.001 .001])
        p = squeeze(p_r(i,:,:));
        %     p = squeeze(ms.room.pfields(i,:,:));
        p = p./max(p(:));
        ii = imagesc(p, [0 .5]);
        ii.AlphaData = ms.room.pfield_alpha;
%         ii = imagesc(p, 'AlphaData', ms.room.pfield_alpha);
%         ii.Parent.ALim = [0 .5];
        axis square off
    end
        p = p*NaN;
        subplot_tight(2, ds+1, [ds+1, 2*(ds+1)], [.1 .001]); cla
        ii = imagesc(p, [0 .5]);
        ii.AlphaData = NaN*ms.arena.pfield_alpha;
        colorbar
        axis off
    figure(11); clf;
    set(gcf, 'Name', 'ca traces (blk) and inferred spks (red)')
    hold on; 
    stacked_traces(c_sub, .9, {'k-'})
    stacked_traces(s_sub, .9, {'r-'})
    ts = ms.timestamps./1000;
    mts = mod(ts, 60.001);
    xd = find(mts(1:end-1) > mts(2:end));
    xl = abs(round(ts(xd)));
    set(gca, 'XTickLabel', xl, 'XTick', xd);
    xlabel('Time (sec)')
    