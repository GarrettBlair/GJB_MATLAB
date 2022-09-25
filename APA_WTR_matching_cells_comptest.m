figure(5); clf
all_animal_name = {'Hipp18240'};%animals{animalLoop};

pop_sec_res = 60;
cell_sec_res = 3;
pall = cell(2,1);

rng('default')
for aLoop = 1:2
animal_name = all_animal_name{aLoop};%animals{animalLoop};
experiment_folder = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\';
dir_list_fname = '_directory_list.csv';

dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
processedDir        = sprintf('%s%s/processed_files/', experiment_folder, animal_name);
contourDir          = sprintf('%s%s/matching_contours/', experiment_folder, animal_name);
matchingfname = sprintf('%smatching_matrix.mat', contourDir);
DAT_Dir             = sprintf('%sDAT_files/', experiment_folder);

AnimalDir = setup_imaging_Sessionfiles(animal_name, dir_file, DAT_Dir, processedDir, contourDir);
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
for s1 = 2:nsess
    for s2 = s1:nsess
        if s1~=s2
            fname1 = AnimalDir.processedFile{s1};
            fname2 = AnimalDir.processedFile{s2};
            
            f1 = load(fname1);
            f2 = load(fname2);
            
%           matched = sum(cellmap>0,2)==size(cellmap,2); %
            matched = cellmap(:,s1)>0 & cellmap(:,s2)>0;
%             sum(matched)
            spks1 = f1.ms.neuron.S_matw(cellmap(matched,s1), :);
            spks2 = f2.ms.neuron.S_matw(cellmap(matched,s2), :);
            
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
            
            [popcorr1, ~, ~] = Fenton_pop_stability(spks1, pop_sec_res, f1.ms.timestamps./1000, false);
            [popcorr2, ~, ~] = Fenton_pop_stability(spks2, pop_sec_res, f2.ms.timestamps./1000, false);
            %             [popcorr3, ~, cellcorr3, ~, ~, ~] = Fenton_pop_stability(cat(2, spks1, spks2), pop_dt1, false);
            
            [cellcorr1, ~, ~] = Fenton_cell_corr(spks1, pcell_dt1, false);
            [cellcorr2, ~, ~] = Fenton_cell_corr(spks2, pcell_dt2, false);
            c1 = triu(cellcorr1,1);
            c2 = triu(cellcorr2,1);
            lowerinds = c1==0 & c2==0;
            sess_corr(s1,s2) = corr(c1(~lowerinds), c2(~lowerinds));
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
% fdir = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\datfiles\V4Min\';
% afname = sprintf('%s%s', 'GRAU_V4Miniscope_TestMouse2_CaImaging_Pretraining_Arena.dat');
% rfname = sprintf('%s%s', 'GRAU_V4Miniscope_TestMouse2_CaImaging_Pretraining_Room.dat');
clear
behav = [];
fdir = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\datfiles\V4Min\';
% afname = sprintf('%s%s', fdir, 'GRAU_V4Miniscope_TestMouse1_CaImaging_Pretraining_Arena.dat');
% rfname = sprintf('%s%s', fdir, 'GRAU_V4Miniscope_TestMouse1_CaImaging_Pretraining_Room.dat');
% behav.msFileName = "C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\TestMouse1\processed_files\2022_07_05___17_09_41_ms_placecells_data.mat";
% afname = sprintf('%s%s', fdir, 'GRAU_V4Miniscope_TestMouse1_CaImaging_TR2-1_Arena.dat');
% rfname = sprintf('%s%s', fdir, 'GRAU_V4Miniscope_TestMouse1_CaImaging_TR2-1_Room.dat');
% behav.msFileName = "C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\TestMouse1\processed_files\2022_07_13___15_12_07_ms_placecells_data.mat";
afname = sprintf('%s%s', fdir, 'GRAU_V4Miniscope_TestMouse1_CaImaging_TR2-RET_Arena.dat');
rfname = sprintf('%s%s', fdir, 'GRAU_V4Miniscope_TestMouse1_CaImaging_TR2-RET_Room.dat');
behav.msFileName = "C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\TestMouse1\processed_files\2022_07_18___09_27_18_ms_placecells_data.mat";

% afname = sprintf('%s%s', fdir, 'GRAU_V4Miniscope_TestMouse2_CaImaging_RET_Arena.dat');
% rfname = sprintf('%s%s', fdir, 'GRAU_V4Miniscope_TestMouse2_CaImaging_RET_Room.dat');
% behav.msFileName = "C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\TestMouse2\processed_files\2022_07_18___08_51_01_ms_placecells_data.mat";
% afname = sprintf('%s%s', fdir, 'GRAU_V4Miniscope_TestMouse2_CaImaging_TR1_Arena.dat');
% rfname = sprintf('%s%s', fdir, 'GRAU_V4Miniscope_TestMouse2_CaImaging_TR1_Room.dat');
% behav.msFileName = "C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\TestMouse2\processed_files\2022_07_12___11_15_53_ms_placecells_data.mat";

behav.roomFname = rfname;
behav.arenaFname = afname;
%%
[arena, arena_params] = read_APA_csv(afname, []);
[room , room_params] = read_APA_csv(rfname, []);

behav.timestamps = arena.timestamps;
behav_params.arena_radius = 20;
% % % [xx,yy] = ginput(3);
% % % [R,xcyc] = fit_circle_through_3_points([xx yy])
behav_params.center = [130 130];
behav_params.radius = 112;
behav_params.behav_fps = round(1/median(abs(diff(behav.timestamps./1000))));
behav_params.behav_smoothing_interval = .5;
behav_params.num_partitions           = 2;
behav_params.max_spd_thresh           = 50;
behav_params.min_spd_thresh           = 2;
behav_params.min_samples              = 10;

% ARENA FRAME
nanind = arena.rawx==0 & arena.rawy == 0;
x = arena.rawx;
y = arena.rawy;
x(nanind) = NaN;
y(nanind) = NaN;
t = arena.timestamps;
if length(unique(t)) ~= length(t)
    ind = find(diff(t)==0)+1;
    t(ind) = t(ind)+1;
end
ts = t./1000; % convert from ms to seconds
nanind = (isnan(x) & isnan(y));
xn = interp1(ts(~nanind), x(~nanind), ts(nanind), 'linear');
yn = interp1(ts(~nanind), y(~nanind), ts(nanind), 'linear');
x(nanind) = xn; y(nanind) = yn;

x = behav_params.arena_radius*(x-behav_params.center(1))./behav_params.radius;
y = behav_params.arena_radius*(y-behav_params.center(2))./behav_params.radius;
ksize = round(behav_params.behav_fps*behav_params.behav_smoothing_interval);
kern = ones(ksize, 1); kern = kern./sum(kern(:));
arena.x = conv(x, kern, 'same');
arena.y = conv(y, kern, 'same');

% ROOM FRAME
nanind = room.rawx==0 & room.rawy == 0;
x = room.rawx;
y = room.rawy;
x(nanind) = NaN;
y(nanind) = NaN;
t = room.timestamps;
if length(unique(t)) ~= length(t)
    ind = find(diff(t)==0)+1;
    t(ind) = t(ind)+1;
end
ts = t./1000; % convert from ms to seconds
nanind = (isnan(x) & isnan(y));
xn = interp1(ts(~nanind), x(~nanind), ts(nanind), 'linear');
yn = interp1(ts(~nanind), y(~nanind), ts(nanind), 'linear');
x(nanind) = xn; y(nanind) = yn;

x = behav_params.arena_radius*(x-behav_params.center(1))./behav_params.radius;
y = behav_params.arena_radius*(y-behav_params.center(2))./behav_params.radius;
ksize = round(behav_params.behav_fps*behav_params.behav_smoothing_interval);
kern = ones(ksize, 1); kern = kern./sum(kern(:));
room.x = conv(x, kern, 'same');
room.y = conv(y, kern, 'same');

%
behav_dt = [median(diff(behav.timestamps)); diff([behav.timestamps])]/1000;
behav.dt = behav_dt;
room.speed  =  sqrt(diff([room.x(1); room.x]).^2   + diff([room.y(1); room.y]).^2)./behav_dt;
arena.speed =  sqrt(diff([arena.x(1); arena.x]).^2 + diff([arena.y(1); arena.y]).^2)./behav_dt;

behav.room = room; behav.arena = arena;


temp = load(behav.msFileName);
ms = temp.ms;
params = temp.params;

ms.room.x = interp1(behav.timestamps,  behav.room.x,  ms.timestamps, 'linear');
ms.room.y = interp1(behav.timestamps,  behav.room.y,  ms.timestamps, 'linear');
ms.arena.x = interp1(behav.timestamps, behav.arena.x, ms.timestamps, 'linear');
ms.arena.y = interp1(behav.timestamps, behav.arena.y, ms.timestamps, 'linear');

dt = ms.dt_corrected;
ms.room.speed = sqrt(diff([ms.room.x(1); ms.room.x]).^2 + diff([ms.room.y(1); ms.room.y]).^2)./dt;
ms.arena.speed = sqrt(diff([ms.arena.x(1); ms.arena.x]).^2 + diff([ms.arena.y(1); ms.arena.y]).^2)./dt;
ms.room.speed_smooth = conv(ms.room.speed,   kern, 'same');
ms.arena.speed_smooth = conv(ms.arena.speed, kern, 'same');
% is_moving = ms.arena.speed_smooth>params_sub.min_spd_thresh;

[~, speed_epochs] = get_speed_epochs(ms.arena.speed_smooth, behav_params);

ms.speed_epochs = speed_epochs;
is_moving = speed_epochs;

%
params.pos_bins = [-20:2:20];
params.occupancy_thresh = .00; % in seconds
N = params.pfield_kernel_radius*2+1; alpha = 2.5; % default
params.pfield_gauss_std = ((N-1)/(2*alpha))*abs(params.pos_bins(1)-params.pos_bins(2));

spks = normalize_rows(ms.neuron.S_matw);
[ms.room]   = construct_place_maps_2D(ms.room,  ms.room.x(is_moving),  ms.room.y(is_moving),  ms.dt(is_moving), spks(:, is_moving), params.pos_bins, params);
[ms.arena]  = construct_place_maps_2D(ms.arena, ms.arena.x(is_moving), ms.arena.y(is_moving), ms.dt(is_moving), spks(:, is_moving), params.pos_bins, params);
%
[ms.arena.pcell_stats] = place_cell_stats(spks, ms.arena.pfields, ms.arena.spkmap, ms.arena.vmap);
[ms.room.pcell_stats] = place_cell_stats(spks, ms.room.pfields,  ms.room.spkmap,  ms.room.vmap);

ms.arena.pfield_alpha = ~isnan(ms.arena.vmap);
ms.room.pfield_alpha = ~isnan(ms.room.vmap);
%%
rng(100);
nsegs = size(spks,1);
segs = 1:nsegs; 
ds = 30;
[~, ord] = sort(rand(1,nsegs));
segs2use = segs(ord(1:ds));
% [~, spord] = sort(ms.arena.pcell_stats.spkRate, 'descend');
% segs2use = segs(spord(1:ds));

% as = ms.arena.pcell_stats.infoPerSpike(segs2use);
infspk = ms.arena.pcell_stats.infoPerSpike(segs2use);
% spkr = ms.arena.pcell_stats.spkRate(segs2use(1:ds));
as = infspk;
[~, arena_info_ord] = sort(as, 'descend');
p_a = ms.arena.pfields_smooth(segs2use,:,:);
p_a = p_a(arena_info_ord,:,:);

s_sub = spks(segs2use,:);
s_sub = s_sub(arena_info_ord,:);
% c_sub = ms.neuron.C(segs2use,:);
c_sub = ms.neuron.C(segs2use,:) + ms.neuron.YrA(segs2use,:);
c_sub = c_sub(arena_info_ord,:);
c_sub = normalize_rows(c_sub);
p_r = ms.room.pfields_smooth(segs2use,:,:);
p_r = p_r(arena_info_ord,:,:);
% ss = get(0,'screensize');
%     ar = mean([1, ss(3)/ss(4)])-1;
%     n = ceil(sqrt(nsegs/ds));
    nc = 2; % round(n*(1-ar))+1;
    nr = ds; % round(n*(1+ar));
    
    smoothing_kern = gausswin(7);
% smoothing_kern = gausswin(ksize+1);
smoothing_kern = smoothing_kern*smoothing_kern';
smoothing_kern = smoothing_kern./(sum(smoothing_kern(:)));

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
    