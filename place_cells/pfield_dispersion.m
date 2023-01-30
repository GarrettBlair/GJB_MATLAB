function [runs_z, alldiff_av] = pfield_dispersion(filename, plotting)
% function [expected_rate] = expected_firing_rate(ms_struct)
load(filename, 'ms')
% load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\2022_09_16_H17_41_09_WTR10_@placecells.mat')
% load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\2022_09_16_H17_03_07_TR9_@placecells.mat')
%% calculating expected spikes
bins = ms.params.pos_bins;
is_moving       = ms.speed_epochs;
spks = normalize_rows(ms.neuron.S_matw);
% spks = (ms.neuron.S_matw);
% spks = spks(: , is_moving);
spks(: , ~is_moving) = 0;

ms_frame = ms.room;
p = ms_frame.pfields_smooth;
% p = ms_frame.pfields;

nsegs = size(p,1);
numSamples = length(ms_frame.x);

[~, ~, ~, xbin, ybin] = histcounts2(ms_frame.x, ms_frame.y, bins, bins);
nanind = xbin==0 & ybin==0;
xbin(nanind) = interp1(find(~nanind), xbin(~nanind), find(nanind), 'nearest');
ybin(nanind) = interp1(find(~nanind), ybin(~nanind), find(nanind), 'nearest');

dt = ms.dt;
room_expected_FR = zeros(nsegs, length(ms_frame.x));
for segLoop = 1:nsegs
    p_e = squeeze(p(segLoop,:,:));
    for posLoop = 1:numSamples
        room_expected_FR(segLoop, posLoop) = p_e(ybin(posLoop), xbin(posLoop))*dt(posLoop);
    end
%     room_expected_FR(room_expected_FR<nanmean(p_e(:))) = NaN;
end
room_expected_FR(isnan(room_expected_FR)) = 0;
% room_expected_FR(: , ~is_moving) = 0;
% room_expected_FR = room_expected_FR(: , is_moving);

%
ms_frame = ms.arena;
p = ms_frame.pfields_smooth;
% p = ms_frame.pfields;

[~, ~, ~, xbin, ybin] = histcounts2(ms_frame.x, ms_frame.y, bins, bins);
nanind = xbin==0 & ybin==0;
xbin(nanind) = interp1(find(~nanind), xbin(~nanind), find(nanind), 'nearest');
ybin(nanind) = interp1(find(~nanind), ybin(~nanind), find(nanind), 'nearest');

arena_expected_FR = zeros(nsegs, length(ms_frame.x));
for segLoop = 1:nsegs
    p_e = squeeze(p(segLoop,:,:));
    for posLoop = 1:numSamples
        arena_expected_FR(segLoop, posLoop) = p_e(ybin(posLoop), xbin(posLoop))*dt(posLoop);
    end
%     arena_expected_FR(segLoop, arena_expected_FR(segLoop,:)<nanmean(p_e(:))) = NaN;
end
arena_expected_FR(isnan(arena_expected_FR)) = 0;
% arena_expected_FR(: , ~is_moving) = 0;
% arena_expected_FR = arena_expected_FR(: , is_moving);

% spike matrix, normalized to peak per cell
deltaT = 3;
t = ms.timestamps(is_moving)./1000;
% observed rate per cell
spks_bin_obs = bin_spks_time(spks, deltaT, ms.timestamps./1000, false); % t, false); % ms.timestamps./1000, false);
% average rate per cell
numBinSamples = size(spks_bin_obs, 2);
meanFR = sum(spks,2)./(ms.timestamps(end)./1000);
meanFRmat = meanFR*ones(1,numBinSamples);


room_expected_bin = bin_spks_time(room_expected_FR,  deltaT, ms.timestamps./1000, false);
bad_inds = (room_expected_bin-meanFRmat) <0;
room_dispersion = (spks_bin_obs - room_expected_bin)./sqrt(room_expected_bin);
room_dispersion(bad_inds) = NaN;%0;

arena_expected_bin = bin_spks_time(arena_expected_FR, deltaT, ms.timestamps./1000, false);
bad_inds = (arena_expected_bin-meanFRmat) < 0;
arena_dispersion = (spks_bin_obs - arena_expected_bin)./sqrt(arena_expected_bin);
arena_dispersion(bad_inds) = NaN;%0;
%
alldiff = nanmean(room_dispersion - arena_dispersion,1);
alldiff = alldiff./max(abs(alldiff));
alldiff_av = nanmedian(alldiff);
if plotting
scale =1;
figure(29); clf; 
seg = 22;%30;
ex_pfield = squeeze(ms.room.pfields_smooth(seg,:,:));
subplot_tight(2,4,1, [.01 .01]); hold on
imagesc(ex_pfield, [-.2, .3])
axis image off
colormap hot
colorbar

subplot_tight(2,4,2, [.1 .01]); hold on
plot(arena_expected_FR(seg,:), 'b')
plot(room_expected_FR(seg,:), 'r')
plot(spks(seg,:), 'k')
axis tight
set(gca, 'YTick', [])

subplot_tight(2,2,2, [.01 .075]); hold on
title('ROOM frame', 'Color', 'r')
stacked_traces(room_expected_bin, scale, {'r'})
stacked_traces(arena_expected_bin, scale, {'b'})
stacked_traces(spks_bin_obs, scale, {'k-'})
axis tight
ylabel('Cell #')
set(gca, 'XTick', [])


subplot_tight(2,2,4, [.1 .075]); cla; hold on
xm = nanmedian(alldiff);
plot([0 length(alldiff)], [xm xm], 'g', 'LineWidth', 2)
plot(alldiff, 'k-')
ylabel('Room - Arena')

axis tight
set(gca, 'XTick', 0:5*60/deltaT:size(room_dispersion,2), 'XTickLabel', (0:5*60/deltaT:size(room_dispersion,2))/20)
set(gca, 'YTick', -1:.5:1, 'YLim', [-1.2 1.2])
xlabel('Time (min)')

% histbins = [-.05:.005:.05];
histbins = [0:.25:4];

subplot_tight(2,4,5); hold on
% histogram(mean(room_dispersion, 2),  [-1.5:.075:1.5], 'FaceColor', 'r', 'FaceAlpha', .25)
histogram(nanvar(room_dispersion, [], 2),  histbins, 'FaceColor', 'r', 'FaceAlpha', .25)
xx = max(abs(get(gca, 'XTick')));
yy = max(abs(get(gca, 'YTick')));
% plot([0 0], [-1 20], 'k:', 'LineWidth', 3)
xm = mean(nanvar(room_dispersion, [], 2));
plot([xm xm], [-1 ceil(1.25*yy)], 'r', 'LineWidth', 2)
axis([-1, 1.5*xx, 0, 1.1*yy])


subplot_tight(2,4,6); hold on
% histogram(mean(arena_dispersion, 2), [-1.5:.075:1.5], 'FaceColor', 'b', 'FaceAlpha', .25)
histogram(nanvar(arena_dispersion, [], 2), histbins, 'FaceColor', 'b', 'FaceAlpha', .25)
xx = max(abs(get(gca, 'XTick')));
yy = max(abs(get(gca, 'YTick')));
% plot([0 0], [-1 20], 'k:', 'LineWidth', 3)
xm = mean(nanvar(arena_dispersion, [], 2));
plot([xm xm], [-1 ceil(1.25*yy)], 'b', 'LineWidth', 2)
axis([-1, 1.5*xx, 0, 1.1*yy])
end
runs_z=runtest(alldiff');
% % figure;
% % histbins = [-.05:.005:.05];
% % 
% % histogram2(mean(arena_dispersion, 2), mean(room_dispersion, 2), histbins, histbins, 'FaceColor', 'b', 'FaceAlpha', .25)
% % 
% % [theta1,rho1] = cart2pol(ms.arena.x,ms.arena.y);
