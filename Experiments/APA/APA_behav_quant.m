%% Read in the extracted position from APA_extract_pos_batch.ipynb
topdir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\2022_06_10\18_44_02';
% topdir = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/Hipp16942/2022_06_18/16_59_59';
fname_position = sprintf('%s/experiment/behav_position_data.csv', topdir);
fname_params = sprintf('%s/experiment/behav_ext_params.json', topdir);
caimanFilename = sprintf('%s/MiniLFOV/caiman_cnmfe_out.mat', topdir);
msTSFile = sprintf('%s/MiniLFOV/timeStamps.csv', topdir);
msCropFile = sprintf('%s/MiniLFOV/Crop_params.mat', topdir);
msOriFile = sprintf('%s/MiniLFOV/headOrientation.csv', topdir);
% Read in the coresponding files
[behav, params] = read_APA_csv(fname_position, fname_params);
caiman_data = load(caimanFilename);
TS_data = readtable(msTSFile);
crop_params = load(msCropFile);
ORI_Data = readtable(msOriFile);

%%
% plot3(behav.xc, behav.yc, behav.timestamps)
arena_radius = 40; % in cm
pos_bins = -40:4:40; % in cm, x and y
speed_thresh = 5; % speed thresh in cm/sec
smoothing_interval = .5; % in seconds, length of smoothing kernel


rotate_behav = true;
nan_interp = true;
nanind = behav.x==0 & behav.y == 0;
%
x = behav.x;
y = behav.y;
x(nanind) = NaN;
y(nanind) = NaN;
t = behav.timestamps;
if length(unique(t)) ~= length(t)
    ind = find(diff(t)==0)+1;
    t(ind) = t(ind)+1;
end
ts = t./1000;
if nan_interp && any(~nanind)
    nanind = (isnan(x) & isnan(y));
    xn = interp1(ts(~nanind), x(~nanind), ts(nanind), 'linear');
    yn = interp1(ts(~nanind), y(~nanind), ts(nanind), 'linear');
    x(nanind) = xn; y(nanind) = yn;
end
behav.wasnan = nanind;

x = arena_radius*(x-params.center(1))./params.radius;
y = arena_radius*(y-params.center(2))./params.radius;
% x = movmedian(x, round(params.behav_fps/2));
% y = movmedian(y, round(params.behav_fps/2));
ksize = round(params.behav_fps*smoothing_interval);
kern = ones(ksize, 1); kern = kern./sum(kern(:));
x = conv(x, kern, 'same');
y = conv(y, kern, 'same');
% distsLED = sqrt((ax - params.center(1)).^2 + (ay - params.center(2)).^2);
% distsLED = sqrt((behav.arena_x - params.center(1)).^2 + (behav.arena_y - params.center(2)).^2);
% distsRAT = sqrt((behav.x - params.center(1)).^2 + (behav.y - params.center(2)).^2);

ax = behav.arena_x;
ay = behav.arena_y;
nanind = behav.arena_x==0 & behav.arena_y == 0;
ax(nanind) = NaN;
ay(nanind) = NaN;
ax = arena_radius*(ax-params.center(1))./params.radius;
ay = arena_radius*(ay-params.center(2))./params.radius;

distsRAT = sqrt((x -  0).^2  + (y - 0).^2);
distsLED = sqrt((ax - 0).^2 + (ay - 0).^2);
nanind = distsLED<=50 & distsLED<=-50; % 100 cm ring threshold for excluding arena LED
ax(nanind) = NaN;
ay(nanind) = NaN;
%     t = behav.timestamps;

%
if nan_interp && any(~nanind)
    nanind = (isnan(ax) & isnan(ay));
    xn = interp1(ts(~nanind), ax(~nanind), ts(nanind), 'spline');
    yn = interp1(ts(~nanind), ay(~nanind), ts(nanind), 'spline');
    ax(nanind) = xn; ay(nanind) = yn;
end
% ax = conv(ax, kern, 'same');
% ay = conv(ay, kern, 'same');
% ax = movmedian(ax, ksize);
% ay = movmedian(ay, ksize);



[theta1,rho1] = cart2pol(x,y);
[theta2,rho2] = cart2pol(ax,ay);
rho2 = movmedian(rho2, ksize);
% ay = movmedian(ay, ksize);
rhodiff = rho2-median(rho2); % small fluctation of maze position from being offcenter
[xx, yy] = pol2cart(theta1-theta2, rho1-rhodiff);
[x, y] = pol2cart(theta1, rho1-rhodiff);
[axx, ayy] = pol2cart(theta2-theta2, rho2-rhodiff);
[ax, ay] = pol2cart(theta2, rho2-rhodiff);

% arena.speed
% arena.x
% arena.y
% arena.dt

behav_dt = [median(diff(behav.timestamps)); diff([behav.timestamps])]/1000;
behav.spd_roomframe = sqrt(diff([x(1); x]).^2 + diff([y(1); y]).^2)./behav_dt;
behav.spd_arenaframe = sqrt(diff([xx(1); xx]).^2 + diff([yy(1); yy]).^2)./behav_dt;
% adt = [median(diff(behav.timestamps)); diff([behav.timestamps])]/1000;
% axspd = sqrt(diff([ax(1); ax]).^2 + diff([ay(1); ay]).^2)./adt;
% axxspd = sqrt(diff([axx(1); axx]).^2 + diff([ayy(1); ayy]).^2)./adt;
if true
    figure(3); clf;
    subplot(2,2,1)
    hold on;
    scatter3(behav.arena_x, behav.arena_y, ts, 1, 'k.')
    scatter3(behav.x, behav.y, ts, 1, 'b.')
    title('Extracted')
    set(gca, 'View', [-60 60])
    axis square
    
    subplot(2,2,2)
    hold on;
    plot3(ax, ay, ts, 'k-')
    plot3(x, y, ts, 'b-')
    title('Interp & Scaled')
    set(gca, 'View', [-60 60])
    axis square
    
    subplot(2,2,3)
    hold on;
    plot3(x, y, ts, 'b-')
    plot3(xx, yy, ts, 'r-')
    title('Room (b) vs. Arena (r)')
    set(gca, 'View', [-60 60])
    axis square
    
    subplot(4,2,6)
    % hold on;
    plot(ts, behav.spd_roomframe, 'b-')
    title('Room speed (b)')
    axis tight
    subplot(4,2,8)
    plot(ts, behav.spd_arenaframe, 'r-')
    title('Arena speed (r)')
    % set(gca, 'View', [-60 60])
    axis tight
end
%% Cull contours
[good_inds, bad_inds, ~, bounds] = Draw_contour_bounding(caiman_data.fullA, ...
    caiman_data.dims, caiman_data.maxFrame, [], false);
neuron = remove_segments(caiman_data, [caiman_data.idx_components_bad, bad_inds'], true);
%%
ms = [];
ms.spatialDownsample = crop_params.spatialDownSample;
ms.temporalDownsample = crop_params.temporalDownSample;
ms.fileName = crop_params.tiffStackout;
[ms.height, ms.width] = size(imread(ms.fileName,1));
ms.frameNum = TS_data.FrameNumber(1:ms.temporalDownsample:end);
ms.timestamps = TS_data.TimeStamp_ms_(1:ms.temporalDownsample:end);
% ms.buffer = TS_data.BufferIndex;

ms.room.x = interp1(t, x, ms.timestamps, 'linear');
ms.room.y = interp1(t, y, ms.timestamps, 'linear');
ms.room.led_x = interp1(t, ax, ms.timestamps, 'linear');
ms.room.led_y = interp1(t, ay, ms.timestamps, 'linear');
ms.arena.x = interp1(t, xx, ms.timestamps, 'linear');
ms.arena.y = interp1(t, yy, ms.timestamps, 'linear');
ms.arena.led_x = interp1(t, axx, ms.timestamps, 'linear');
ms.arena.led_y = interp1(t, ayy, ms.timestamps, 'linear');
ms.neuron = neuron;
ms_dt = [median(diff(ms.timestamps)); diff([ms.timestamps])]/1000;
ms.dt = ms_dt;
ms.room.speed = sqrt(diff([ms.room.x(1); ms.room.x]).^2 + diff([ms.room.y(1); ms.room.y]).^2)./ms_dt;
ms.arena.speed = sqrt(diff([ms.arena.x(1); ms.arena.x]).^2 + diff([ms.arena.y(1); ms.arena.y]).^2)./ms_dt;

nbins = length(pos_bins)-1;
[vmap_room, vmap_room_counts]   = make_occupancymap_2D(ms.room.x, ms.room.y, ms.timestamps/1000, pos_bins, pos_bins);
[vmap_arena, vmap_arena_counts] = make_occupancymap_2D(ms.arena.x, ms.arena.y, ms.timestamps/1000, pos_bins, pos_bins);
vmap_room(vmap_room<1) = NaN;
vmap_arena(vmap_arena<1) = NaN;
nsegs = size(ms.neuron.S, 1);
pfield_room = NaN(nsegs, nbins, nbins);
pfield_arena = NaN(nsegs, nbins, nbins);
pfield_room_smooth = NaN(nsegs, nbins, nbins);
pfield_arena_smooth = NaN(nsegs, nbins, nbins);
spks = normalize_rows(ms.neuron.S);
smoothing_kern = ones(5,5); smoothing_kern = smoothing_kern./(sum(smoothing_kern(:)));

ms.room.speed_smooth = conv(ms.room.speed,   kern, 'same');
ms.arena.speed_smooth = conv(ms.arena.speed, kern, 'same');
is_moving = ms.arena.speed_smooth>speed_thresh;
for i = 1:nsegs
    [smap, ~] = make_occupancymap_2D(ms.room.x(is_moving),  ms.room.y(is_moving),  spks(i,(is_moving)), pos_bins, pos_bins);
    smap(vmap_room<1) = NaN;
    pfield_room(i,:,:) = smap./vmap_room;
    smap(isnan(smap)) = 0;
    sm_smap = conv2(smap, smoothing_kern, 'same');
    sm_smap(isnan(vmap)) = NaN;
    pfield_room_smooth(i,:,:) = sm_smap./vmap_room;
    
    [smap, ~] = make_occupancymap_2D(ms.arena.x(is_moving), ms.arena.y(is_moving), spks(i,(is_moving)), pos_bins, pos_bins);
    smap(vmap_arena<1) = NaN;
    pfield_arena(i,:,:) = smap./vmap_arena;
    smap(isnan(smap)) = 0;
    sm_smap = conv2(smap, smoothing_kern, 'same');
    sm_smap(isnan(vmap)) = NaN;
    pfield_arena_smooth(i,:,:) = sm_smap./vmap_arena;
end
%%
nc = ceil(sqrt(nsegs));
figure(10); clf; 
set(gcf, 'Name', 'ARENA FRAME')
colormap(plasma)
figure(11); clf
set(gcf, 'Name', 'ROOM FRAME')
colormap(viridis)

arena_alpha = ~isnan(vmap_arena);
room_alpha = ~isnan(vmap_room);
imAlpha(isnan(vmap_room))=0;
    figure(10);
for i = 1:nsegs
    subplot_tight(nc, nc, i, [.01 .01])
    p = squeeze(pfield_arena_smooth(i,:,:));
    p = p./max(p(:));
    imagesc(p, 'AlphaData', arena_alpha);
    axis square off
end
    figure(11);
for i = 1:nsegs
    subplot_tight(nc, nc, i, [.01 .01])
    p = squeeze(pfield_room_smooth(i,:,:));
    p = p./max(p(:));
    imagesc(p, 'AlphaData', room_alpha);
    axis square off
end

%%
spks = normalize_rows(ms.neuron.S)>0;
% spks = normalize_rows(ms.neuron.C+ms.neuron.YrA);
spks(isnan(spks)) = 0;
% [pco_tau, pco_prob, sumspks] = Fenton_pco(spks, 11*5, false, 'Kendall');
samplerate = round(1/median(ms_dt));
[pco_tau, pco_prob, sumspks] = Fenton_pco(spks, samplerate*10, false, 'Kendall');
figure(4); clf;
subplot(1,3,1)
imagesc(normalize_rows(sumspks));
subplot(1,3,2)
imagesc(pco_tau, [0 .3]);
subplot(1,3,3)
imagesc(pco_prob);

%%
% figure(9); clf; hold on
% spks = normalize_rows(ms.neuron.C);
% spks2 = normalize_rows(ms.neuron.S);
% scale = .8;
% for i = 1:size(spks,1)
%     s = spks2(i,:)*scale;
%     plot(s + i - 1, 'k')
%     s = spks(i,:)*scale;
%     plot(s + i - 1, 'r')
% end
ori_ts = ORI_Data.TimeStamp_ms_;
q = [ORI_Data.qw, ORI_Data.qx, ORI_Data.qy, ORI_Data.qx];
[~, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
ms.ori = ORI_Data;
ms.time = ms.timestamps;
[ori_str] = LFOV_ori_extract(q, [0:30:360]);


ms.roll = roll; ms.pitch = pitch; ms.yaw = yaw;
figure; hold on; plot(roll, 'r.'); plot(pitch, 'm.'); plot(yaw, 'y.')
%%
figure(3); clf;
subplot(2,2,1)
hold on; 
scatter3(behav.arena_x, behav.arena_y, ts, 1, 'k.')
scatter3(behav.x, behav.y, ts, 1, 'b.')
title('Extracted')
set(gca, 'View', [-60 60])
axis square

subplot(2,2,2)
hold on;
plot3(ax, ay, ts, 'k-')
plot3(x, y, ts, 'b-')
title('Interp & Scaled')
set(gca, 'View', [-60 60])
axis square

subplot(2,2,3)
hold on;
plot3(x, y, ts, 'b-')
plot3(xx, yy, ts, 'r-')
title('Room (b) vs. Arena (r)')
set(gca, 'View', [-60 60])
axis square

subplot(4,2,6)
% hold on;
plot(ts, behav.spd_roomframe, 'b-')
title('Room speed (b)')
axis tight
subplot(4,2,8)
plot(ts, behav.spd_arenaframe, 'r-')
title('Arena speed (r)')
% set(gca, 'View', [-60 60])
axis tight