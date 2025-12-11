clear
% d_load = 'C:\Users\gjb326\Desktop\TRACKER DOCS\Tracker Videos\';
% d_save = 'C:\Users\gjb326\Desktop\SLEAP\test2\';
d_load = 'C:\Users\gjb326\Desktop\TRACKER DOCS\Tracker Videos\HPCACC34990';
dat_dir = 'C:\Users\gjb326\Desktop\TRACKER DOCS\DAT Files\';
d_save = 'C:\Users\gjb326\Desktop\SLEAP\SMART task imaging\sleap input avis\';
vids = dir([d_load '*.avi']);

runvids = false(length(vids),1);
arena_datfiles = cell(length(vids),1);
room_datfiles = cell(length(vids),1);
fprintf('\n\tVideos in  %s:', vids(1).folder)
for ii = 1:length(vids)
    fname = [vids(ii).folder '\' vids(ii).name];
    fsave_name = [d_save '\' vids(ii).name];
%     fprintf('\n%3d.  %s...', ii, fname)
    fprintf('\n%3d.  %s...', ii, vids(ii).name)
    us = find('_' == vids(ii).name);
    arenadat = [dat_dir vids(ii).name(1:us(2)-1) '_Arena.dat'];
    roomdat = [dat_dir vids(ii).name(1:us(2)-1) '_Room.dat'];
    if isfile(arenadat) && isfile(roomdat) 
        arena_datfiles{ii} = arenadat;
        room_datfiles{ii}  = roomdat;
    end
    if isfile(fname) && ~isfile(fsave_name)
        runvids(ii) = true;
        fprintf(' <-- ')
    end
end
%%
fprintf('\n\tResaving...\n')
vids = vids(runvids);
room_datfiles = room_datfiles(runvids);
arena_datfiles = arena_datfiles(runvids);
for ii = 1:length(vids)
    fname = [vids(ii).folder '\' vids(ii).name];
    fsave_name = [d_save '\' vids(ii).name];
    if isfile(fname) && ~isfile(fsave_name)
%         try
        fprintf('\n~~~%s\n', vids(ii).name);
        [Y, vt, badcount] = load_tracker_avi(fname, fsave_name, [], room_datfiles{ii}, arena_datfiles{ii});
%         catch
%             warning(sprintf('\n~~~%s\n', vids(ii).name));
%             badvids(ii) = true;
%         end
        fprintf('%d/%d  done    %2.1f%%\n', ii, length(vids), 100*ii/length(vids));
    else
%         fprintf('~ Skip %s \n', vids(ii).name);
    end
end

% need to rerun in sleap 24513_IL3
%%
if false
sleap_fname = "C:\Users\gjb326\Desktop\SLEAP\predicted outputs\SmartyRats 521 predict.014_NAPA_24510_ILCON11_TrackerVideo.analysis.h5";
fname = "C:\Users\gjb326\Desktop\TRACKER DOCS\Tracker Videos\NAPA_24510_ILCON11_TrackerVideo.avi";
ts = csvread("C:\Users\gjb326\Desktop\SLEAP\test2\NAPA_24510_ILCON11_TrackerVideo_timeStamps.csv");

% sleap_fname = "C:\Users\gjb326\Desktop\SLEAP\predicted outputs\SmartyRats 521 predict.014_NAPA_24510_ILCON11_TrackerVideo.analysis.h5";
% fname = "C:\Users\gjb326\Desktop\TRACKER DOCS\Tracker Videos\NAPA_24510_ILCON11_TrackerVideo.avi";
% ts = csvread("C:\Users\gjb326\Desktop\SLEAP\test2\NAPA_24510_ILCON11_TrackerVideo_timeStamps.csv");

% [Y, vt, badcount] = load_tracker_avi(fname, []);
% varnames = {'track_names' 'node_names' 'edge_names' 'edge_inds' 'tracks' 'track_occupancy'...
%     'point_scores' 'instance_scores' 'tracking_scores' 'labels_path' 'video_path' 'video_ind' 'provenance'};
varnames = {'node_names' 'edge_names' 'edge_inds' 'tracks' 'track_occupancy'...
    'point_scores' 'instance_scores' 'labels_path' 'video_path' 'video_ind'};
sleap = [];
for ii = 1:length(varnames)
    structstr = sprintf('sleap.%s = h5read(sleap_fname, ''/%s'');', varnames{ii}, varnames{ii});
    eval(structstr)
end
n_nodes = length(sleap.node_names);
n_edges = length(sleap.edge_names);
mxy = squeeze(mean(sleap.tracks,2));
%%
ts_dt = 1/mean(abs(diff(ts)));
ksize = 0; % round(ts_dt*.5);
kern = gausswin(ksize, 1); kern = kern./sum(kern(:));
for ii = 1:n_nodes
    x = sleap.tracks(:, ii, 1);
    y = sleap.tracks(:, ii, 2);
    nanind = isnan(x) | isnan(y);
    x(nanind) = interp1(find(~nanind), x(~nanind), find(nanind), 'linear');
    y(nanind) = interp1(find(~nanind), y(~nanind), find(nanind), 'linear');
    sleap.tracks(:, ii, 1) = x;
    sleap.tracks(:, ii, 2) = y;
    if ksize > 0
        xn = conv(x, kern, 'same');
        yn = conv(y, kern, 'same');
        xn(1:floor(ksize/2)+1) = xn(floor(ksize/2)+2);
        xn(length(x)-floor(ksize/2):length(x)) = xn(length(x)-floor(ksize/2)-1);
        yn(1:floor(ksize/2)+1) = yn(floor(ksize/2)+2);
        yn(length(x)-floor(ksize/2):length(x)) = yn(length(x)-floor(ksize/2)-1);
        sleap.tracks(:, ii, 1) = xn;
        sleap.tracks(:, ii, 2) = yn;
    end
end
bx = sleap.tracks(:, contains(sleap.node_names,'body_center'), 1);
by = sleap.tracks(:, contains(sleap.node_names,'body_center'), 2);
hx = sleap.tracks(:, contains(sleap.node_names,'head_base'), 1);
hy = sleap.tracks(:, contains(sleap.node_names,'head_base'), 2);
nx = sleap.tracks(:, contains(sleap.node_names,'nose'), 1);
ny = sleap.tracks(:, contains(sleap.node_names,'nose'), 2);

head_ang = rad2deg(atan2(nx-hx, ny-hy));
body_ang = rad2deg(atan2(hx-bx, hy-by));
%         track_names: 0
%         node_names: 10
%         edge_names: 20
%         edge_inds: 20
%         tracks: (30876, 10, 2, 1)
%         track_occupancy: (1, 30876)
%         point_scores: (30876, 10, 1)
%         instance_scores: (30876, 1)
%         tracking_scores: (30876, 1)
%         labels_path: C:/Users/gjb326/Desktop/SLEAP/SmartyRats 521 predict.slp
%         video_path: C:/Users/gjb326/Desktop/SLEAP/test2/NAPA_24522_ILCON13_TrackerVideo.avi
%         video_ind: 74
%         provenance: {}
%

i_score = sleap.instance_scores/3;
nanind = isnan(i_score);
i_score(nanind) = 1;
if ksize>0
i_score = ceil(conv(i_score, kern, 'same'));
else
i_score = ceil(i_score);
end
i_cmap = stoplight(max(i_score));
mxy = squeeze(mean(sleap.tracks,2));

save_vid = true;
if save_vid
    vv = VideoWriter('C:\Users\gjb326\Desktop\TRACKER DOCS\figures\sleap_ex.avi');
    vv.Quality=80;
    vv.FrameRate = 30;
    vv.open();
end

figure(1); colormap bone
set(gcf, 'Position', [60 302  654  492], 'Color', 'k')
ds=1;
win=60;
[h,w,nf] = size(Y);
for idx = 40000:ds:42000 % 31:ds:nf
    %%
    clf;
    subplot_tight(1,1,1,[0,0])
    im = squeeze(Y(:,:,idx));
    hold on
    imagesc(im, [20, 250])
    for e1 = 1:n_edges
        xs = [sleap.tracks(idx, sleap.edge_inds(1, e1)+1, 1), sleap.tracks(idx, sleap.edge_inds(2, e1)+1, 1)]; % zero indexed
        ys = [sleap.tracks(idx, sleap.edge_inds(1, e1)+1, 2), sleap.tracks(idx, sleap.edge_inds(2, e1)+1, 2)]; % zero indexed
        plot(xs, ys, 'm-')
        scatter(xs, ys, 5, 'Marker', 'o', 'MarkerFaceColor', i_cmap(i_score(idx), :), 'MarkerEdgeColor', 'none')
    end
    axis image off
    plot(mxy(idx-30:idx,1), mxy(idx-30:idx,2), 'w-')
    scatter(mxy(idx,1), mxy(idx,2), 30, 'Marker', 'd', 'MarkerFaceColor', i_cmap(i_score(idx), :), 'MarkerEdgeColor', 'k')
    
    subplot_tight(4,5,1,[0,0])
    c1 = max(1, round(mxy(idx,1)-win));
    c2 = min(w, round(mxy(idx,1)+win));    
    r1 = max(1, round(mxy(idx,2)-win));
    r2 = min(h, round(mxy(idx,2)+win));    
    im2 = im(r1:r2,c1:c2);
    im2 = imrotate(im2, -body_ang(idx)+180, 'bilinear', 'crop');
    imagesc(im2, [20, 250])
    axis image off
%     text(180, 40, sprintf('head ang %3.0f', round(head_ang(idx))), 'FontSize', 20, 'Color', 'r')
    
    drawnow()
%     pause(.2)
    if save_vid==true
        temp = getframe(gcf);
        writeVideo(vv, temp)
    end
end
if save_vid==true
    vv.close()
end

%%
close all
figure(2); clf; hold on;
plot(sleap.instance_scores)
figure(3); clf; hold on;
for ii = 1:n_nodes
    plot(ii + sleap.point_scores(:,ii))
end
%%
end
