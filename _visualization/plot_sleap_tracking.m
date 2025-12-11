function [viddata, vidtime] = plot_sleap_tracking(sleap, frames, tail, mean_flag, vidName, cam_info, room)
np = length(sleap.node_names);
nf = length(sleap.instance_scores);
save_video_fname = []; % 'C:\Users\gjb326\RAWDATA\TRACKER\Garrett_SMART_task\SLEAP\34990_IL5.avi';
%%%%%%%
sleap.tracks = room.sleap_pos;
sleap.tracks(:,:,2) = sleap.tracks(:,:,2)*-1;
%%%%%%%
x = squeeze(sleap.tracks(:, :, 1));
y = squeeze(sleap.tracks(:, :, 2));
% figure; plot(x,y)
if ~isempty(cam_info)
    
%     x = (x+rescale_xys(1))*rescale_xys(3);
%     y = (y+rescale_xys(2))*rescale_xys(3);
    
    x = (cam_info.sleap_arena_size_px.*x./40)+cam_info.sleap_pos_center_X;
    y = (cam_info.sleap_arena_size_px.*y./40)+cam_info.sleap_pos_center_Y;
end
if isempty(frames)
    frames = 1:nf;
end
% hold on; plot(x,y)
mx = mean(x,2);
my = mean(y,2);

if ~isempty(vidName)
[Y, vt, badcount] = load_tracker_avi(vidName, [], max(frames), room.fileName);
% [Y, vt, badcount] = load_tracker_avi(vidName, [], [], room_datfilename);
[h2,w2, nf_vid] = size(Y);
h1=0; w1 = 0;
else
    Y = 0;
    h1 = min(my)*1.1;
    h2 = max(my)*1.1;
    w1 = min(mx)*1.1;
    w2 = max(mx)*1.1;
%     y =y*-1;
end

bx = sleap.tracks(:, contains(sleap.node_names,'body_center'), 1);
by = sleap.tracks(:, contains(sleap.node_names,'body_center'), 2);
hx = sleap.tracks(:, contains(sleap.node_names,'head_base'), 1);
hy = sleap.tracks(:, contains(sleap.node_names,'head_base'), 2);
nx = sleap.tracks(:, contains(sleap.node_names,'nose'), 1);
ny = sleap.tracks(:, contains(sleap.node_names,'nose'), 2);

head_ang = rad2deg(atan2(nx-hx, ny-hy));
body_ang = rad2deg(atan2(hx-bx, hy-by));
mxy = squeeze(mean(sleap.tracks,2));
% mxy(:,2) = mxy(:,2)*-1;
%%
clrs = magma(np*3);
win = 0;
% clrs = clrs(floor(np/2):floor(np/2)+np, :);
clrs = clrs(np*2:-1:np*1, :);
figure(); 
% set(gcf, 'Position', [300 300 100+(w2-w1)/2 100+(h2-h1)/2])
set(gcf, 'Position', [300 300 100+(w2-w1)/2 100+(h2-h1)/2])
hold on
im=0;
if ~isempty(save_video_fname)
    vv = VideoWriter(save_video_fname);
    vv.Quality=80;
    vv.FrameRate = 30;
    vv.open();
end
node2tail={'nose' 'body_center'};
viddata = cell(length(frames),1);
vidtime = vt(frames);
for i = 1:length(frames)
    %%
    clf; %hold on
    f = frames(i);
    ftail = f-tail:f;
    ftail = ftail(ftail>0 & ftail<=nf);
    if any(Y)
        subplot_tight(1,1,1,[0,0])
        im = squeeze(Y(:,:,f));
%         im = conv2(im, ones(5,5)./25, 'same');
        imagesc(im, [0 150]);
        colormap gray
        set(gca, 'XTick', [], 'YTick', [])
    end
    hold on
    for ni = 1:length(sleap.edge_inds)
        a = sleap.edge_inds(1, ni)+1;
        b = sleap.edge_inds(2, ni)+1;
        plot([x(f, a) x(f, b)], [y(f, a) y(f, b)], 'Color', 'k', 'LineWidth', 2)
    end
    for ni = 1:np
        if mean_flag==true; plot(mx, my, 'Color', [.8 .8 .8]); end
        if any(contains(sleap.node_names{ni}, node2tail))
        plot(x(ftail, ni), y(ftail,ni), 'Color', clrs(ni,:))
        scatter(x(f, ni), y(f,ni), 30, '.', 'MarkerEdgeColor', clrs(ni,:))
        end
    end
    axis([w1 w2 h1 h2])
    
    if any(im(:)) && win>0
        subplot_tight(4,5,1,[0,0])
        c1 = max(1, round(mx(f)-win));
        c2 = min(w2, round(mx(f)+win));
        r1 = max(1, round(my(f)-win));
        r2 = min(h2, round(my(f)+win));
        im2 = im(r1:r2,c1:c2);
        im2 = imrotate(im2, -head_ang(f)+180, 'bilinear', 'crop');
        imagesc(im2, [0, 150])
        axis image off
    end
    temp = getframe(gcf);
    viddata{i} = temp.cdata;
    drawnow
    if ~isempty(save_video_fname)
        writeVideo(vv, temp)
    end

%     pause((1/30))
end
if ~isempty(save_video_fname)
    vv.close()
end
