function make_diss_decoding_vid_orientation(ms, behav, spk_im, contours, contour_sig, colors, inds, save_flag, filename);

% beh_file = behav.tiff_filename

% Makes a video of the data from signal and spatial
%       Inputs:
%   signal       = matrix of calcium acitivty, [segments x samples]
%   spatial      = matrix of segment ROI shapes, linearized pixels [pixels x segments]
%   spatial_dims = dimensions of analyzed video
%   ms           = miniscope structure for plotting the motion corrected raw video
%   behav        = behavior structure with x and y position coords for plotting
%   ds           = downsampling rate of the analyzed data
%   save_flag    = logical flag for whether or not to save video output
%   filename     = file name of output video and figure title
%
%       -- Garrett Blair 8/14/2018

mov = figure(1); 
clf
r = 3;
c = 3;
figure_dims = [r*200, c*300]; % figure dimensions to be plotted and saved
% set(mov, 'Name', filename, 'Position', [100 100 figure_dims(2)+100 figure_dims(1)+100])
set(mov, 'Name', filename, 'Position', [74 440 989 867])
set(mov,'color','k')

signal = spk_im;%neuron.C;
spatial = contours;%neuron.A;
ds = ms.temporalDownsample;
d1s = ms.neuron.dims(1);
d2s = ms.neuron.dims(2);
nsegs = size(signal,1);

% number of samples to plot before and after current sample
back_window = 22*10;%22*4;%
back_window_max = back_window;
forward_window_min = 60;%22*4;%
forward_window = 22*45;%22*4;%

% msframes = 1:ds:(ms.numFrames-ds+1);
msframes = 1:ms.numFrames;%ms.frameNum;
hSmall = fspecial('average', 2);

% Making pretty colors for spatial shapes
ca = .3;
rr = linspace(ca, 1-ca, nsegs);
gg = ones(nsegs,1)*.1;
bb = linspace(1-ca, ca, nsegs);
% colors = [rr' gg bb'];
% colors = colormap;%plasma(nsegs);
light_colors = colors;%[rr'+.2 gg bb'+.2];

cc = reshape(contours, d1s, d2s, nsegs);
%% threshold cell bodies
%% Interpolate behavior data for miniscope timestamps
% bx = behav.x;
% by = behav.y;
% bt = behav.time;
t = ms.time;
x = ms.x;
y = ms.y;
%%
% figure;
% imagesc(ms.maxFrame{1}, [0 255])
% axis image
% fprintf('\nClick the edges of the average cell...\n')
% [w,v] = ginput(2);
% kern = round(sqrt((w(1) - w(2)).^2+(v(1) - v(2)).^2));
kern = 12;
se = strel('disk', 3);

    a1 = signal.*(colors(:,1)*ones(1,size(signal,2)));
    a1 = conv2(1, [0 0 0 1 1 1 1]./4, a1, 'same');
    a2 = signal.*(colors(:,2)*ones(1,size(signal,2)));
    a2 = conv2(1, [0 0 0 1 1 1 1]./4, a2, 'same');
    a3 = signal.*(colors(:,3)*ones(1,size(signal,2)));
    a3 = conv2(1, [0 0 0 1 1 1 1]./4, a3, 'same');
    signal_im = cat(3, a1, a2, a3);
%     spikes_video_rgb{1} = a1;
%     spikes_video_rgb{2} = a2;
%     spikes_video_rgb{3} = a3;

    a1 = spatial.*(ones(size(spatial,1),1)*colors(:,1)');
    a2 = spatial.*(ones(size(spatial,1),1)*colors(:,2)');
    a3 = spatial.*(ones(size(spatial,1),1)*colors(:,3)');
    
%     contours_video_rgb{1} = a1;
%     contours_video_rgb{2} = a2;
%     contours_video_rgb{3} = a3;
if save_flag

    ftif = Fast_Tiff_Write(filename);
end
figure(mov);
b1 = min(x);
b2 = max(y)*1.05;
startf = 1+back_window;
zoom = false(d1s, d2s);
zoom(242:506, 178:538) = true;
% zoom(58:187, 77:251) = true;
numframes = startf+min([10000 size(signal,2)-forward_window_min-1]);
bgf =  filter2(ones(5,5)./25, ms.neuron.meanFrame)/1.3;
shock_f = 2927;

[o] = LFOV_ori_extract(ms);
yy = o.yaw;
rr = o.roll;
pp = o.pitch;
yy = movmedian(yy, 5);
rr = movmedian(rr, 5);
pp = movmedian(pp, 5);

tiffname = 'H:\Miniscope Data\Hipp8\2020_03_23\15_28_19_linear10\Miniscope\msCam_MC_v21.tiff';%ms.tiff_filename
for frame = startf:2:startf+3000%numframes
    %% Behavior plotting
%     subplot_tight(6,4, 1:2,[0, .1]);
    subplot_tight(3,2,1,[0, .05]);
    set(gca, 'Color', [.2 .2 .2], 'XColor', [.9 .9 .9], 'YColor', [.9 .9 .9]); cla
    hold on
    plot(x, y, '-', 'Color', [.2 .2 .2], 'LineWidth', 5)
    plot(x(frame-floor(back_window/10):frame), y(frame-floor(back_window/10):frame), 'r-');%, 'Color', [.7 .7 .7])
    plot(x(frame), y(frame), 'r.', 'MarkerSize', 20)
    set(gca, 'Color', 'k')
    hold off
%     axis([min(x)-.05*min(x) max(x)+.05*max(x) min(y)-.05*min(y) max(y)+.05*max(y)])
    axis([-200 200 -40 200])
%     axis square
    axis off

%     subplot_tight(6,4,5:6, [0.025, .15]); cla
% %     subplot(r,c,1);
%     set(gca, 'Color', [.2 .2 .2], 'XColor', [.9 .9 .9], 'YColor', [.9 .9 .9])
%     plot(ms.real_bin, '-', 'Color', [.7 .7 .7], 'LineWidth', 3)
%     hold on
%     plot(frame, ms.real_bin(frame), 'r.', 'MarkerSize', 30)
%     if ~isnan(ms.real_bin(frame))
%     plot(frame, ms.pred_bin(frame), 'yo', 'MarkerSize', 15)
%     end
%     patch([frame-25 frame-25 frame+25 frame+25], [1 85 85 1],'white','FaceAlpha',.2, 'EdgeColor', 'none');
%     set(gca, 'Color', 'k')
%     hold off
% %     axis([min(x)-.05*min(x) max(x)+.05*max(x) min(y)-.05*min(y) max(y)+.05*max(y)])
%     xlim([frame-back_window, frame+forward_window])
%     ylim([ 0 85])
%     axis off

    % imagesc of signal activity
    if frame + forward_window > size(signal,2)
        forward_window = forward_window - 1;
    end
%     traces_ax = subplot_tight(r,c,[4:9]); 
    % Raw miniscope image, stabilized
    subplot_tight(3,5, [6 7], [.05, 0])
%     set(gca, 'Color', 'none'); axis off
    raw_frame = imread(tiffname, frame);
    raw_frame = filter2(ones(5,5)./25, double(raw_frame));
    raw_frame = raw_frame - bgf;
    raw_frame = uint8(raw_frame);
    im = raw_frame(150:580, 100:590);
    imshow((im*3));%(reshape(raw_frame(zoom), 187-58+1, 251-77+1));
    axis image
    axis off
    title('DF/f', 'Color', [1 .5 .5], 'FontSize', 16)
    daspect([1 1 1])


    
    % spatial contours of signal
%     subplot_tight(r,4,3:4, [.05, 0]); 
    subplot_tight(3,5, [11 12], [.05, 0])
%     subplot(r,c,3); 
    cla
%     act = act(end:-1:1,:);
%     imagesc(reshape(act(zoom), 187-58+1, 251-77+1), [0 20]); colormap viridis
    act1  = reshape(a1*contour_sig(:,frame), d1s, d2s);
    act2  = reshape(a2*contour_sig(:,frame), d1s, d2s);
    act3  = reshape(a3*contour_sig(:,frame), d1s, d2s);
    act_im = double(cat(3, act1, act2, act3));

%     imagesc(act, [0 20]); colormap colors
    image(act_im(150:580, 100:590, :).*5);
    axis image
    axis off
    title('Spatial components', 'Color', [1 .5 .5], 'FontSize', 16)

    
%     subplot_tight(3,4,[2:3], [0.1, .05]); cla
    subplot_tight(3,5,[8 9 10 13 14 15], [0.1, .05]); cla
    cla
    activity = signal_im(:, [frame-back_window:frame+forward_window],:);
    image(activity);
%     axis image
    hold on
    patch([back_window-25 back_window-25 back_window+25 back_window+25], [1 nsegs nsegs 1],'white','FaceAlpha',.2, 'EdgeColor', 'none');
    hold off
    set(gca, 'YColor', [.9 .9 .9])
    set(gca, 'XColor', [.9 .9 .9])
    set(gca, 'YTick', round(linspace(1, nsegs, 10)), 'XTick', [1+back_window], 'XTickLabel', sprintf('%i sec', round(ms.time([frame])./1000)), 'FontSize', 12)
%     text(30, -3, 'Normalized calcium traces', 'Color', [1 .5 .5], 'FontSize', 16)
    axis tight
    xlabel('Session time')


        subplot_tight(3,2,2, [0.1, .05]); cla
    hold on
    backsub = round(back_window/4);
    set(gca, 'Color', 'none', 'YColor', 'w', 'XColor', 'w', 'YTick', [-1*pi pi], 'YTickLabel', {'-pi' 'pi'}, 'XTick', [1+backsub], 'XTickLabel', sprintf('%i sec', round(ms.time([frame])./1000)), 'FontSize', 12)
    patch([backsub-5 backsub-5 backsub+5 backsub+5], [1.4*pi -1.4*pi -1.4*pi 1.4*pi],'white','FaceAlpha',.3, 'EdgeColor', 'none');
    plot(rr(frame-backsub:frame+backsub), 'c', 'LineWidth', 2)
    plot(pp(frame-backsub:frame+backsub), 'y', 'LineWidth', 2)
    plot(yy(frame-backsub:frame+backsub), 'm', 'LineWidth', 2)
    axis tight
    ylim([-1.4*pi 1.4*pi])
    title('Head orientation', 'Color', [1 .5 .5], 'FontSize', 16)
    text(1+backsub+10, 1.3*pi, 'roll', 'Color', 'c', 'FontSize', 16)
    text(1+backsub+25, 1.3*pi, 'pitch', 'Color', 'y', 'FontSize', 16)
    text(1+backsub+45, 1.3*pi, 'yaw', 'Color', 'm', 'FontSize', 16)
    
    drawnow
    %%
    if save_flag
        temp = getframe(gcf);
%         temp = imresize(temp.cdata, figure_dims);
%         avi_file.writeVideo(temp);
        im = temp.cdata;
        im = imrotate(im, -90);
        im = im(:, end:-1:1, :);
        ftif.WriteIMG(im);
    end
    if back_window < back_window_max
        back_window = back_window + 1;
    end


end
if save_flag
%     avi_file.close();
    ftif.close;
end