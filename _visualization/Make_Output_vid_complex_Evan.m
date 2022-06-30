function Make_Output_vid_complex_Evan(signal, spatial, spatial_dims, ms, behav, ds, save_flag, filename)
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
r = 4;
c = 3;
figure_dims = [r*200, c*300]; % figure dimensions to be plotted and saved
set(mov, 'Name', filename, 'Position', [100 100 figure_dims(2)+100 figure_dims(1)+100])
set(mov,'color','k')

d1s = spatial_dims(1);
d2s = spatial_dims(2);
nsegs = size(signal,1);

% number of samples to plot before and after current sample
back_window = 30;
back_window_max = 1*30;
forward_window_min = 60;
forward_window = 5*30;

msframes = 1:ds:(ms.numFrames-ds+1);
hSmall = fspecial('average', 2);

% Making pretty colors for spatial shapes
ca = .3;
rr = linspace(ca, 1-ca, nsegs);
gg = ones(nsegs,1)*.1;
bb = linspace(1-ca, ca, nsegs);
colors = [rr' gg bb'];
light_colors = [rr'+.2 gg bb'+.2];


%% threshold cell bodies
thresh_spatial = zeros(size(spatial,1), nsegs);
norm = ones(1, nsegs);
for seg = 1:nsegs
    cutoff = max(spatial(spatial(:,seg)>0, seg))*.5;
    cutoff_max = max(spatial(spatial(:,seg)>0, seg))*.6;

    bad = spatial(:,seg) < cutoff;
    thresh_spatial(:, seg) = spatial(:,seg);
    thresh_spatial(bad, seg) = 0;

    bad = spatial(:,seg) < cutoff | spatial(:,seg) > cutoff_max;
    outline(:, seg) = spatial(:,seg);
    outline(bad, seg) = 0;
    norm(seg) = max(signal(seg,:));
end
norm(:) = max(norm)/2;
%% Interpolate behavior data for miniscope timestamps
%%
% figure;
% imagesc(ms.maxFrame{1}, [0 255])
% axis image
% fprintf('\nClick the edges of the average cell...\n')
% [w,v] = ginput(2);
% kern = round(sqrt((w(1) - w(2)).^2+(v(1) - v(2)).^2));
kern = 12;
se = strel('disk', ceil(kern/2));

%%
max_f = imresize(ms.maxFrame{1}-ms.minFrame{1}, 1/ms.spatialDownSample)/255;
mean_f =  filter2(hSmall, ms.meanFrame{1}/255);
min_f =  filter2(hSmall, ms.minFrame{1}/255);
max_range = (max(ms.maxFrame{1}(:))  -  max(ms.minFrame{1}(:)))/255;

neframe = anisodiff2D(max_f-imresize(min_f, 1/ms.spatialDownSample), 5, .25, .5, 1);
bg = imopen(neframe, se);
ne_frame = neframe-bg;
max_ne_range = max(ne_frame(:));
if save_flag
    avi_file = VideoWriter(filename);
    avi_file.Quality = 100;
    avi_file.open();
end
outline_base = squeeze(sum(reshape(full(outline), d1s,d2s, nsegs), 3));
figure(mov);
for frame = 5000+back_window:7000;%size(signal,2)-forward_window_min-1
    %% Behavior plotting
    subplot_tight(r,c,1);
    set(gca, 'Color', [.2 .2 .2], 'XColor', [.9 .9 .9], 'YColor', [.9 .9 .9])
    t = ms.time(msframes(frame));
    near_behav_frame = find(min(abs(behav.time - t)) == abs(behav.time - t));
    behav_frame = msReadFrame(behav, near_behav_frame(1), 0, 0, 0);
    image(behav_frame./255)
% %         title(sprintf('Behavior,   %2.3f sec', frame*33*ds/1000), 'Color', [1 .5 .5], 'FontSize', 8)
    text(round(size(behav_frame,2)/3), size(behav_frame,1)-60, sprintf('15AQ,   %2.1f sec', behav.time(near_behav_frame)/1000), 'Color', [.8 .3 .3], 'FontSize', 8)
    axis image off
    

    %% imagesc of signal activity
    if frame + forward_window > size(signal,2)
        forward_window = forward_window - 1;
    end
    traces_ax = subplot_tight(r,c,[7:12]); 
    cla
    acitivity_norm = (norm'*ones(1, length([frame-back_window:frame+forward_window])));
    activity = signal(:, [frame-back_window:frame+forward_window]);

    imagesc(activity./acitivity_norm, [0, 1]);
    caxis([0 1])
    colormap(traces_ax, jet)
    hold on
    plot([back_window, back_window], [1 nsegs], 'r--')
    hold off
    set(gca, 'YColor', [.9 .9 .9])
    set(gca, 'YTick', [1:nsegs], 'YTickLabel', [1:nsegs], 'FontSize', 5)
    text(0, -3, 'Calcium traces', 'Color', [1 .5 .5], 'FontSize', 8)
    axis tight
    xlabel('Frame number')
    %% spatial contours of signal
    cnmf = zeros(d1s,d2s,3);
    out_signal = zeros(d1s,d2s,3);
    
    for seg = nsegs:-1:1
        n_temp = thresh_spatial(:,seg)*(signal(seg, frame)/norm(seg));
        n_temp = reshape(n_temp,d1s,d2s);
        cnmf(:,:,1) = cnmf(:,:,1) + (n_temp)*colors(seg,1);
        cnmf(:,:,2) = cnmf(:,:,2) + (n_temp)*colors(seg,2);
        cnmf(:,:,3) = cnmf(:,:,3) + (n_temp)*colors(seg,3);
        o_temp = outline(:,seg)*(signal(seg, frame)/norm(seg));
        o_temp = reshape(o_temp,d1s,d2s);
        out_signal(:,:,1) = out_signal(:,:,1) + (o_temp)*colors(seg,1);
        out_signal(:,:,2) = out_signal(:,:,2) + (o_temp)*colors(seg,2);
        out_signal(:,:,3) = out_signal(:,:,3) + (o_temp)*colors(seg,3);
    end
    %% Raw miniscope image, stabilized
    subplot_tight(r,c,2)
    raw_frame = msReadFrame(ms, msframes(frame), false, true, false)/255;
    raw_frame = filter2(hSmall, raw_frame);
    imshow(raw_frame)
    axis image
    axis off
    title(sprintf('Raw, %2.1f sec', ms.time(msframes(frame))/1000), 'Color', [1 .5 .5], 'FontSize', 8)
    daspect([1 1 1])

    subplot_tight(r,c,3)
    imshow(raw_frame-min_f)
    caxis([0 max_range*.75])
    axis image
    axis off
    title('Df/f', 'Color', [1 .5 .5], 'FontSize', 8)
    daspect([1 1 1])
    
    ne = subplot_tight(r,c,4);
    [ne_frame] =  quick_NE(raw_frame, se);
%     ne_frame = anisodiff2D(, 5, .25, .5, 1);
%     bg = imopen(ne_frame, se);
%     ne_frame = ne_frame-bg;

    imagesc(ne_frame, [0 max_ne_range*.5]);
%     caxis([0 max_range*.75])
    axis image
    axis off
    colormap(ne, 'gray')
    title('Neural enhance', 'Color', [1 .5 .5], 'FontSize', 8)

    
    %% spatial contours of signal
    subplot_tight(r,c,5); 
    cla
    image(cnmf)
    axis image
    axis off
    title('Spatial components', 'Color', [1 .5 .5], 'FontSize', 8)

    subplot_tight(r,c,6); 
    cla
    out_signal(:,:,2) = outline_base*.05;
    out_hybrid = (out_signal)*.5;
    image(out_hybrid + max_f*2)
%     hold on
%     image(out*.5 + max_im/255)
    axis image
    axis off
    title('Overlay', 'Color', [1 .5 .5], 'FontSize', 8)

    
    drawnow
    
    %%
    if save_flag
        temp = getframe(gcf);
%         temp = imresize(temp.cdata, figure_dims);
        avi_file.writeVideo(temp);
    end
    if back_window < back_window_max
        back_window = back_window + 1;
    end


end
if save_flag
    avi_file.close();
end


