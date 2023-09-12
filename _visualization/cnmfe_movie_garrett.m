%% cnmfe_movie_fernando
% You should change these to match your filenames and cnmfe variables!

% video_outputname = 'test.avi';
% the raw imaging I have is usually combined into one single tiff file that
% matches the frame number and size of the cnmf-e output
% tiff_name = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_12\16_55_19_RET8\HPC_miniscope1\msCam_MC.tiff'; % file can be made in matlab, imageJ, or python easily
tiff_name = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_08_16\12_46_12_CON18\HPC_miniscope1\msCam_MC.tiff'; % file can be made in matlab, imageJ, or python easily

video_outputname = 'E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\2023_08_16\12_46_12_CON18\HPC_miniscope1\test.avi';

frames = 1:8:10000; % frames to display
segs = 1:size(C,1); % which neurons to plot, use this to draw all
traces_scale = 5; % multiplier for the traces to make them easier to see (minimum of 1)
% npix = sum(fullA>0, 1);
% segs = npix<=400;

% spatial = ms.neuron.Yfull; % the spatial contours from cnmfe, usually Y
% activity = ms.neuron.S; % the temporal component to modulate contours by from cnmfe, usually C or S
% spatial = A(:,segs); % the spatial contours from cnmfe, usually A
spatial = fullA(:,segs); % the spatial contours from cnmfe, usually A
activity = C(segs, frames); % the temporal component to modulate contours by from cnmfe, usually C or S
dims = [dims(1) dims(2)]; % dimensions of the movie frame, height x width

%%

Y = double(spatial*activity);

% now we normalize and add an offset for easy plotting
traces = activity;

mint = nanmin(traces, [], 2)*ones(1,size(traces,2));
traces = traces-mint;
maxt = nanmax(traces, [], 2)*ones(1,size(traces,2));
traces = traces./maxt;
traces(isnan(traces)) = 0;
norm_activity = traces;
traces = traces.*traces_scale;
offset = linspace(0,size(traces,1)-1, size(traces,1))' + ones(1,size(traces,2));
traces = traces + offset;

% turning Y into a colored version with r g b defined by jet colormap
cmap = viridis(size(spatial,2));
[~, ord] = sort(rand(size(spatial,2),1));
cmap = cmap(ord,:);

Yr = spatial * (activity .* ( cmap(:,1)*ones(1,size(traces,2)) ));
Yg = spatial * (activity .* ( cmap(:,2)*ones(1,size(traces,2)) ));
Yb = spatial * (activity .* ( cmap(:,3)*ones(1,size(traces,2)) ));
%
y_max = nanmax(Y(:));
y_min = -nanmean(Y(:));

time_win = 90; % number of samples pre and post to plot

%%
figure;
set(gcf, 'Color', 'k', 'Position', [300 300 500 500]) % change this to change the output size
writevid = true;
if writevid==true
v = VideoWriter(video_outputname);
v.Quality = 100;% quality in percent compression
v.FrameRate = 30;
v.open()
end
%
immax = ceil(.85*double(max(max(maxFrame-minFrame))));

for i = 1+time_win:length(frames)-time_win
    %%
    fidx = frames(i);
    twin = i-time_win:i+time_win;
    im = imread(tiff_name, fidx);
    subplot(2,2,1);
    colormap gray
    imagesc(im-uint8(minFrame), [0 immax]) % change this scaling if needed to see raw activity, between 0 and 255
    axis image off
    
    f1 = (reshape(Yr(:,i), dims) + y_min ) ./ y_max;
    f2 = (reshape(Yg(:,i), dims) + y_min ) ./ y_max;
    f3 = (reshape(Yb(:,i), dims) + y_min ) ./ y_max;
    f = cat(3, f1, f2, f3);
    subplot(2,2,2)
    image(f +.2)
    title(sprintf('Frame = %5.0f', fidx), 'Color', 'w')
    axis image off
%     imagesc(f, [y_min y_max])
    subplot(2,2,[3 4]); cla;
    imagesc(1-norm_activity(:,  twin), [-.1 1.1]); hold on;
    plot([time_win+1 time_win+1], [-1*traces_scale length(segs) + 2*traces_scale], '-', 'Color', [.7 .7 .7])
%     for ii = 1:length(segs)
%     plot(traces(ii,  twin)', 'Color', cmap(ii,:)/1.5)
%     end
    axis tight
    ylim([-1*traces_scale length(segs) + 2*traces_scale])
    set(gca, 'Color', 'w', 'XColor', 'w', 'YColor', 'w', 'XTick', []); % time_win+1, 'XTickLabel', fidx)
    
    drawnow; 
    
if writevid==true
    temp = getframe(gcf);
    v.writeVideo(temp.cdata)
end
end
if writevid==true
v.close() % make sure you close it before you try to open/delete/move the file
end
























