function Make_ACC_behavior_vid(pr, prc)
% inputs -
% pr - a data structure where pr (no choice) and prc (choice sess) contains the following
% relevant data: (same for both pr and prc structures)
%       - pr.leverpresses   = logical vector of lever presses (shape - 1 x frames)
%       - pr.headentries    = logical vector of head entry (shape - 1 x frames)
%       - pr.time           = logical vector of head entry (shape - 1 x frames)
%       - pr.S              = deconvolved neural activity (shape - cells x frames)
%
fs = 7.5; % sampling rate of miniscope
vid_length = 80; % seconds to show in both recordings
numf            = round(fs*vid_length); % number of frames to video 
pr_frames       = [490:490+numf-1]; % frames to use, pr
prc_frames      = [4500:4500+numf-1]; % frames to use, prc

video_output = 'D:\testvideo';
pr_recording_dir    = 'D:\20AQ\PR recording\4_13_2019\H9_M31_S38\'; % directory with behavior vids
behav_pr            = msGenerateVideoObj(pr_recording_dir, 'behavCam');
pr_leverpress       = pr.leverpresses(pr_frames); % PR lever press
pr_headentry        = diff([0  pr.leverpresses ]) > 0; % PR head entry
pr_headentry        = pr_headentry(pr_frames);


prc_recording_dir   = 'D:\20AQ\PRC\4_18_2019\H9_M28_S1\'; % directory with behavior vids
behav_prc           = msGenerateVideoObj(prc_recording_dir, 'behavCam');
prc_leverpress      = prc.leverpresses(prc_frames); % PRC lever press
prc_headentry       = diff([0  prc.leverpresses]) > 0; % PR head entry
prc_headentry       = prc_headentry(prc_frames);
     
pr_activity     = 1*(full(pr.S(:, pr_frames))>0);
prc_activity    = 1*(full(prc.S(:, prc_frames))>0);

wins        = 200; % size of the ca imaging
prcrop      = [30 100; 30+wins-1 100+wins-1]; % for cropping the ca imaging, x and y inds of video
prccrop     = [30 100; 30+wins-1 100+wins-1]; % for cropping the ca imaging

pr_miniscope_time   = pr.time(1:4:end); % time vector for pr
pr_miniscope_time   = pr_miniscope_time(pr_frames);
prc_miniscope_time  = prc.time(1:4:end); % time vector for prc
prc_miniscope_time  = prc_miniscope_time(prc_frames);
%
%
pr_behav_vid    = NaN(behav_pr.height, behav_pr.width, numf); % beahvior video matrix
prc_behav_vid   = NaN(behav_prc.height, behav_prc.width, numf); % beahvior video matrix
prc_imaging_vid = NaN(wins, wins, numf); % ca video matrix
pr_imaging_vid  = NaN(wins, wins, numf);% ca video matrix

for i = 1:numf % load in the videos to memory, keeping the timing the same
    %%
    t1 = pr_miniscope_time(i);
    pr_behav_f = find(min(abs(behav_pr.time - t1)) == abs(behav_pr.time - t1), 1);
    
    t2 = prc_miniscope_time(i);
    prc_behav_f = find(min(abs(behav_prc.time - t2)) == abs(behav_prc.time - t2), 1);
    
%     pr_imaging_vid(:,:,i) = imread(pr.tiff_filename, pr_frames(i)); % alternatively load them with imread
    pr_behav_vid(:,:,i) = rgb2gray(msReadFrame(behav_pr, pr_behav_f, 0, 0, 0)./255);
    temp = imread(pr.ms.tiff_filename, pr_frames(i));
    pr_imaging_vid(:,:,i) = temp(prcrop(1,1):prcrop(2,1), prcrop(1,2):prcrop(2,2));
    
%     prc_imaging_vid(:,:,i) = imread(prc.tiff_filename, prc_frames(i)); % alternatively load them with imread
    prc_behav_vid(:,:,i) = rgb2gray(msReadFrame(behav_prc, prc_behav_f, 0, 0, 0)./255);
    temp = imread(prc.ms.tiff_filename, prc_frames(i));
    prc_imaging_vid(:,:,i) = temp(prccrop(1,1):prccrop(2,1), prccrop(1,2):prccrop(2,2));
end
%% Writing the video
figure(46); clf
h = VideoWriter(video_output, 'MPEG-4');
set(gcf, 'Position', [158 401 2040 811], 'Color', [.25 .25 .25])
h.Quality = 100;
h.FrameRate = 15;
h.open;

for i = 1:numf
    %%
    subplot(5, 4, [1 5])
    imshow(pr_behav_vid(:,:,i));
    title('\color{white}Behavior')
    set(gca, 'FontSize', 16)
    axis image off
    text(390, -90, 'No choice session (2x Speed)', 'Color', 'w', 'FontSize', 16)
    
    subplot(5, 4, [2 6])
    imagesc(pr_imaging_vid(:,:,i), [0 100]);
    title('\color{white}Imaging')
    set(gca, 'FontSize', 16)
    axis off image
    
    subplot(5, 4, [3 7])
    imshow(prc_behav_vid(:,:,i));
    title('\color{white}Behavior')
    set(gca, 'FontSize', 16)
    text(420, -90, 'Choice session (2x Speed)', 'Color', 'w', 'FontSize', 16)
    
    subplot(5, 4, [4 8])
    imagesc(prc_imaging_vid(:,:,i), [0 100]);
    title('\color{white}Imaging')
    set(gca, 'FontSize', 16)
    axis off image
    colormap viridis
    
    subplot(5,4, [9 10]); cla
    hold on
    plot(pr_headentry, 'w')
    plot(pr_leverpress, 'g')
    plot([i i], [0 1.25], '--m')
    ylim([0 1.25])
    text(150, 1.25, '\color{green}Lever Presses  \color{white}Head Entries', 'FontSize', 16)
    axis off
    
    subplot(5,4, [11 12]); cla
    hold on
    plot(prc_headentry, 'w')
    plot(prc_leverpress, 'g')
    plot([i i], [0 1.25], '--m')
    text(150, 1.25, '\color{green}Lever Presses  \color{white}Head Entries', 'FontSize', 16)
    ylim([0 1.25])
    axis off
    
    subplot(5,4, [13 14 17 18]); cla
    hold on
    imagesc(pr_activity)
    plot([i i], [0 size(pr_activity,1)], '--m')
    text(size(pr_activity,2)/2 - 100, size(pr_activity,1) + 7, '\color{white}Calcium Transients', 'FontSize', 16)
    axis tight off
    
    subplot(5,4, [15 16 19 20]); cla
    hold on
    imagesc(prc_activity)
    plot([i i], [0 size(prc_activity,1)], '--m')
    text(size(prc_activity,2)/2 - 100, size(prc_activity,1) + 7, '\color{white}Calcium Transients', 'FontSize', 16)
    axis tight off
    
    drawnow
    mov = getframe(gcf);
    h.writeVideo(mov.cdata);
end
h.close;




