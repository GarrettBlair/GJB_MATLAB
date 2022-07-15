function make_cnmfe_movie_demo(input_path, movieOut, stacked_t, time_inds)
% recording_dir = 'C:/Users/gjb326/Desktop/RecordingData/AlejandroGrau/TestMouse1/2022_05_07/17_09_41/';
% processed_file = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\processed_files\2022_06_22___17_41_24_ms_placecells_data.mat';
timestep = 2;
win=60;
kern = gb_kernel(5, 'step3');
normalize_contours = true;

if any(strfind(input_path, '.mat'))
    load(input_path);
    recording_dir = ms.parentDir;
    generate_ms = false;
else
    recording_dir = input_path;
    ms = [];
    ms.parentDir = recording_dir;
    
    generate_ms = true;
end
tiffName = sprintf('%sMiniLFOV/%s', recording_dir, 'msCam_MC.tiff');
I = imread_big(tiffName);
bg = min(I,[],3);
k = 7;
bg = uint8(conv2(bg, ones(k,k)./(k^2), 'same'));

if ~isempty(movieOut) && ~ischar(movieOut)
    error('Movie name must be string')
end
if ~islogical(stacked_t)
    resp = input('\n\nStacked traces [y]   or trace image? ', 's');
    stacked_t = strcmp('y', resp);
end
if length(time_inds)>2 || length(time_inds)<1
    error('time inds must be of size 1 or 2 (units: seconds)')
end
% params.behav_params = behav_params;
%% Make ms structure if a directory was given
if generate_ms
    startframe = 1;
    msTSFile = sprintf('%sMiniLFOV/timeStamps.csv', recording_dir);
    msCropFile = sprintf('%sMiniLFOV/Crop_params.mat', recording_dir);
    msOriFile = sprintf('%sMiniLFOV/headOrientation.csv', recording_dir);
    % Read in the coresponding files
    caimanFilename = sprintf('%s/MiniLFOV/caiman_cnmfe_out.mat', ms.parentDir);
    sminSweepFilename = sprintf('%s/MiniLFOV/deconv_sweep.mat', ms.parentDir);
    caiman_data = load(caimanFilename);
    TS_data = readtable(msTSFile);
    crop_params = load(msCropFile);
    ORI_Data = readtable(msOriFile);
    params.crop_params  = crop_params;
    params.smin_vals = [-50:5:-5];
    
    [nsegs,nframes] = size(caiman_data.C);
    
    [smat, smat_weighted, good_idx, ~] = deconv_sweep_read(sminSweepFilename, params.smin_vals);
    all_good_idx = find(sum(good_idx,1)==size(good_idx,1));
    bad_idx = setdiff(1:size(caiman_data.C,1), all_good_idx);
    caiman_data.idx_components = all_good_idx;
    caiman_data.idx_components_bad = bad_idx;
    temp = sum(smat, 1);
    caiman_data.S_mat = reshape(temp, [nsegs, nframes]);
    temp = sum(smat_weighted, 1);
    caiman_data.S_matw = reshape(temp, [nsegs, nframes]);
    
    neuron = remove_segments(caiman_data, caiman_data.idx_components_bad, false);
    [nsegs,nframes] = size(neuron.C);
    neuron.spks = normalize_rows(neuron.S_matw);
    
    ms.spatialDownsample = crop_params.spatialDownSample;
    ms.temporalDownsample = crop_params.temporalDownSample;
    ms.fileName = crop_params.tiffStackout;
    [ms.height, ms.width] = size(imread(ms.fileName,1));
    % ms.frameNum = TS_data.FrameNumber(1:ms.temporalDownsample:end);
    ms.frameNum = TS_data.FrameNumber(startframe:ms.temporalDownsample:end);
    ms.timestamps = TS_data.TimeStamp_ms_(startframe:ms.temporalDownsample:end);
    tiff_numFrames = size(imfinfo(ms.fileName),1);
    ms.neuron = neuron;
end
%%
% behavTiff = sprintf('%sBehavCam/behavCam.tiff', recording_dir);
% behavFile = sprintf('%sBehavCam/behavCam.tiff', recording_dir);
% B = imread_big(behavTiff);




%%
% stacked_t = false;
if normalize_contours
    a_thresh = 0;
    a = normalize_cols(ms.neuron.fullA);
    a(a<a_thresh) = 0;
else
    a = ms.neuron.fullA;
end
spks = ms.neuron.S_matw;
[nsegs, nframes] = size(spks);
cspks = spks*0;


for i = 1:nsegs
    cspks(i,:) = conv(spks(i,:), kern, 'same');
end
for i = 1:nsegs
    cspks(i,:) = conv(spks(i,:), kern, 'same');
end
cspks = normalize_rows(cspks);

rng(348)
step = nsegs/4;
ylabels = [ 1 ,  10*floor([1+step:step:nsegs-step+1]/10), nsegs ];
ts = ms.timestamps./1000;

% movieOut = sprintf('%sMiniLFOV/%s', recording_dir, 'cnmfe_movie.avi');
% movieOut = [];
if ~isempty(movieOut)
    if exist(movieOut, 'file')==2
        resp = input('\n\nDELETE existing movie?  ([N]/y)', 's');
        if any(strfind('yY', resp))
            delete(movieOut)
        else
            kk = datetime;
            movieOG = movieOut;
            avistart = strfind(movieOG, '.avi')-1;
            movieOut = sprintf('%s_%i_%i_%i.avi', movieOG(1:avistart), kk.Hour, kk.Minute, round(kk.Second));
            fprintf('\n\tfile: %s', movieOut)
        end
    end
    if exist('v', 'var')==1
        try
            v.close();
        catch
            v = [];
        end
    end
    v = VideoWriter(movieOut);
    realtime = round(1/median(abs(diff(ts))));
    tstime = timestep*realtime;
    v.FrameRate = 2*tstime;
    v.open()
end
if stacked_t==true
    ccolors = rand(nsegs, 3)/2 + .4;
    set(gcf, 'Position', [200   200   753   482], 'Color', 'k')%[.75 .75 .75]
    sbs = [.1 .1];
    traceScale = 10;
else
    ccolors = viridis(nsegs*2);
    ccolors = ccolors(nsegs+1-floor(nsegs/2):nsegs*2-floor(nsegs/2), :);
    [~,ord] = sort(rand(nsegs,1));
    ccolors = ccolors(ord,:);
    set(gcf, 'Position', [200   200   714   351], 'Color', 'k')%[.75 .75 .75]
    sbs = [.1 .05];
    traceScale = 3;
end
if length(time_inds)==1
    time_inds = [0 time_inds];
end
i1 = find(min(abs(ts-time_inds(1))) == abs(ts-time_inds(1)));
i2 = find(min(abs(ts-time_inds(2))) == abs(ts-time_inds(2)));
if i1 - win <1
    i1 = win+1;
end
if i2+win*4>nframes
    i2 = nframes - win*4;
end

for i = i1:timestep:i2
    %%
    spkSize = i-win:i+win*4;
%     im = 4*(squeeze(I(:,:,i)) - bg);
    im = squeeze(I(:,:,i));
    yr = a*(cspks(:,i).*ccolors(:,1));
    yg = a*(cspks(:,i).*ccolors(:,2));
    yb = a*(cspks(:,i).*ccolors(:,3));
    yr = reshape(yr, ms.neuron.dims);
    yg = reshape(yg, ms.neuron.dims);
    yb = reshape(yb, ms.neuron.dims);
    ycolor = cat(3, yr, yg, yb);
%     yim = uint8(yim*255);
%     spk_sub1 = spks(:, i-win:i+win*4);
    spk_sub2 = cspks(:, spkSize);
    sr = (spk_sub2.*ccolors(:,1));
    sg = (spk_sub2.*ccolors(:,2));
    sb = (spk_sub2.*ccolors(:,3));
    scolor = cat(3, sr, sg, sb);

    figure(1); clf
    subplot_tight(2,3,1, [.01 .01])
    imshow(im); axis image tight
    subplot_tight(2,3,4, [.01 .01])
%     image(ycolor*15)
    image(ycolor*traceScale + .2)
%     imagesc(yim, [0 1])
    axis image off
    subplot_tight(1,3,2:3, sbs)
    cla
    set(gca, 'YColor', 'w', 'XColor', 'w')
    hold on
    t = ts(spkSize);
    x = mod(t,10);
    xx = find(x(1:end-1) > x(2:end));
    if stacked_t
        for j = nsegs:-1:1
            plot(t, j-1+(spk_sub2(j,:)*traceScale), 'LineWidth', 1, 'Color', ccolors(j,:))
        end
        plot(t([win+1, win+1]), [.5 nsegs+traceScale], 'r--', 'LineWidth', 1)
        set(gca, 'XTick', round(t(xx)))
        axis tight
    else
        image(scolor.*traceScale)
        plot([win+1, win+1], [.5 nsegs+.5], 'r--', 'LineWidth', 1)
        
        set(gca, 'XTick', xx, 'XTickLabel', round(t(xx)))
        axis image
    end
    xlabel('Time (sec)')
    ylabel('Components')
    title('Deconvolved calcium activity', 'Color', 'w')
    set(gca, 'YTick', ylabels)
    hold off
    axis tight
    drawnow

    if ~isempty(movieOut)
        temp = getframe(gcf);
        v.writeVideo(temp.cdata)
    end
end
if ~isempty(movieOut)
    v.close()
end









