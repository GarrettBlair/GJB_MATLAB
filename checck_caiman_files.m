topdir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\';
fall = dir(topdir);
%% checking extraction and motion corrections
for i = 1:length(fall)
    %%
    if ~fall(i).isdir
        sess = fall.name;
        fname = [fall(i).folder '/' fall(i).name];
        load(fname);
        
        figure(99); clf
        
        tiffname = ms.neuron.motion_corr_tiff;
        im_min = ms.neuron.minFrame;
        im_mean = ms.neuron.meanFrame;
        im_max = ms.neuron.maxFrame;
        im_pnr = ms.neuron.pnr_im;
        im_corr = ms.neuron.corr_im;
        sma = normalize_rows(ms.neuron.S_matw);
        craw = (ms.neuron.C+ms.neuron.YrA);
        
        
        subplot_tight(3,3,1);
        imshow(im_min, [0 255])
        subplot_tight(3,3,2);
        imagesc(im_pnr); axis image off
        % imshow(im_mean, [0 255])
        subplot_tight(3,3,3);
        imagesc(im_corr.*im_pnr); axis image off
        % imshow(im_max, [0 255])
        subplot_tight(3,1,2);
        imagesc(craw)
        subplot_tight(3,1,3);
        imagesc(sma)
        input('')
    end
end