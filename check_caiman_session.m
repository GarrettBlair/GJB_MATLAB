topdir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\';
matchingfile = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\matching_contours\matching_matrix.mat';
spk_fileloc = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\matched_spk_mats\';
load(matchingfile, 'cellmap', 'file_names');
f_all = dir(topdir);
%% checking extraction and motion corrections
figure(9); clf;
p_all = [];
for i = 1:length(f_all)
    %%
    if ~f_all(i).isdir
        sess = f_all(i).name;
        fname = [f_all(i).folder '/' f_all(i).name];
        load(fname);
        this_sess = contains(file_names, sess);
        cmap_id = cellmap(:, this_sess);
        
        tiffname = ms.neuron.motion_corr_tiff;
        im_min = ms.neuron.minFrame;
        im_mean = ms.neuron.meanFrame;
        im_max = ms.neuron.maxFrame;
        im_pnr = ms.neuron.pnr_im;
        im_corr = ms.neuron.corr_im;
        sma = normalize_rows(ms.neuron.S_matw);
        spk_bin = bin_spks_time(sma, 1, ms.timestamps./1000, false);
        craw = normalize_rows(ms.neuron.C+ms.neuron.YrA);
        
        full_spk_mat = zeros(size(cellmap,1), size(spk_bin,2));%size(ms.neuron.S,2));
        full_spk_mat(cmap_id(cmap_id>0),:) = spk_bin;
        save([spk_fileloc sess], 'full_spk_mat')
        p = corr(craw');
        p(find(eye(size(p))))=NaN;
        i;
        thresh = .35;%nanmean(p(:))+4*nanstd(p(:))
        potential_bad = find(max(p, [], 1)>thresh); %  .35)%
        figure(9); hold on; histogram(p(:))
        p_all = cat(1, p_all, p(:));
        
        figure(99); clf
        subplot_tight(3,3,1);
        imshow(im_min, [0 255])
        subplot_tight(4,3,2);
        imagesc(p); axis image off
%         imagesc(im_pnr); axis image off
        % imshow(im_mean, [0 255])
        subplot_tight(4,3,3);
        imagesc(im_corr.*im_pnr); axis image off
        % imshow(im_max, [0 255])
        subplot_tight(4,1,2);
        imagesc(craw)
        subplot_tight(4,1,3:4);
        stacked_traces(craw, 1, {'k', 'LineWidth', 1.5})
        stacked_traces(sma, 1, {'r'})
%         imagesc(sma)
%         input('')
    end
end