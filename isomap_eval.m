topdir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\';
f_all = dir(topdir);
%% checking extraction and motion corrections
XtotT = cell(length(f_all),1);
i_isfile = false(length(f_all),1);

dims2use = 1:30;
bin_res = 5;
for i = 1:length(f_all)
    %%
    if ~f_all(i).isdir
        i_isfile(i) = true;
        sess = f_all(i).name;
        fname = [f_all(i).folder '/' f_all(i).name];
        load(fname);
        
        figure(i); clf
        
        tiffname = ms.neuron.motion_corr_tiff;
        im_min = ms.neuron.minFrame;
        im_mean = ms.neuron.meanFrame;
        im_max = ms.neuron.maxFrame;
        im_pnr = ms.neuron.pnr_im;
        im_corr = ms.neuron.corr_im;
        spks = (ms.neuron.S_matw);
        craw = normalize_rows(ms.neuron.C+ms.neuron.YrA);
        
        [popcorr, ~, ~] = Fenton_pop_stability(spks, bin_res, ms.timestamps./1000, false);
        
        spd = ms.arena.speed'>=5;
        man = ms.neuron.dataIso(:, dims2use)';
%         [mancorr, ~, ~] = Fenton_pop_stability(normalize_rows(man), bin_res, ms.timestamps./1000, false);
        [mancorr, ~, ~] = Fenton_pop_stability((man), bin_res, ms.timestamps./1000, false);
%         man = man(:, spd);
        man_bin = bin_spks_time(man, 1, ms.timestamps./1000, false);
%         spdbin = bin_spks_time(spd, bin_res, ms.timestamps./1000, false);
%         man_spd = man(1:20, spdbin>5);
%         man_stl = man(:, spdbin<=5);
        XtotT{i} = man_bin';
        subplot_tight(3,3,1);
        imshow(im_min, [0 255])
        subplot_tight(3,3,2);
        imagesc(popcorr); axis image off
%         imagesc(im_pnr); axis image off
        % imshow(im_mean, [0 255])
        subplot_tight(3,3,3);
        imagesc(mancorr); axis image off
%         imagesc(im_corr.*im_pnr); axis image off
        % imshow(im_max, [0 255])
        subplot_tight(3,1,2);
        imagesc(spks)
        subplot_tight(3,1,3);
        imagesc(man_bin)
        drawnow
%         input('')
    else
        i_isfile(i) = false;
    end
end
XtotT = XtotT(i_isfile);
%%
options.kappa       = 0; % default=0;
options.n_t         = 200;%500; % Default: n_t=500. If number of samples per manifold M>500, then set at M.
options.flag_NbyM   = false; % default = false
[output] = manifold_analysis_outputanchor(XtotT, options);