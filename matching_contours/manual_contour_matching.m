function manual_contour_matching(top_dir, save_dir, animal_name, data_file_regexp, record_name, sessions2match, ref_session)

if ~exist('sessions2match', 'var')
    sessions2match = [];
end
    
currentDir = pwd;
warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')
% cd(top_dir)
%%
ftot = dir([top_dir '\' data_file_regexp]);
numfiles = length(ftot);%size(aligned_data_struct.spatial_footprints, 2);
ms_file_names = cell(numfiles, 1);
ms_sess_name = cell(numfiles, 1);
filesFolder = ftot(1).folder;
sess_num = 0;
for i = 1:numfiles
    sessFileName = ftot(i).name;
    fname = [filesFolder '\' sessFileName];
    try
        temp = load(fname, 'params');
        sess_num = sess_num+1;
        ms_file_names{sess_num} = fname;
        ms_sess_name{sess_num} = sessFileName;
    catch
        warning('\n%s found, no params and assumed not an analysis file', fname);
    end
end
fprintf('\n%d files found:', sess_num);
ms_file_names = ms_file_names(1:sess_num);
ms_sess_name = ms_sess_name(1:sess_num);
for i = 1:sess_num
    fprintf('\n\t%s', ms_sess_name{sess_num});
end

if ~isempty(sessions2match)
    ms_file_names = ms_file_names(sessions2match);
else
    sessions2match = 1:length(ms_file_names);
end

% ms_file_names = ms_file_names(end-3:end);
% warning('subsampling for test!')

numsess = length(ms_file_names);


% sessNum = sessions2match(1);
cname = sprintf('%s\\%s\\', save_dir, animal_name);
if ~isfolder(cname)
    fprintf('\nSave dir created: \n\t%s\n', cname)
    mkdir(cname);
end
%%
% load(aligned_fname, 'aligned_data_struct');
% numsess = size(aligned_data_struct.spatial_footprints, 2);
% ref_session = aligned_data_struct.reference_session_index;
if ~exist('ref_session', 'var') || isempty(ref_session)
    ref_session = ceil(numsess/2); %aligned_data_struct.reference_session_index;
end
target_indices = setdiff(1:numsess, ref_session);

% bg = imresize(aligned_data_struct.overlapping_FOV, 1.5);

% footprint_projs = aligned_data_struct.footprints_projections;
% footprint_projs_corrected = aligned_data_struct.footprints_projections_corrected;
spatial_footprints = cell(numsess,1); %aligned_data_struct.spatial_footprints;
footprints_projections = cell(numsess,1);
% centroid_locations = aligned_data_struct.centroid_locations;
corrected_spatial_footprints = cell(numsess,1);
corrected_footprints_projections = cell(numsess,1);

down_right_corrections = zeros(numsess, 2);
scaling_corrections = ones(numsess, 1);


% max_raw = cell(1, numsess);
% min_frame = cell(1, numsess);
% 
% max_ne = cell(1, numsess);
% max_ne_neg = cell(1, numsess);
%%
fprintf('BG image - Session: ');
fields2use = {'minFrame', 'meanFrame', 'corr_im', 'pnr_im'};
fieldsOut = {'min_frame', 'max_raw', 'max_ne', 'max_ne_neg'};
bg_size = [1 1]; %imresize(aligned_data_struct.overlapping_FOV, 1.5);
for ii = 1:4
    eval(sprintf('%s = cell(1,numsess);', fieldsOut{ii}));
end
for i = 1:numsess
    fname = ms_file_names{i};
    fprintf('\n\t%s.. ', fname);
    %     fprintf('%d.. ', i);
    temp =  load(fname, 'ms');
    
    for ii = 1:4
        eval(sprintf(' im = temp.ms.neuron.%s;', fields2use{ii}));
        imn = im - quantile(im(:), .01);
        imn = imn./quantile(imn(:), .999);
        eval(sprintf('%s{1,i} = imn;', fieldsOut{ii}));
    end    
%     imne = imn - anisodiff2D(imn, 50, 1/7, 30, 2);
%     imne([1:10, end-10:end], :) = 0;
%     imne(:, [1:10, end-10:end]) = 0;
%     imne = (imne-min(imne(:)))./max(max(imne-min(imne(:))));
%     max_ne{1,i} = imne;%.^2;
%     max_ne_neg{1,i} = (1 - imne).^2;
%     max_ne{1,i} = temp.ms.neuron.corr_im;%.^2;
%     max_ne_neg{1,i} = temp.ms.neuron.pnr_im;
    
    bg_size(1) = max(bg_size(1), temp.ms.neuron.dims(1));
    bg_size(2) = max(bg_size(2), temp.ms.neuron.dims(2));
    
end
bg_size = round(bg_size*1.1);
bg = zeros(bg_size); %imresize(aligned_data_struct.overlapping_FOV, 1.5);
fprintf(' Done! \n');
%%
fprintf('Contours - Session: ');
for i = 1:numsess
    %
    fname = ms_file_names{i};
    fprintf('\n\t%s.. ', fname);
    temp =  load(fname, 'ms');
    
    if isfield(temp.ms.neuron,  'fullA')
        a = temp.ms.neuron.fullA;
    else
        a = temp.ms.neuron.A;
    end
    aa = normalize_cols(a);
    aa(aa<.5) = 0;
    nsegs = size(a,2);
%     all_conts = zeros(nsegs, temp.ms.neuron.dims(1), temp.ms.neuron.dims(2));
%     thresh_conts = zeros(nsegs, temp.ms.neuron.dims(1), temp.ms.neuron.dims(2));
    all_conts = zeros(nsegs, bg_size(1), bg_size(2));
    thresh_conts = false(nsegs, bg_size(1), bg_size(2));
    
    targ_diffs = size(bg) - double(temp.ms.neuron.dims);%size(targ_contour_projection);
    targ_inds = [1+floor(targ_diffs(1)/2), 1+floor(targ_diffs(2)/2)];
    h = temp.ms.neuron.dims(1);
    w = temp.ms.neuron.dims(2);
    for seg = 1:nsegs
%         aa = a(:,seg);
%         aa = a(:,seg);
%         at = aa;
%         at(at<=max(at)*.7)= 0;
        all_conts(seg, targ_inds(1):(targ_inds(1)+h-1), targ_inds(2):(targ_inds(2)+w-1)) = reshape(a(:,seg), [h, w]);
        thresh_conts(seg, targ_inds(1):(targ_inds(1)+h-1), targ_inds(2):(targ_inds(2)+w-1)) = reshape(aa(:,seg)>0, [h, w]);
    end
    bg1 = bg;
    bg1(targ_inds(1):(targ_inds(1)+h-1), targ_inds(2):(targ_inds(2)+w-1)) = max_raw{1,i};
    max_raw{1,i} = bg1;
    
    bg1 = bg;
    bg1(targ_inds(1):(targ_inds(1)+h-1), targ_inds(2):(targ_inds(2)+w-1)) = max_ne{1,i};
    max_ne{1,i} = bg1;
    
    bg1 = bg;
    bg1(targ_inds(1):(targ_inds(1)+h-1), targ_inds(2):(targ_inds(2)+w-1)) = max_ne_neg{1,i};
    max_ne_neg{1,i} = bg1;
    
    spatial_footprints{i} = all_conts;
    footprints_projections{i} = squeeze(sum(thresh_conts,1)>0);
end
fprintf(' Done! \n');

%%
s = 0;
cm_fig = figure;
set(cm_fig, 'Position', [100 100 1539 821]);
set(gcf, 'Name', 'Arrows to move pink, space toggle shift size, q-confirm, r-reset, a/z resize up/down, c change background, v toggle all corrected contours')


ref_contour_projection = footprints_projections{ref_session};
max_contour_inds = [Inf 0 Inf 0];
sess_colors = jet(numsess);
% sess_colors = viridis(numsess);
sess_colors = .5+(sess_colors/2);
raw_option = 2;
resize_step_size = .01;
scatter_pnr = 0;

if ~isempty(record_name)
    v = VideoWriter(record_name);
    v.Quality = 100;
    v.FrameRate = 20;
    v.open;
end
targ_contour_projection = [];
main_denom = numsess/3;%1.2;
other_denom = numsess;
    tc = clock;
for i = 1:numsess
    %%
    starting_guess = true;
    click_count = 0;
    click_shift = false;
    xx = NaN(2,1); yy = NaN(2,1);
    if i ~= ref_session
        target_idx = i;%target_indices(i);
        if  isempty(targ_contour_projection)
            targ_contour_projection = footprints_projections{target_idx};
            corrected_target_contours = footprints_projections{target_idx};
        else
            targ_contour_projection = corrected_target_contours;
        end
        
        if mod(scatter_pnr,2) == 1
                all_proj_r = ref_contour_projection/main_denom; % .*((1 - sess_colors(thissess,1))/main_denom);%/main_denom; % ref session is grey
                all_proj_g = ref_contour_projection/main_denom; % .*((1 - sess_colors(thissess,2))/main_denom);%/main_denom; % ref session is grey
                all_proj_b = ref_contour_projection/main_denom; % .*((1 - sess_colors(thissess,3))/main_denom);%/main_denom; % ref session is grey
                corr_sessns = setdiff(1:i, [ref_session i]);
                for ii2 = 1:length(corr_sessns)
                    thissess = corr_sessns(ii2);
                    all_proj_r = all_proj_r + footprints_projections{thissess}.*(sess_colors(thissess,1)/other_denom); % previous sessions are dull color
                    all_proj_g = all_proj_g + footprints_projections{thissess}.*(sess_colors(thissess,2)/other_denom); % previous sessions are dull color
                    all_proj_b = all_proj_b + footprints_projections{thissess}.*(sess_colors(thissess,3)/other_denom); % previous sessions are dull color
                end
                all_proj_r = all_proj_r + corrected_target_contours.*(sess_colors(target_idx,1)/main_denom); % current session is bright color
                all_proj_g = all_proj_g + corrected_target_contours.*(sess_colors(target_idx,2)/main_denom); % current session is bright color
                all_proj_b = all_proj_b + corrected_target_contours.*(sess_colors(target_idx,3)/main_denom); % current session is bright color
                all_proj = cat(3, all_proj_r, all_proj_g, all_proj_b);
                im = 1.2*all_proj./max(all_proj(:));% make_color_im(bg*0, corrected_target_contours*.2, ref_contour_projection*.2, corrected_target_contours*.2);
        else
            im = make_color_im(bg*0, targ_contour_projection, ref_contour_projection, targ_contour_projection);
        end
%         if starting_guess == true
%             xc = xcorr2(corrected_target_contours*1, 1*ref_contour_projection);
%             [~, mxc] = max(xc, [], 1);
%             [~, mi] = max(xc(:));
%             [mr, mc] = ind2sub([size(xc)], mi);
%             shifth = size(corrected_target_contours,1) - mr;
%             shiftw = size(corrected_target_contours,2) - mc;
%             corrected_target_contours = circshift(corrected_target_contours, shifth, 1);
%             corrected_target_contours = circshift(corrected_target_contours, shiftw, 2);
%             
%             starting_guess = false;
%         end
%         down_right_shift = [shifth shiftw];
%         resize_val = scaling_corrections(target_idx);
        down_right_shift = down_right_corrections(target_idx,:);
        resize_val = scaling_corrections(target_idx);
        
        switch mod(raw_option, 4)
            case 0
                max_image = max_raw;
            case 1
                max_image = max_ne;
            case 2
                max_image = max_ne_neg;
            case 3
                max_image = min_frame;
        end
        targ_im = max_image{1,target_idx};
        raw_im = make_color_im(bg*0, targ_im, max_image{1,ref_session}, targ_im);
        
        button = 1;
        %%
        while button~=113 % q
            %%
            
            im_lims = [find(any(sum(im,3),1), 1)-20, find(any(sum(im,3),1), 1, 'last')+20,...
                find(any(sum(im,3),2), 1)-20, find(any(sum(im,3),2), 1, 'last')+20];
            figure(cm_fig); clf
            subplot_tight(1,2,1); cla
            image(im)
            if click_shift == false && ~isnan(xx(1))
                hold on
                scatter(xx(1), yy(1), 50, 'ro')
                hold off
            end
%             text(im_lims(1)+5, im_lims(3)+5, sprintf('Matching: %s\nREF: %s', ms_file_names{target_idx}(end-30:end), ms_file_names{ref_session}(end-30:end)), 'Color', 'w', 'FontSize', 16);
%             text(im_lims(1)+5, im_lims(3)+5, ['REF - ' ms_sess_name{ref_session}], 'Color', 'w', 'FontSize', 16);
%             text(im_lims(1)+5, im_lims(3)-5, ['TRG - ' ms_sess_name{target_idx}],  'Color', 'w', 'FontSize', 16);
            text(im_lims(1)+5, im_lims(3)+10, ['REF - ' ms_sess_name{ref_session}], 'Color', 'g', 'FontSize', 10);
            text(im_lims(1)+5, im_lims(3)-10, ['TRG - ' ms_sess_name{target_idx}],  'Color', 'm', 'FontSize', 10);
            axis image off
            axis(im_lims)%ms_sess_name{i}
            
            subplot_tight(1,2,2); cla
            image(raw_im);
            axis image off
            axis(im_lims)
            text(im_lims(1)+5, im_lims(3)+5, sprintf('Shifts: [x:%d y:%d]  to  scale: %1.2f',...
                down_right_shift(2), down_right_shift(1), resize_val), 'Color', 'w', 'FontSize', 16);
            if ~isempty(record_name)
                temp = getframe(gcf);
                v.writeVideo(temp.cdata);
            end
            [opx, opy, button] = ginput(1);
            step_size = round(9*abs(sin(s)) + 1);
            if ~isempty(button)
                switch button
                    case {3, 99} % right click OR 'c'
                        raw_option = raw_option + 1;
                        switch mod(raw_option, 4)
                            case 0
                                max_image = max_raw;
                            case 1
                                max_image = max_ne;
                            case 2
                                max_image = max_ne_neg;
                            case 3
                                max_image = min_frame;
                        end
                    case 114 % r, reset shifts
                        down_right_shift = [0 0];
                        resize_val = 1;
%                     case 98 % b, back targ sess
%                         if i > 1
%                             i = i-1
%                             button = 113; % q
%                         end
                    case 28 % left
                        down_right_shift(2) = down_right_shift(2) - step_size;
                    case 29 % right
                        down_right_shift(2) = down_right_shift(2) + step_size;
                    case 30 % up
                        down_right_shift(1) = down_right_shift(1) - step_size;
                    case 31 % down
                        down_right_shift(1) = down_right_shift(1) + step_size;
                    case 32 % space, toggle shift size
                        s = s+.5*pi; % steps in either 20 or 1 picel; toggles
                    case 97 % a, resize up
                        resize_val = resize_val + resize_step_size;
                    case 118 % v, toggle snr scatter
                        scatter_pnr = scatter_pnr+1;
                    case 122 % a, resize down
                        resize_val = resize_val - resize_step_size;
                    case 1 % LClick, clicking to move
                        click_count = click_count+1;
                        click_ind = mod(click_count+1, 2)+1;
                        xx(click_ind) = round(opx);
                        yy(click_ind) = round(opy);
                        click_shift = false;
                        if click_ind == 2
                            click_shift = true;
                            dx = xx(1) - xx(2);
                            dy = yy(1) - yy(2);
                            down_right_shift(2) = down_right_shift(2) + dx;
                            down_right_shift(1) = down_right_shift(1) + dy;
                        end
                end
            end
            down_right_shift = round(down_right_shift);
            
            % RESIZING CONTOURS
            temp_targ_cs = imresize(targ_contour_projection, resize_val, 'nearest');
            if size(temp_targ_cs, 1) > size(bg, 1)
                ad =  abs(size(bg, 1) - size(temp_targ_cs, 1))/2;
                temp_targ_cs = temp_targ_cs(floor(ad)+1:end-ceil(ad), :);
            end
            if size(temp_targ_cs, 2) > size(bg, 2)
                ad =  abs(size(bg, 2) - size(temp_targ_cs, 2))/2;
                temp_targ_cs = temp_targ_cs(:, floor(ad)+1:end-ceil(ad));
            end
            % SHIFTING CONTOURS
            corrected_target_contours = make_color_im(bg*0, temp_targ_cs);
            corrected_target_contours = circshift(squeeze(corrected_target_contours(:,:,1)), down_right_shift);
            
            if mod(scatter_pnr,2) == 1
                all_proj_r = ref_contour_projection/main_denom; % .*((1 - sess_colors(thissess,1))/main_denom);%/main_denom; % ref session is grey
                all_proj_g = ref_contour_projection/main_denom; % .*((1 - sess_colors(thissess,2))/main_denom);%/main_denom; % ref session is grey
                all_proj_b = ref_contour_projection/main_denom; % .*((1 - sess_colors(thissess,3))/main_denom);%/main_denom; % ref session is grey
                corr_sessns = setdiff(1:i, [ref_session i]);
                for ii2 = 1:length(corr_sessns)
                    thissess = corr_sessns(ii2);
                    all_proj_r = all_proj_r + footprints_projections{thissess}.*(sess_colors(thissess,1)/other_denom); % previous sessions are dull color
                    all_proj_g = all_proj_g + footprints_projections{thissess}.*(sess_colors(thissess,2)/other_denom); % previous sessions are dull color
                    all_proj_b = all_proj_b + footprints_projections{thissess}.*(sess_colors(thissess,3)/other_denom); % previous sessions are dull color
                end
                all_proj_r = all_proj_r + corrected_target_contours.*(sess_colors(target_idx,1)/main_denom); % current session is bright color
                all_proj_g = all_proj_g + corrected_target_contours.*(sess_colors(target_idx,2)/main_denom); % current session is bright color
                all_proj_b = all_proj_b + corrected_target_contours.*(sess_colors(target_idx,3)/main_denom); % current session is bright color
                all_proj = cat(3, all_proj_r, all_proj_g, all_proj_b);
                im = 1.2*all_proj./max(all_proj(:));% make_color_im(bg*0, corrected_target_contours*.2, ref_contour_projection*.2, corrected_target_contours*.2);
            else
                im = make_color_im(bg*0, corrected_target_contours, ref_contour_projection, corrected_target_contours);
            end
            % RESIZING IMAGE
            targ_im = imresize(max_image{1,target_idx}, resize_val);
            if size(targ_im, 1) > size(bg, 1)
                ad =  abs(size(bg, 1) - size(targ_im, 1))/2;
                targ_im = targ_im(floor(ad)+1:end-ceil(ad), :);
            end
            if size(targ_im, 2) > size(bg, 2)
                ad =  abs(size(bg, 2) - size(targ_im, 2))/2;
                targ_im = targ_im(:, floor(ad)+1:end-ceil(ad));
            end
            % SHIFTING IMAGE
            corrected_target_image = make_color_im(bg*0, targ_im);
            corrected_target_image = circshift(squeeze(corrected_target_image(:,:,1)), down_right_shift);
            
            raw_im = make_color_im(bg*0, corrected_target_image, max_image{1,ref_session}, corrected_target_image);
        end
        
        down_right_corrections(target_idx,:) = down_right_shift;
        scaling_corrections(target_idx) = resize_val;
        footprints_projections{i} = corrected_target_contours>0;
    else
        target_idx = ref_session;
    end
    stemp = spatial_footprints{target_idx};
    stemp = shiftdim(stemp,1); % make it height x width x cells
    stemp = imresize(stemp, scaling_corrections(target_idx));
    
    if size(stemp, 1) > size(bg, 1)
        ad =  abs(size(bg, 1) - size(stemp, 1))/2;
        stemp = stemp(floor(ad)+1:end-ceil(ad), :, :);
    end
    if size(stemp, 2) > size(bg, 2)
        ad =  abs(size(bg, 2) - size(stemp, 2))/2;
        stemp = stemp(:, floor(ad)+1:end-ceil(ad), :);
    end
    targ_diffs = size(bg) - [size(stemp, 1), size(stemp,2)];
    targ_inds = [1+floor(targ_diffs(1)/2), 1+floor(targ_diffs(2)/2)];
    contour_bg = zeros(size(bg, 1), size(bg, 2), size(stemp,3));
    h = size(stemp,1); w = size(stemp,2);
    contour_bg(targ_inds(1):(targ_inds(1)+h-1), targ_inds(2):(targ_inds(2)+w-1), :) = stemp;
    
    stemp = circshift(contour_bg, [down_right_corrections(i,:), 0]);
    stemp = shiftdim(stemp, 2); % make it cells x height x width
    
    cproj = squeeze(sum(stemp, 1)>0);
    h1 = find(any(cproj,2), 1, 'first');
    h2 = find(any(cproj,2), 1, 'last');
    w1 = find(any(cproj,1), 1, 'first');
    w2 = find(any(cproj,1), 1, 'last');
    max_contour_inds(1) = min(max_contour_inds(1), h1);
    max_contour_inds(2) = max(max_contour_inds(2), h2);
    max_contour_inds(3) = min(max_contour_inds(3), w1);
    max_contour_inds(4) = max(max_contour_inds(4), w2);

    corrected_spatial_footprints{target_idx} = single(stemp); % for Hipp8 or large files
    corrected_footprints_projections{target_idx} = cproj;
    targ_contour_projection = [];

%     sessNum = sessions2match(i);
%     cname = sprintf('%s\\%s\\contours_shifted_linear%d.mat', save_dir, animal_name, sessNum);
%     fprintf('\nSaving: \n\t%s', cname);
%     contours_shifted = corrected_spatial_footprints{i}(:, max_contour_inds(1):max_contour_inds(2),  max_contour_inds(3):max_contour_inds(4));
%     save(cname, 'contours_shifted', '-v7.3');
%     newfname = sprintf('%s\\%s\\%s_manual_match_%d_%d_%d_H%d_M%d.mat', save_dir, animal_name, animal_name, tc(1), tc(2), tc(3), tc(4), tc(5));
%     save(newfname, 'ms_file_names', 'sessions2match', 'down_right_corrections', 'scaling_corrections', '-v7.3');
%     
%     corrected_spatial_footprints{i} = [];
%     spatial_footprints{i} = [];
end
%%
if ~isempty(record_name)
    v.close();
end

% corrected_centroid_locations = compute_centroid_locations(corrected_spatial_footprints', 1);
%

for i = 1:numsess
%     sessNum = sessions2match(i);
    cname = sprintf('%s\\%s\\%s.mat', save_dir, animal_name, ms_sess_name{i});
    fprintf('\nSaving: \n\t%s', cname);
    contours_shifted = corrected_spatial_footprints{i}(:, max_contour_inds(1):max_contour_inds(2),  max_contour_inds(3):max_contour_inds(4));
    save(cname, 'contours_shifted', '-v7.3');
%     corrected_spatial_footprints{i} = [];
%     spatial_footprints{i} = [];
    
end

% tc = clock;
newfname = sprintf('%s\\%s\\%s_manual_match_%d_%d_%d_H%d_M%d.mat', save_dir, animal_name, animal_name, tc(1), tc(2), tc(3), tc(4), tc(5));
fprintf('%s\n', newfname)
% save(newfname, 'ms_file_names', 'sessions2match', 'down_right_corrections', 'scaling_corrections', 'corrected_footprints_projections', 'corrected_centroid_locations', '-v7.3');
save(newfname, 'ms_file_names', 'ms_sess_name', 'sessions2match', 'down_right_corrections', 'scaling_corrections', 'max_contour_inds', '-v7.3');
cd(currentDir);
fprintf('\nDONE!\n')

