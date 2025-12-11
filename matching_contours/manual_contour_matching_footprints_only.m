function [contours_shifted, projections_shifted, down_right_corrections, scaling_corrections] = manual_contour_matching_footprints_only(contours_cell, sessions2match, ref_session, nsegs2show, record_name, bg_images)



numsess = length(contours_cell);
if ~exist('sessions2match', 'var') || isempty(sessions2match)
    sessions2match = 1:numsess;
end
if ~exist('nsegs2show', 'var') || isempty(nsegs2show)
    nsegs2show = 0;
end
%%

if ~exist('ref_session', 'var') || isempty(ref_session)
    ref_session = ceil(numsess/2); %aligned_data_struct.reference_session_index;
end
target_indices = setdiff(1:numsess, ref_session);

spatial_footprints = cell(numsess,1); %aligned_data_struct.spatial_footprints;
footprints_projections = cell(numsess,1);
% centroid_locations = aligned_data_struct.centroid_locations;
corrected_spatial_footprints = cell(numsess,1);
corrected_footprints_projections = cell(numsess,1);

down_right_corrections = zeros(numsess, 2);
scaling_corrections = ones(numsess, 2);
%%
bg_size = [size(contours_cell,2) size(contours_cell,3)];
for i = 2:numsess
    
    conts = contours_cell{i};
    [ns, h, w] = size(conts);

    
    bg_size(1) = max(bg_size(1), h);
    bg_size(2) = max(bg_size(2), w);
    
end
bg_size = round(bg_size*1.1);
bg = zeros(bg_size); %imresize(aligned_data_struct.overlapping_FOV, 1.5);
fprintf(' Done! \n');
%%
fprintf('Contours - Session: ');
max_image = cell(1, numsess);
for i = 1:numsess
    %%
    conts = contours_cell{i};
    [nsegs, h, w] = size(conts);
    
    all_conts = zeros(nsegs, bg_size(1), bg_size(2));
%     thresh_conts = zeros(nsegs, bg_size(1), bg_size(2));
    thresh_conts = zeros(nsegs, bg_size(1), bg_size(2));
    
    targ_diffs = size(bg) - [h, w];%size(targ_contour_projection);
    targ_inds = [1+floor(targ_diffs(1)/2), 1+floor(targ_diffs(2)/2)];

    all_conts(:, targ_inds(1):(targ_inds(1)+h-1), targ_inds(2):(targ_inds(2)+w-1)) = conts;
%     all_conts = permute(all_conts, [3 1 2]);
    thresh_conts(:, targ_inds(1):(targ_inds(1)+h-1), targ_inds(2):(targ_inds(2)+w-1)) = conts;
    thresh_conts = uint8(normalize_matrix(thresh_conts)*255);
%     thresh_conts = permute(thresh_conts, [3 1 2]);
    
    spatial_footprints{i} = all_conts;
    if nsegs2show==0
        footprints_projections{i} = squeeze(sum(thresh_conts,1));%>0;
    else
        if nsegs2show>nsegs
            showsegs = floor(nsegs/2);
        else
            showsegs = nsegs2show;
        end
        [~, randord] = sort(rand(nsegs,1));
        im1 = squeeze(sum(thresh_conts(randord(1:showsegs), :, :),1));%>0;
        im2 = squeeze(sum(thresh_conts(randord(showsegs+1:end), :, :),1))./4;%>0;
        footprints_projections{i} = im1+im2;%>0;
    end
    if isempty(bg_images)
        max_image{i} = bg;
    else
        max_image{i} = bg_images{i};
    end
end
fprintf(' Done! \n');

%%
s = 0;
cm_fig = figure;
set(cm_fig, 'Position', [100 100 1539 821]);
set(gcf, 'Name', 'Arrows to move pink, space toggle shift size, [q]Confirm, [r]Reset, [a/z]Resize all, [s/x]Rows, [d/c]Cols, [e]Change bg, [v]Toggle all corrected contours')


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
i = 1;

for i = 1:numsess % while i<numsess % 
    %%
    starting_guess = true;
    click_count = 0;
    click_shift = false;
    xx = NaN(2,1); yy = NaN(2,1);
    main_denom = i+1;%1.2;

    if i ~= ref_session
        target_idx = i;%target_indices(i);
        if  isempty(targ_contour_projection)
            targ_contour_projection = footprints_projections{target_idx};
            corrected_target_contours = footprints_projections{target_idx};
        else
            targ_contour_projection = corrected_target_contours;
        end
                
        down_right_shift = down_right_corrections(target_idx,:);
        resize_val = scaling_corrections(target_idx,:);
%         resize_val_rc = scaling_corrections(target_idx);
        
%         switch mod(raw_option, 4)
%             case 0
%                 max_image = max_raw;
%             case 1
%                 max_image = max_ne;
%             case 2
%                 max_image = max_ne_neg;
%             case 3
%                 max_image = min_frame;
%         end
%         max_image = bg;
        targ_im = bg ; % max_image{1,target_idx};
%         raw_im = make_color_im(bg*0, targ_im, max_image{1,ref_session}, targ_im);
        raw_im = make_color_im(bg*0, targ_im, bg, targ_im);
        im = make_color_im(bg*0, targ_contour_projection, ref_contour_projection, targ_contour_projection);
        
        button = 1;
        %%
        while button~=113 % q
            im_lims = [find(any(sum(im,3),1), 1)-20, find(any(sum(im,3),1), 1, 'last')+20,...
                find(any(sum(im,3),2), 1)-20, find(any(sum(im,3),2), 1, 'last')+20];
            figure(cm_fig); clf
            subplot_tight(1,2,1); cla
            image(im./300)
            if click_shift == false && ~isnan(xx(1))
                hold on
                scatter(xx(1), yy(1), 50, 'ro')
                hold off
            end
            
            text(im_lims(1)+5, im_lims(3)+10, ['REF - ' string(ref_session)], 'Color', 'g', 'FontSize', 10);
            text(im_lims(1)+5, im_lims(3)-10, ['TRG - ' string(target_idx )], 'Color', 'm', 'FontSize', 10);
            axis image off
            axis(im_lims)%ms_sess_name{i}
            
            subplot_tight(1,2,2); cla
            image(raw_im);
            axis image off
            axis(im_lims)
            text(im_lims(1)+5, im_lims(3)+5, sprintf('Shifts: [x:%d y:%d]  to  scale: [r-%1.2f, c-%1.2f]',...
                down_right_shift(2), down_right_shift(1), resize_val(1), resize_val(2)), 'Color', 'w', 'FontSize', 16);
            if ~isempty(record_name)
                temp = getframe(gcf);
                v.writeVideo(temp.cdata);
            end
            [opx, opy, button] = ginput(1);
            step_size = round(9*abs(sin(s)) + 1);
            if ~isempty(button)
                switch button
%                     case {3, 101} % right click OR 'e'
%                         raw_option = raw_option + 1;
%                         switch mod(raw_option, 4)
%                             case 0
%                                 max_image = max_raw;
%                             case 1
%                                 max_image = max_ne;
%                             case 2
%                                 max_image = max_ne_neg;
%                             case 3
%                                 max_image = min_frame;
%                         end
                    case 114 % r, reset shifts
                        down_right_shift = [0 0];
                        resize_val = [1 1];
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
                    case 115 % s, resize rows up
                        resize_val(1) = resize_val(1) + resize_step_size;
                    case 120 % x, resize rows down
                        resize_val(1) = resize_val(1) - resize_step_size;
                    case 100 % d, resize cols up
                        resize_val(2) = resize_val(2) + resize_step_size;
                    case 99 % c, resize cols down
                        resize_val(2) = resize_val(2) - resize_step_size;
                    case 118 % v, toggle snr scatter
                        scatter_pnr = scatter_pnr+1;
                    case 122 % a, resize down
                        resize_val = resize_val - resize_step_size;
                    case 27 % esc - quit
                        [~, ~, button2] = ginput(1);
                        contours_shifted = [];
                        projections_shifted = [];
                        close(gcf)
                        if button2 == 27
                            return;
                        end
%                     case 8 % bkspc, previous contours
%                         if i-1 == ref_session
%                             i = i-3;
%                         else
%                             i = i-2;
%                         end
%                         button = 113;
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
            resize_target = round(size(targ_contour_projection).*resize_val);
            temp_targ_cs = imresize(targ_contour_projection, resize_target, 'nearest');
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
            resize_target = round(size(max_image{1,target_idx}).*resize_val);
            targ_im = imresize(max_image{1,target_idx}, resize_target);
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
        scaling_corrections(target_idx,:) = resize_val;
%         scaling_corrections_RC(target_idx, :) = resize_val;
        footprints_projections{i} = corrected_target_contours>0;
    else
        target_idx = ref_session;
    end
    stemp = spatial_footprints{target_idx};
    stemp = shiftdim(stemp,1); % make it height x width x cells
    resize_target = round([size(stemp,1) size(stemp,2)].*scaling_corrections(target_idx,:));

    stemp = imresize(stemp, resize_target);
    
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
    
    cproj = squeeze(sum(stemp, 1));%>0);
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
    
%     i = i+1

end
%%
if ~isempty(record_name)
    v.close();
end

% corrected_centroid_locations = compute_centroid_locations(corrected_spatial_footprints', 1);
%
contours_shifted = cell(numsess,1);
projections_shifted = cell(numsess,1);
for i = 1:numsess
%     sessNum = sessions2match(i);
%     cname = sprintf('%s\\%s', save_dir, ms_sess_name{i});
%     fprintf('\nSaving: \n\t%s', cname);
    contours_shifted{i} = corrected_spatial_footprints{i}(:, max_contour_inds(1):max_contour_inds(2),  max_contour_inds(3):max_contour_inds(4));
    projections_shifted{i} = corrected_footprints_projections{i}(max_contour_inds(1):max_contour_inds(2),  max_contour_inds(3):max_contour_inds(4));
%     save(cname, 'contours_shifted', '-v7.3');
%     corrected_spatial_footprints{i} = [];
%     spatial_footprints{i} = [];
    
end

% % tc = clock;
% newfname = sprintf('%s\\%s_manual_match_%d_%d_%d_H%d_M%d.mat', save_dir, animal_name, tc(1), tc(2), tc(3), tc(4), tc(5));
% fprintf('%s\n', newfname)
% % save(newfname, 'ms_file_names', 'sessions2match', 'down_right_corrections', 'scaling_corrections', 'corrected_footprints_projections', 'corrected_centroid_locations', '-v7.3');
% save(newfname, 'ms_file_names', 'ms_sess_name', 'sessions2match', 'down_right_corrections', 'scaling_corrections', 'max_contour_inds', '-v7.3');
% cd(currentDir);
% fprintf('\nDONE!\n')

end
% function [im] = make_session_projection(ref_contour_projection, ref_session, target_idx, corrected_target_contours, footprints_projections, main_denom, sess_colors, i)
%     all_proj_r = ref_contour_projection/main_denom; % .*((1 - sess_colors(thissess,1))/main_denom);%/main_denom; % ref session is grey
%     all_proj_g = ref_contour_projection/main_denom; % .*((1 - sess_colors(thissess,2))/main_denom);%/main_denom; % ref session is grey
%     all_proj_b = ref_contour_projection/main_denom; % .*((1 - sess_colors(thissess,3))/main_denom);%/main_denom; % ref session is grey
%     corr_sessns = setdiff(1:i, [ref_session i]);
%     for ii2 = 1:length(corr_sessns)
%         thissess = corr_sessns(ii2);
%         all_proj_r = all_proj_r + footprints_projections{thissess}.*(sess_colors(thissess,1)/other_denom); % previous sessions are dull color
%         all_proj_g = all_proj_g + footprints_projections{thissess}.*(sess_colors(thissess,2)/other_denom); % previous sessions are dull color
%         all_proj_b = all_proj_b + footprints_projections{thissess}.*(sess_colors(thissess,3)/other_denom); % previous sessions are dull color
%     end
%     all_proj_r = all_proj_r + corrected_target_contours.*(sess_colors(target_idx,1)/main_denom); % current session is bright color
%     all_proj_g = all_proj_g + corrected_target_contours.*(sess_colors(target_idx,2)/main_denom); % current session is bright color
%     all_proj_b = all_proj_b + corrected_target_contours.*(sess_colors(target_idx,3)/main_denom); % current session is bright color
%     all_proj = cat(3, all_proj_r, all_proj_g, all_proj_b);
%     im = 1.2*all_proj./max(all_proj(:)).*300;% make_color_im(bg*0, corrected_target_contours*.2, ref_contour_projection*.2, corrected_target_contours*.2);
% end
