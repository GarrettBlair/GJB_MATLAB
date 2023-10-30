function [momentary_pos_info, sub_str, ensemble_prob_min, ensemble_prob_av] = Fenton_ipos(ms, integration_time, frame_string, params)
% integration_time = .1; % time in seconds to bin activity
% spks_all = normalize_rows(ms.neuron.S_matw);
plotting = false;
tic
if any(strfind(frame_string, 'polar'))
    binsx = params.theta_bins;
    binsy = params.rho_bins;
    if strcmp(frame_string, 'room')
        xo = ms.room.x;
        yo = ms.room.y;
    elseif strcmp(frame_string, 'arena')
        xo = ms.arena.x;
        yo = ms.arena.y;
    end
    [xo, yo] = cart2pol(xo,yo);
else % euclidean
    binsx = params.pos_bins;
    binsy = params.pos_bins;
    if strcmp(frame_string, 'room')
        xo = ms.room.x;
        yo = ms.room.y;
    elseif strcmp(frame_string, 'arena')
        xo = ms.arena.x;
        yo = ms.arena.y;
    end
end
dt = ms.dt;
to = ms.timestamps./1000;

spks = ms.spks>0;
[spks_bin, group] = bin_spks_time(spks, integration_time, ms.timestamps./1000, false);

spks_bin = round(spks_bin);
max_n_spks1 = 0;%ceil(size(spks, 2)/size(spks_bin,2));
max_n_spks2 = max(spks_bin(:));
max_n_spks = max(max_n_spks1, max_n_spks2);
% [~, speed_epochs] = get_speed_epochs(ms.arena.speed_smooth, params);

% params.pos_bins = [-45, -36:4:36, 45]; % in cm, x and y
if isempty(params) % use default
params.pos_bins = [-45, -36:4:36, 45]; % in cm, x and y
params.occupancy_thresh = integration_time; % min in seconds
params.skip_ensemble = true;
end
% ksize = params.pfield_kernel_radius;
occupancy_thresh = params.occupancy_thresh;


x_av = average_spks_time(xo', integration_time, ms.timestamps./1000, false, 'mean');
y_av = average_spks_time(yo', integration_time, ms.timestamps./1000, false, 'mean');
t_av = average_spks_time(to', integration_time, ms.timestamps./1000, false, 'mean');
[dt_av, grouped_time] = bin_spks_time(dt', integration_time, ms.timestamps./1000, false);

% ipos_timestamps = t_av;

% nbins = length(bins)-1;
nbinsx = length(binsx)-1;
nbinsy = length(binsy)-1;

[vmap, countmap, xbin, ybin]        = make_occupancymap_2D(x_av, y_av, dt_av, binsx, binsy);
% [vmap_full, countmap_full, xb, yb]  = make_occupancymap_2D(xo, yo, dt, binsx, binsy);
% xbin(xbin==0) = interp1(find(xbin~=0), xbin(xbin~=0), find(xbin==0), 'nearest', 'extrap');
% ybin(ybin==0) = interp1(find(ybin~=0), ybin(ybin~=0), find(ybin==0), 'nearest', 'extrap');
% [vmap_o, ~, xbino, ybino]   = make_occupancymap_2D(xo, yo, dt, bins, bins);
[nsegs, nsamples] = size(spks_bin);

vmap(countmap<1) = NaN;
vmap(vmap<occupancy_thresh) = NaN;
prob_vmap = countmap;
prob_vmap = prob_vmap./nansum(prob_vmap(:));

sub_str = [];
sub_str.x       = x_av';
sub_str.y       = y_av';
sub_str.xbin    = xbin;
sub_str.ybin    = ybin;
sub_str.t       = t_av';
sub_str.dt      = dt_av';
sub_str.spks_bin= spks_bin;
sub_str.vmap    = vmap;
sub_str.spks_bin_group    = group;
%%
momentary_pos_info = NaN(size(spks_bin));
valid_frames = 1:size(spks_bin,2); % ms.good_frames


plotting=false;%true
if plotting==true
    figure(99); clf;
    set(gcf, 'Position', [200   359   900   341])
end
% 93

for i = 1:nsegs
    %%
    % generate the conditional probability maps for each n spikes
    smap_prob = NaN(max_n_spks+1, nbinsy, nbinsx);
    smb_prob = NaN(max_n_spks+1, 1);
    if plotting==true
        figure(99); clf
    end
    for spikeval = 0:max_n_spks
        sm = spks_bin(i, :);
%         [smap] = make_summap_2D(x_av,  y_av,  sm, bins, bins);
%         [mi, ir] = infoSkaggs(smap, vmap)
        if any(sm==spikeval)
            smb = sm==spikeval;%(sm~=ii) = NaN;
            smb_prob(spikeval+1) = sum(smb)/length(smb);
            [smap] = make_probabilitymap_2D(x_av,  y_av,  smb, binsx, binsy);
            smap(isnan(vmap)) = NaN;
            smap_prob(spikeval+1,:,:) = smap;

        elseif plotting==true
            subplot(2,max_n_spks+1,spikeval+1); cla
            axis off
        end
    end
    % normalize smap for all position
    smap_prob = smap_prob./ nansum(smap_prob,1);
    for spikeval = 0:max_n_spks
        smap = squeeze(smap_prob(spikeval+1,:,:));
%         smap = smap./nansum(smap(:));
%         smap_prob(spikeval+1,:,:) = smap;
            if plotting==true
                figure(99);
                subplot(2,max_n_spks+1,spikeval+1); cla
%                 imagesc(smap, [-.1*nanmax(smap_prob(:)) nanmax(smap_prob(:))]); % colorbar; axis square
                imagesc(smap, [0 1]); % colorbar; axis square
            end
    end
    % extract the probability of every time sample based on location and
    % spike number
    for j = 1:length(xbin)
        nspks_ind = spks_bin(i,j) + 1;
        p_ix = smap_prob(nspks_ind, ybin(j), xbin(j));
        ttt = smap_prob(:, ybin(j), xbin(j));
        ttt= ttt/nansum(ttt);
        p_ix = ttt(nspks_ind);
        p_i = smb_prob(nspks_ind);
        px = prob_vmap(ybin(j), xbin(j));
        momentary_pos_info(i, j) = px * p_ix * log2( p_ix / p_i  );
    end
    if plotting==true % for plotting and checking
        subplot(2,1,2); cla
        yyaxis('left'); cla
        plot(t_av, spks_bin(i,:), 'LineWidth', 2)
        ylim([-1 max_n_spks+1])
        hold on
        yyaxis('right'); cla
        plot(t_av, momentary_pos_info(i,:))
%         ylim([-1 1.2+max(momentary_pos_info(i,:))])
        xlim([-5 t_av(end)+5])
        drawnow
%         input('')
    end
end
%%
plotting = false;
% ensemble_dist_mat = squareform(pdist(spks_bin', 'cosine'));
% ensemble_dist_mat = squareform(pdist(spks_bin', 'mahalanobis'));
ensemble_dist_mat = squareform(pdist(spks_bin', 'euclidean'));
ensemble_dist_mat(find(eye(nsamples))) = NaN;
ensemble_prob_min = NaN(nsamples, 1);
ensemble_prob_av = NaN(nsamples, 1);
if params.skip_ensemble ~= true
current_bin_prob = NaN(nbinsy, nbinsx);

for i=1:nbinsx
    for j = 1:nbinsy
        if isnan(vmap(j,i))==false % any(j==ybin & i==xbin)
            currentbin = i==xbin & j==ybin;
            current_bin_prob(j,i) = sum(currentbin)/length(currentbin);
            if sum(currentbin)>2
%                 ensemble_dist_mat_pos = squareform(pdist(spks_bin(:, currentbin)', 'cosine')); % euclidean
                ensemble_dist_mat_pos = squareform(pdist(spks_bin(:, currentbin)', 'euclidean')); % euclidean
                ensemble_dist_mat_pos(find(eye(sum(currentbin)))) = NaN;
                binInd = find(currentbin);
                
                for binLoop = 1:sum(currentbin)
                    alldists = ensemble_dist_mat(binInd(binLoop),:);
                    alldists = alldists(setdiff(1:size(spks_bin, 2), binInd(binLoop)));
                    min_dist_all = nanmin(alldists);
                    av_dist_all  = nanmedian(alldists);
                    
                    thisdists = ensemble_dist_mat_pos(binLoop,:);
                    thisdists = thisdists(setdiff(1:sum(currentbin), binLoop));
                    
                    min_dist_thisbin = nanmin(thisdists);
                    av_dist_thisbin  = nanmedian(thisdists);
                    
                    px = current_bin_prob(j,i);
                    ps_x = sum(min_dist_thisbin<=ensemble_dist_mat_pos(:)) ./ (sum(currentbin)^2);
                    ps = sum(min_dist_all<=ensemble_dist_mat(:)) ./ (size(spks_bin, 2)^2);
                    ensemble_prob_min(binInd(binLoop)) = px * ps_x * log2( ps_x / ps );
                    
                    ps_x = sum(av_dist_thisbin<=ensemble_dist_mat_pos(:)) ./ (sum(currentbin)^2);
                    ps = sum(av_dist_all<=ensemble_dist_mat(:)) ./ (size(spks_bin, 2)^2);
                    ensemble_prob_av(binInd(binLoop))  = px * ps_x * log2( ps_x / ps );
                end
            end
        end
    end
end

end
% figure(4); clf; hold on;
% yyaxis('left')
% % plot(nanmean(momentary_pos_info,1));
% mtrace = nansum(abs(momentary_pos_info),1);
% mtrace = nanmean((momentary_pos_info),1);
% plot(mtrace);
% yyaxis('right')
% plot(ensemble_prob)
% 
% 
% figure; plot3(xbin, ybin, mtrace);
% figure; plot3(xbin, ybin, ensemble_prob);
% % drawnow
