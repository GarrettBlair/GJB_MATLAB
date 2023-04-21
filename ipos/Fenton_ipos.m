function [momentary_pos_info, sub_str, ensemble_prob_x] = Fenton_ipos(ms, integration_time, frame_string, params)
% integration_time = .1; % time in seconds to bin activity
% spks_all = normalize_rows(ms.neuron.S_matw);
plotting = false;
tic
if strcmp(frame_string, 'room')
    xo = ms.room.x;
    yo = ms.room.y;
elseif strcmp(frame_string, 'arena')
    xo = ms.arena.x;
    yo = ms.arena.y;
end
dt = ms.dt;
to = ms.timestamps./1000;

[spks_bin, group] = bin_spks_time(ms.neuron.S_matw>0, integration_time, ms.timestamps./1000, false);
max_n_spks1 = ceil(size(ms.neuron.S_matw,2)/size(spks_bin,2));
max_n_spks2 = max(spks_bin(:));
max_n_spks = max(max_n_spks1, max_n_spks2);
% [~, speed_epochs] = get_speed_epochs(ms.arena.speed_smooth, params);

% params.pos_bins = [-45, -36:4:36, 45]; % in cm, x and y
if isempty(params) % use default
params.pos_bins = [-45, -36:4:36, 45]; % in cm, x and y
params.occupancy_thresh = integration_time; % min in seconds
end
bins = params.pos_bins;
% ksize = params.pfield_kernel_radius;
occupancy_thresh = params.occupancy_thresh;


x_av = average_spks_time(xo', integration_time, ms.timestamps./1000, false, 'mean');
y_av = average_spks_time(yo', integration_time, ms.timestamps./1000, false, 'mean');
t_av = average_spks_time(to', integration_time, ms.timestamps./1000, false, 'mean');
dt_av = bin_spks_time(dt', integration_time, ms.timestamps./1000, false);

% ipos_timestamps = t_av;

nbins = length(bins)-1;
[vmap, countmap, xbin, ybin]   = make_occupancymap_2D(x_av, y_av, dt_av, bins, bins);
% [vmap_o, ~, xbino, ybino]   = make_occupancymap_2D(xo, yo, dt, bins, bins);
[nsegs, nsamples] = size(spks_bin);

vmap(countmap<1) = NaN;
vmap(vmap<occupancy_thresh) = NaN;

sub_str = [];
sub_str.x       = x_av;
sub_str.y       = y_av;
sub_str.xbin    = xbin;
sub_str.ybin    = ybin;
sub_str.t       = t_av;
sub_str.dt      = dt_av;
sub_str.spks_bin= spks_bin;
sub_str.vmap    = vmap;
sub_str.spks_bin_group    = group;
%%
momentary_pos_info = NaN(size(spks_bin));
if plotting==true
    figure(99); clf;
    set(gcf, 'Position', [200   359   900   341])
end
% 93


for i = 1:nsegs
    %%
    % generate the conditional probability maps for each n spikes
    smap_prob = NaN(max_n_spks+1, nbins, nbins);
    smb_prob = NaN(max_n_spks+1, 1);
    for spikeval = 0:max_n_spks
        sm = spks_bin(i, :);
%         [smap] = make_summap_2D(x_av,  y_av,  sm, bins, bins);
%         [mi, ir] = infoSkaggs(smap, vmap)
        if any(sm==spikeval)
            smb = sm==spikeval;%(sm~=ii) = NaN;
            smb_prob(spikeval+1) = sum(smb)/length(smb);
            [smap] = make_probabilitymap_2D(x_av,  y_av,  smb, bins, bins);
            smap(isnan(vmap)) = NaN;
            smap_prob(spikeval+1,:,:) = smap;
            if plotting==true
                figure(99);
                subplot(2,max_n_spks+1,spikeval+1); cla
                imagesc(smap, [0 1]); % colorbar; axis square
%                 surf(bins(1:end-1), bins(1:end-1), smap)
            end
        elseif plotting
            subplot(2,max_n_spks+1,spikeval+1); cla
            axis off
        end
    end
    % extract the probability of every time sample based on location and
    % spike number
    for j = 1:length(xbin)
        nspks_ind = spks_bin(i,j) + 1;
        p_ix = smap_prob(nspks_ind, ybin(j), xbin(j));
        momentary_pos_info(i, j) = p_ix * log2( p_ix/smb_prob(nspks_ind)  );
    end
    if plotting==true % for plotting and checking
        subplot(2,1,2); cla
        yyaxis('left'); cla
        plot(t_av, spks_bin(i,:))
        ylim([-1 max_n_spks+1])
        hold on
        yyaxis('right'); cla
        plot(t_av, momentary_pos_info(i,:), 'LineWidth', 2)
        ylim([-1 1.2+max(momentary_pos_info(:))])
        xlim([-5 t_av(end)+5])
        drawnow
%         input('')
    end
end
%%
plotting = false;

ensemble_dist_mat = squareform(pdist(spks_bin', 'euclidean'));
% ensemble_dist_mat = squareform(pdist(spks_bin', 'mahalanobis'));
% ensemble_dist_mat = squareform(pdist(spks_bin', 'euclidean'));
ensemble_dist_mat(find(eye(nsamples))) = NaN;
ensemble_prob = NaN(nsamples, 1);
ensemble_prob_x = NaN(nsamples, 1);
nbins = length(bins)-1;
current_bin_prob = NaN(nbins, nbins);
% prob_map = NaN(nbins, nbins);
for i=1:nbins
    for j = 1:nbins
        if isnan(vmap(j,i))==0 % any(j==ybin & i==xbin)
            currentbin = i==xbin & j==ybin;
            current_bin_prob(j,i) = sum(currentbin)/length(currentbin);
            if sum(currentbin)>2
                
                ensemble_dist_mat_pos = squareform(pdist(spks_bin(:, currentbin)', 'euclidean'));
                ensemble_dist_mat_pos(find(eye(sum(currentbin)))) = NaN;
                binInd = find(currentbin);
                
                for binLoop = 1:sum(currentbin)
                    thisdists = ensemble_dist_mat_pos(binLoop,:);
                    thisdists = thisdists(setdiff(1:sum(currentbin), binLoop));
                    min_dist_thisbin = min(thisdists);
                    
                    ps_x = sum(min_dist_thisbin>=ensemble_dist_mat_pos(:)) ./ (sum(currentbin)^2);
                    
                    px = current_bin_prob(j,i);
                    
                    alldists = ensemble_dist_mat(binInd(binLoop),:);
                    alldists = alldists(setdiff(1:size(spks_bin, 2), binInd(binLoop)));
                    min_dist_all = min(alldists);
                    ps = sum(min_dist_all>=ensemble_dist_mat(:)) ./ (size(spks_bin, 2)^2);
                    ensemble_prob(binInd(binLoop)) = ps_x*log2(ps_x/ps);
                    ensemble_prob_x(binInd(binLoop)) = px*ps_x*log2(ps_x/ps);
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
