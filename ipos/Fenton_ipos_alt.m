function [momentary_pos_info, sub_str, ensemble_prob_min, ensemble_prob_av, p_x] = Fenton_ipos(ms, integration_time, frame_string, params)
% integration_time = .1; % time in seconds to bin activity
% spks_all = normalize_rows(ms.neuron.S_matw);
plotting = false;
tic

if any(strfind(frame_string, 'polar'))
    binsx = params.yaw_bins;
    binsy = params.rho_bins;
    if contains(frame_string, 'room')
        xo = ms.room.x;
        yo = ms.room.y;
    elseif contains(frame_string, 'arena')
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
spks = ms.spks>0;
dt = ms.dt;
to = ms.timestamps./1000;
if isfield(ms, 'goodFrames') && any(~ms.goodFrames)
    if (length(xo)~=sum(ms.goodFrames) || size(spks,2)~=sum(ms.goodFrames))
    fprintf('\n Ipos - downsample frames to good frames \n\t\t\t(removing %d) \n', sum(~ms.goodFrames))
    f = ms.goodFrames;
    xo = xo(f); yo = yo(f); to = to(f);
    spks = spks(:, f);
    dt = dt(f); 
    end    
%     valid_frames = 1:size(spks,2); % ms.good_frames
end
baddt = (dt./median(dt)) >= 3;
dt(baddt) = 3*median(dt);

% [spks_bin, group] = bin_spks_time(spks, integration_time, to, false);
[spks_bin, group] = average_spks_time(spks, integration_time, to, false, 'sum');

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


x_av = average_spks_time(xo', integration_time, to, false, 'mean');
y_av = average_spks_time(yo', integration_time, to, false, 'mean');
t_av = average_spks_time(to', integration_time, to, false, 'mean');
% [dt_av, grouped_time] = bin_spks_time(dt', integration_time, to, false);
[dt_av, ~] = average_spks_time(dt', integration_time, to, false, 'sum');

% ipos_timestamps = t_av;

% nbins = length(bins)-1;
nbinsx = length(binsx)-1;
nbinsy = length(binsy)-1;

[vmap, countmap, xbin, ybin]        = make_occupancymap_2D(x_av, y_av, dt_av, binsx, binsy);
% [vmap_full, countmap_full, xb, yb]  = make_occupancymap_2D(xo, yo, dt, binsx, binsy);
% xbin(xbin==0) = interp1(find(xbin~=0), xbin(xbin~=0), find(xbin==0), 'nearest', 'extrap');
% ybin(ybin==0) = interp1(find(ybin~=0), ybin(ybin~=0), find(ybin==0), 'nearest', 'extrap');
% [vmap_o, ~, xbino, ybino]   = make_occupancymap_2D(xo, yo, dt, bins, bins);
[nsegs, numSamples] = size(spks_bin);

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

plotting=false;%true
if plotting==true
    figure(99); clf;
    set(gcf, 'Position', [200   359   900   341])
end
for i = 1:nsegs
    %%
    % generate the conditional probability maps for each n spikes
    if plotting==true
        figure(99); clf
    end
    sm = spks_bin(i, :);
    for spikeval = 0:max_n_spks
        smb = sm==spikeval;%(sm~=ii) = NaN;
        smb_idx = find(smb);
        p_i  = mean(smb); % overall prob
        for j = 1:length(smb_idx)
            pos_idx = find(ybin == ybin(smb_idx(j)) & xbin == xbin(smb_idx(j)));
            smb_pos = sm(pos_idx) == spikeval;
            p_ix = mean(smb_pos); % prob at position
            momentary_pos_info(i, smb_idx(j)) = p_ix * log2( p_ix / p_i  );
        end
    end
    if any(isnan(momentary_pos_info(i,:)))
        figure;
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
ensemble_dist_mat(find(eye(numSamples))) = NaN;
ensemble_prob_min = NaN(numSamples, 1);
ensemble_prob_av = NaN(numSamples, 1);
p_x = NaN(numSamples, 1);
if params.skip_ensemble ~= true
current_bin_prob = NaN(nbinsy, nbinsx);

for i=1:nbinsx
    for j = 1:nbinsy
        if isnan(vmap(j,i))==false % any(j==ybin & i==xbin)
            currentbin = i==xbin & j==ybin;
            current_bin_prob(j,i) = sum(currentbin)/numSamples;
            numThisBin = sum(currentbin);
            if numThisBin>2
%                 ensemble_dist_mat_pos = squareform(pdist(spks_bin(:, currentbin)', 'cosine')); % euclidean
                ensemble_dist_mat_pos = squareform(pdist(spks_bin(:, currentbin)', 'euclidean')); % euclidean
                ensemble_dist_mat_pos(find(eye(numThisBin))) = NaN;
                binInd = find(currentbin);
                
                for binLoop = 1:numThisBin
                    alldists = ensemble_dist_mat(binInd(binLoop),:);
                    alldists = alldists(setdiff(1:numSamples, binInd(binLoop)));
                    min_dist_all = nanmin(alldists);
                    av_dist_all  = nanmedian(alldists);
                    
                    thisdists = ensemble_dist_mat_pos(binLoop,:);
                    thisdists = thisdists(setdiff(1:numThisBin, binLoop));
                    
                    min_dist_thisbin = nanmin(thisdists);
                    av_dist_thisbin  = nanmedian(thisdists);
                    
                    px = 1; % current_bin_prob(j,i);
                    ps_x = sum(min_dist_thisbin<=ensemble_dist_mat_pos(:)) ./ (numThisBin^2);
                    ps = sum(min_dist_all<=ensemble_dist_mat(:)) ./ (numSamples^2);
                    ensemble_prob_min(binInd(binLoop)) = px * ps_x * log2( ps_x / ps );
                    
                    ps_x = sum(av_dist_thisbin<=ensemble_dist_mat_pos(:)) ./ (numThisBin^2);
                    ps = sum(av_dist_all<=ensemble_dist_mat(:)) ./ (numSamples^2);
                    ensemble_prob_av(binInd(binLoop))  = px * ps_x * log2( ps_x / ps );
                    
                    p_x(binInd(binLoop)) = px;
                end
            end
        end
    end
end

end

