function [momentary_pos_info, sub_str, ensemble_prob_min, ensemble_prob_av, p_x] = Fenton_ipos_testing(ms, integration_time, frame_string, params)
% integration_time = .1; % time in seconds to bin activity
% spks_all = normalize_rows(ms.neuron.S_matw);
plotting = false;
tic
max_n_spks = 10;% # bins for spiking

if any(strfind(frame_string, 'polar'))
    binsx = params.yaw_bins;
    binsy = params.rho_bins;
else % euclidean
    binsx = params.pos_bins;
    binsy = params.pos_bins;
end

[sub_str] = bin_session_data(ms, integration_time, frame_string, params);
if any(strfind(frame_string, 'polar'))
%     x_av = sub_str.theta;
%     y_av = sub_str.rho;
    xbin = sub_str.thetabin;
    ybin = sub_str.rhobin;
else
%     x_av = sub_str.x;
%     y_av = sub_str.y;
    xbin = sub_str.xbin;
    ybin = sub_str.ybin;
end
t_av = sub_str.t;
spks_bin = sub_str.spks_bin;
vmap    = sub_str.vmap;

%%
[nsegs, numSamples] = size(spks_bin);
nbinsx = length(binsx)-1;
nbinsy = length(binsy)-1;
momentary_pos_info = NaN(size(spks_bin));


plotting=false;%true
if plotting==true
    figure(99); clf;
    set(gcf, 'Position', [200   359   900   341])
end
for i = 1:nsegs
    %%
    % generate the conditional probability maps for each n spikes
    sm = spks_bin(i, :);
    if any(sm)
    smap_prob = NaN(max_n_spks, nbinsy, nbinsx);
    [smb_prob, sm_edge, s_bin] = histcounts(sm, max_n_spks, 'Normalization', 'probability');
    if plotting==true
        figure(99); clf
    end
    for xx=1:nbinsx
        for yy = 1:nbinsy
            if isnan(vmap(yy,xx))==false % any(j==ybin & i==xbin)
                currentbin = xx==xbin & yy==ybin;
                smap_prob(:, yy, xx) = histcounts(sm(currentbin), max_n_spks, 'Normalization', 'probability');                
            end
        end
    end
    % normalize smap for all position
    smap_prob = smap_prob./ nansum(smap_prob,1);
    for spikeval = 0:max_n_spks-1
        smap = squeeze(smap_prob(spikeval+1,:,:));
%         smap = smap./nansum(smap(:));
%         smap_prob(spikeval+1,:,:) = smap;
            if plotting==true
                figure(99);
                subplot(2,max_n_spks,spikeval+1); cla
%                 imagesc(smap, [-.1*nanmax(smap_prob(:)) nanmax(smap_prob(:))]); % colorbar; axis square
                imagesc(smap, [0 1]); % colorbar; axis square
            end
    end
    % extract the probability of every time sample based on location and
    % spike number
    for j = 1:numSamples
        if ~isnan(vmap(ybin(j), xbin(j)))
            nspks_ind = s_bin(j);
            p_ix = smap_prob(nspks_ind, ybin(j), xbin(j));
            p_i = smb_prob(nspks_ind);
            px = 1; % p_x has a large effect and is already in p_ix 
            momentary_pos_info(i, j) = px * p_ix * log2( p_ix / p_i  );
        end
    end
    end
    if plotting==true % for plotting and checking
        subplot(2,1,2); cla
        yyaxis('left'); cla
        plot(t_av, spks_bin(i,:), 'LineWidth', 2)
        ylim([-.1 1.1])
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
if params.skip_ensemble == false
ensemble_dist_mat = squareform(pdist(spks_bin', 'cosine'));
ensemble_dist_mat(find(eye(numSamples))) = NaN;
ensemble_prob_min = NaN(numSamples, 1);
ensemble_prob_av = NaN(numSamples, 1);
p_x = NaN(numSamples, 1);
current_bin_prob = NaN(nbinsy, nbinsx);

for i=1:nbinsx
    for j = 1:nbinsy
        if isnan(vmap(j,i))==false % any(j==ybin & i==xbin)
            currentbin = i==xbin & j==ybin;
            current_bin_prob(j,i) = sum(currentbin)/numSamples;
            numThisBin = sum(currentbin);
            if numThisBin>2
%                 ensemble_dist_mat_pos = squareform(pdist(spks_bin(:, currentbin)', 'cosine')); % euclidean
                ensemble_dist_mat_pos = squareform(pdist(spks_bin(:, currentbin)', 'cosine'));
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

if isfield(ms, 'goodFrames') && any(~ms.goodFrames)  
    momentary_pos_info(:, ~ms.goodFrames)   = NaN;
    ensemble_prob_av(~ms.goodFrames)        = NaN;
    ensemble_prob_min(~ms.goodFrames)       = NaN;
end

