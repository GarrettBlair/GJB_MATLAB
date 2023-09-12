function [str] = pfield_split_bayesian_decoding(ms, x, y, spks, pos_bins, cells2use, n_rand)
%%
% spks = ms.spks;
% x = ms.room.x;
% y = ms.room.y;
% spd = ms.arena.speed_smooth;
% pos_bins = params.pos_bins;

t = ms.timestamps./1000;
[spks_bin, ~]        = bin_spks_time(spks,   .25, ms.timestamps./1000, false);
[x_average, ~]       = average_spks_time(x', .25, ms.timestamps./1000, false, 'mean');
[y_average, ~]       = average_spks_time(y', .25, ms.timestamps./1000, false, 'mean');
[time_average, ~]    = average_spks_time(t', .25, ms.timestamps./1000, false, 'mean');
[~,~,~,x_bin, y_bin] = histcounts2(x_average, y_average, pos_bins, pos_bins);
% [spdav, ~]       = average_spks_time(spd', .25, ms.timestamps./1000, false, 'median');

[pfield_splits, ~]   = average_spks_time(ms.room.split_vec', .25, ms.timestamps./1000, false, 'mean');

pfield_splits = double(pfield_splits>=median(pfield_splits)) + 1;

spks_bin(isnan(spks_bin)) = 0;
n_cells_spiking = sum(spks_bin>0,1);


%%
goodsegs = cells2use; % ms.room.split_p<=1 & ms.room.pcell_stats.peakRate>0;
% n_rand = 100;

nbins = length(pos_bins)-1;
x_decoded = NaN(length(x_average),1);
y_decoded = NaN(length(y_average),1);
x_decoded_shuffle = NaN(length(x_average),n_rand);
y_decoded_shuffle = NaN(length(y_average),n_rand);
posterior_probability_map = NaN(nbins, nbins, length(y_average));
for split_loop = 1:length(unique(pfield_splits))
    if split_loop == 1
        placefields = ms.room.pfields_smooth_split2;
    elseif split_loop==2
        placefields = ms.room.pfields_smooth_split1;
    else 
        placefields=[];
    end
%     aaxi = split_loop;
    validinds = find(n_cells_spiking>0 & pfield_splits==split_loop);
    
    spike_matrix = spks_bin(goodsegs, validinds);
    [nsegs,nsamples1] = size(spike_matrix);
    placefields = placefields(goodsegs,:,:);
    
    linear_pfields = reshape(placefields, [nsegs, nbins^2]);
    linear_pfields(isnan(linear_pfields)) = 0;
    
    posterior_prob = linear_pfields'*spike_matrix;
    posterior_prob = normalize_cols(posterior_prob);
    % posterior_prob_submap = posterior_prob;
    posterior_prob_submap = reshape(posterior_prob, [nbins, nbins, nsamples1]);
    [~, pm] = max(posterior_prob, [], 1);
    [pmy, pmx] = ind2sub([nbins, nbins], pm);

    posterior_probability_map(:,:,validinds) = posterior_prob_submap;
    x_decoded(validinds) = pmx;
    y_decoded(validinds) = pmy;
    
    %%%% RANDOM ID SHUFFLE n times
    for j = 1:n_rand
        [~, randord] = sort(rand(nsegs,1));
        spks1_rand = spike_matrix(randord,:);
        
        posterior_prob = linear_pfields'*spks1_rand;
        posterior_prob = normalize_cols(posterior_prob);
        % posterior_prob_submap = posterior_prob;
        % posterior_prob_submap = reshape(posterior_prob, [nbins, nbins, nsamples1]);
        [~, pm] = max(posterior_prob, [], 1);
        [randy, randx] = ind2sub([nbins, nbins], pm);
        
        x_decoded_shuffle(validinds, j) = randx;
        y_decoded_shuffle(validinds, j) = randy;
    end
end
decode_dist = sqrt( (x_decoded-x_bin').^2 + (y_decoded-y_bin').^2);
decode_dist_shuffle = NaN(length(decode_dist), n_rand);
decode_dist_shuffle_median = NaN(n_rand,1);
decode_dist_shuffle_quant25 = NaN(n_rand,1);
decode_dist_shuffle_quant75 = NaN(n_rand,1);
for j = 1:n_rand
    decode_dist_shuffle(:, j) = ( sqrt( (x_decoded_shuffle(:,j)-x_bin').^2 + (y_decoded_shuffle(:,j)-y_bin').^2) ) ;
    d = sqrt( (x_decoded_shuffle(:,j)-x_bin').^2 + (y_decoded_shuffle(:,j)-y_bin').^2);
    decode_dist_shuffle_median(j) = nanmedian( d ) ;
    decode_dist_shuffle_quant25(j) = quantile(d, .25);
    decode_dist_shuffle_quant75(j) = quantile(d, .75);
end

str.decode_dist                     = decode_dist;
str.decode_dist_shuffle_median      = decode_dist_shuffle_median;
str.decode_dist_shuffle_quant25     = decode_dist_shuffle_quant25;
str.decode_dist_shuffle_quant75     = decode_dist_shuffle_quant75;
str.xbin                            = x_bin;
str.ybin                            = y_bin;
str.xdecoded                        = x_decoded;
str.ydecoded                        = y_decoded;
str.posterior_probability_map       = posterior_probability_map;
str.spks_bin                        = spks_bin;
str.x_average                       = x_average;
str.y_average                       = y_average;
str.time_average                    = time_average;
str.pfield_splits                   = pfield_splits;
%%
plotting = false;

if plotting == true
figure(1078); clf;
subplot(1,2,1)
histogram2(decode_dist, n_cells_spiking', [0:20], [0:17], 'Normalization', 'probability', 'FaceColor', 'b')
subplot(1,2,2)
histogram2(nanmean(decode_dist_shuffle,2), n_cells_spiking', [0:20], [0:17], 'Normalization', 'probability', 'FaceColor', 'k')
% histogram2(decode_dist', spdav, [0:20], [0:5:40], 'Normalization', 'probability', 'FaceColor', 'b'); 
% subplot(1,2,2)
% histogram2(decode_dist_shuffle(:)', spdav, [0:20], [0:5:40], 'Normalization', 'probability', 'FaceColor', 'k'); 
figure(1077); clf;
hold on
% histogram(d(spdav<=5), [0:20], 'Normalization', 'probability'); 
% histogram(d(spdav>5), [0:20], 'Normalization', 'probability'); 
histogram(decode_dist_shuffle(:), [0:20], 'Normalization', 'probability', 'FaceColor', 'k')
plot([nanmedian(decode_dist_shuffle_median) nanmedian(decode_dist_shuffle_median)], [0 .2], 'k')
histogram(decode_dist, [0:20], 'Normalization', 'probability', 'FaceColor', 'b'); 
plot([nanmedian(decode_dist) nanmedian(decode_dist)], [0 .2], 'b')
%%
figure(17); clf; axis square
for ii = 5:2:size(spks_bin,2)
    a = squeeze(posterior_probability_map(:,:,ii));
    imagesc(a)
%     imagesc(a, [-.1 .25])
    set(gca, 'YDir', 'normal')
    hold on
    plot(x_bin(ii-4:ii),y_bin(ii-4:ii), 'r-')
    plot(x_bin(ii),y_bin(ii), 'mo')
    plot(x_decoded(ii),y_decoded(ii), 'kx')
    axis square
    drawnow
    pause(.1)
end

end