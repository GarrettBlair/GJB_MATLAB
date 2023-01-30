function [popcorr, popcorr_prob, pop_sumspks] = Fenton_pop_stability(spks, time_res, time_vec, sliding_method)
%% Compute the binned activity correlation within spike matrix
% INPUTS:
% 'spks' - cell by time activity matrix (M x N)
% 'window_size' - time widow size in samples
% 'sliding_method' - (bool) use sliding window to compute summed spikes; else use fixed bins
% OUTPUT: 
% 'pco_tau'  - Kendall tau for correlation (N x N)
% 'pco_prob' - Kendall tau probability (N x N)

% [pop_sumspks] = bin_spks(spks, pop_window_size, sliding_method);
[pop_sumspks] = bin_spks_time(spks, time_res, time_vec, sliding_method);
[~, nsamples] = size(pop_sumspks);

% next calc the correlation between time bins i and j
[popcorr,  popcorr_prob] = corr(pop_sumspks, 'type', 'Pearson');
popcorr((eye(nsamples)==1)) = NaN;
popcorr_prob((eye(nsamples)==1)) = NaN;
if true
figure; 
subplot(211);
stacked_traces(normalize_rows(pop_sumspks), .5); axis square
colormap viridis
% imagesc(pop_sumspks); axis square
subplot(212);
imagesc(popcorr, [-.3 1]); axis square
drawnow
end
end