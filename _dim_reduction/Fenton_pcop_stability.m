function [popcorr, popcorr_prob, pop_sumspks] = Fenton_pcop_stability(spks, pop_window_size, sliding_method)
%% Compute the binned activity correlation within spike matrix
% INPUTS:
% 'spks' - cell by time activity matrix (M x N)
% 'window_size' - time widow size in samples
% 'sliding_method' - (bool) use sliding window to compute summed spikes; else use fixed bins
% OUTPUT: 
% 'pco_tau'  - Kendall tau for correlation (N x N)
% 'pco_prob' - Kendall tau probability (N x N)

[pop_sumspks] = bin_spks(spks, pop_window_size, sliding_method);
[ncells, nsamples] = size(pop_sumspks);

% next calc the correlation between time bins i and j
popcorr = NaN(nsamples);
popcorr_prob = NaN(nsamples);
[popcorr,  popcorr_prob] = corr(pop_sumspks, 'type', 'Pearson');
popcorr((eye(nsamples)==1)) = NaN;
popcorr_prob((eye(nsamples)==1)) = NaN;

% % next calc the correlation between cells i and j
% [cell_sumspks] = bin_spks(spks, cell_window_size, sliding_method);
% cellcorr = NaN(ncells);
% cellcorr_prob = NaN(ncells);
% [cellcorr,  cellcorr_prob] = corr(cell_sumspks', 'type', 'Kendall');
% cellcorr((eye(ncells)==1)) = NaN;
% cellcorr_prob((eye(ncells)==1)) = NaN;
end