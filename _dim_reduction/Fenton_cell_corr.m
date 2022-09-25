function [cellcorr, cellcorr_prob, cell_sumspks] = Fenton_cell_corr(spks, time_res, time_vec, sliding_method)
%% Compute the binned activity correlation within spike matrix
% INPUTS:
% 'spks' - cell by time activity matrix (M x N)
% 'window_size' - time widow size in samples
% 'sliding_method' - (bool) use sliding window to compute summed spikes; else use fixed bins
% OUTPUT: 
% 'pco_tau'  - Kendall tau for correlation (N x N)
% 'pco_prob' - Kendall tau probability (N x N)
[ncells, ~] = size(spks);

% next calc the correlation between cells i and j
% [cell_sumspks] = bin_spks(spks, cell_window_size, sliding_method);
[cell_sumspks] = bin_spks(spks, time_res, time_vec, sliding_method);
[cellcorr,  cellcorr_prob] = corr(cell_sumspks', 'type', 'Kendall');
cellcorr((eye(ncells)==1)) = NaN;
cellcorr_prob((eye(ncells)==1)) = NaN;
end