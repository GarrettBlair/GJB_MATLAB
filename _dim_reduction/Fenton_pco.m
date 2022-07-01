function [pco_tau, pco_prob, sumspks] = Fenton_pco(spks, window_size, sliding_method, corr_method)
%% Compute the binned activity correlation within spike matrix
% INPUTS:
% 'spks' - cell by time activity matrix (M x N)
% 'window_size' - time widow size in samples
% 'sliding_method' - (bool) use sliding window to compute summed spikes; else use fixed bins
% OUTPUT: 
% 'pco_tau'  - Kendall tau for correlation (N x N)
% 'pco_prob' - Kendall tau probability (N x N)

[sumspks] = bin_spks(spks, window_size, sliding_method);
[~, nsamples] = size(sumspks);

% next calc the correlation between time bins i and j
[pco_tau,  pco_prob] = corr(sumspks, 'type', corr_method);
pco_tau((eye(nsamples)==1)) = NaN;
pco_prob((eye(nsamples)==1)) = NaN;
end