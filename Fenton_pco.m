function [pco_tau, pco_prob, sumspks] = Fenton_pco(spks, window_size, sliding_method, corr_method)
%% Compute the binned activity correlation within spike matrix
% INPUTS:
% 'spks' - cell by time activity matrix (M x N)
% 'window_size' - time widow size in samples
% 'sliding_method' - (bool) use sliding window to compute summed spikes; else use fixed bins
% OUTPUT: 
% 'pco_tau'  - Kendall tau for correlation (N x N)
% 'pco_prob' - Kendall tau probability (N x N)

[nsegs, nsamples] = size(spks);
if sliding_method == true
    %%%% Compute summation with sliding bins
    sumspks = NaN(nsegs, nsamples-window_size);
    % first calc the number of spikes within the time resolution window time_res
    for i = 1:nsamples-window_size
        sumspks(:, i) = sum(spks(:, i:i+window_size), 2);
    end
else
    %%%% Compute summation with fixed bins, no sliding
    nsub = floor(nsamples/window_size)-1;
    tvec = 1:window_size:nsamples;
    sumspks = NaN(nsegs, nsub);
    % first calc the number of spikes within the time resolution window time_res
    for i = 1:nsub
        sumspks(:, i) = sum(spks(:, tvec(i):tvec(i+1)), 2);
    end
end
[~, nsamples] = size(sumspks);

% next calc the correlation between time bins i and j
[pco_tau,  pco_prob] = corr(sumspks, 'type', corr_method);
pco_tau((eye(nsamples)==1)) = NaN;
pco_prob((eye(nsamples)==1)) = NaN;
end