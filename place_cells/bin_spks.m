function [sumspks] = bin_spks(spks, window_size, sliding_method)

[nsegs, nsamples] = size(spks);
if sliding_method == true
    %%%% Compute summation with sliding bins
    sumspks = NaN(nsegs, nsamples);
    % first calc the number of spikes within the time resolution window time_res
    for i = 1:nsamples
        inds = i-window_size:i+window_size;
        valid = inds>0 & inds<=nsamples;
        sumspks(:, i) = sum(spks(:, inds(valid)), 2);
    end
else
    %%%% Compute summation with fixed bins, no sliding
    nsub = floor(nsamples/window_size)+1;
    tvec = [1:window_size:nsamples, nsamples];
    % first calc the number of spikes within the time resolution window time_res
    if tvec(end)-tvec(end-1) < window_size/4
%         warning('Last bin too small to include, skipping')
        nsub = nsub-1;
    end
    sumspks = NaN(nsegs, nsub);
    for i = 1:nsub-1
        sumspks(:, i) = sum(spks(:, tvec(i):tvec(i+1)), 2);
    end
end