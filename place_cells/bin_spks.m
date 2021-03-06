function [sumspks] = bin_spks(spks, window_size, sliding_method)

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