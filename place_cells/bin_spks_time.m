function [sumspks] = bin_spks_time(spks, time_res, time_vec, sliding_method)

[nsegs, nsamples] = size(spks);
if sliding_method == true
    %%%% Compute summation with sliding bins
    sumspks = NaN(nsegs, nsamples);
    % first calc the number of spikes within the time resolution window time_res
    for i = 1:nsamples
        a1 = time_vec(i) - round(time_res/2);
        a2 = time_vec(i) + round(time_res/2);
        ind = time_vec>a1 & time_vec<=a2;
        sumspks(:, i) = sum(spks(:, ind), 2);
    end
else
    %%%% Compute summation with fixed bins, no sliding
    grouped_time = cumsum(diff(mod(time_vec, time_res))<0);
    nsub = max(grouped_time); % floor(nsamples/window_size)+1;
    % first calc the number of spikes within the time resolution window time_res
    last_group = find(nsub == grouped_time);
    last_group_duration = time_vec(last_group(end)) - time_vec(last_group(1));
    % Check how long the final bin is. If < half the time res, omit it
    if last_group_duration < time_res/2 % tvec(end)-tvec(end-1) < window_size/2
        %         warning('Last bin too small to include, skipping')
        nsub = nsub-1;
    end
    sumspks = NaN(nsegs, nsub);
    for i = 1:nsub
        sumspks(:, i) = sum(spks(:, grouped_time==i), 2);
    end
end