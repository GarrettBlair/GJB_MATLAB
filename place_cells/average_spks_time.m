function [averagespks, grouped_time] = average_spks_time(spks, time_res, time_vec, sliding_method, average_method)

[nsegs, nsamples] = size(spks);
[d1, d2] = size(time_vec);
if nsegs>nsamples
%     warning('Input dimensions appear wrong; expected segs x time, but time>segs')
%     spks = spks';
%     [nsegs, nsamples] = size(spks);
end
% if d2 ~= nsamples && d1 == nsamples 
%     time_vec = time_vec';
% end
if sliding_method == true
    %%%% Compute summation with sliding bins
    averagespks = NaN(nsegs, nsamples);
    % first calc the number of spikes within the time resolution window time_res
    for i = 1:nsamples
        a1 = time_vec(i) - round(time_res/2);
        a2 = time_vec(i) + round(time_res/2);
        ind = time_vec>a1 & time_vec<=a2;
        averagespks(:, i) = sum(spks(:, ind), 2);
    end
else
    %%%% Compute summation with fixed bins, no sliding
    time_vec = time_vec - time_vec(1);
    nt = ceil(time_vec(end)/time_res);
    grouped_time = time_vec*NaN;
    for i = 1:nt
        inds = find(time_vec>=(i-1)*time_res & time_vec<(i)*time_res); 
        grouped_time(inds) = i;
    end

%     try
%     grouped_time = [0; cumsum(diff(mod(time_vec, time_res))<=0)];
%     catch
%     grouped_time = [0  cumsum(diff(mod(time_vec, time_res))<=0)];
%     end
    nsub = nt; %max(grouped_time); % floor(nsamples/window_size)+1;
    averagespks = NaN(nsegs, nsub);
    for i = 1:nsub
        switch average_method
            case 'sum'
                averagespks(:, i) = sum(spks(:, grouped_time==i), 2);
            case 'mean'
                averagespks(:, i) = mean(spks(:, grouped_time==i), 2);
            case 'average'
                averagespks(:, i) = mean(spks(:, grouped_time==i), 2);
            case 'median'
                averagespks(:, i) = median(spks(:, grouped_time==i), 2);
            case 'mode'
                averagespks(:, i) = mode(spks(:, grouped_time==i), 2);
            otherwise
                error('Uknown averaging method')
        end
    end
end