function [vmap_time, vmap_counts] = make_occupancymap_1D(pos_x, timestamps, bins_x)
%%
%%%% INPUT:
% pos_x, pos_y   = x,y values of position samples, should be speed filtered already
% timestamps = values of time in seconds at each position, used to
%                   determine occupancy duration
% bins_x, bins_y = x,y bin edges (or number) used in histcounts2
% minsamples = occupancy inclusion threshold (in time) for spatial bins; 
%          bins occuppied for less than mintime will be excluded from the rate map
%
%%%% OUTPUT:
% vmap_counts = number of samples in a bin
% vmap_counts = ammount of time (seconds) in a bin

if isrow(pos_x)
    pos_x = pos_x';
end
if isrow(timestamps)
    timestamps = timestamps';
end

temp_ts = cat(1, timestamps, timestamps(end));
dt = abs(diff(temp_ts));
if median(dt)>1 % probably sampled in milliseconds
    warning('Variable ''timestamps'' should be given in seconds, and it appears to be smaller')
end

nx_bins = length(bins_x)-1;
[vmap_counts, ~, xbin] = histcounts(pos_x, bins_x);
vmap_time = NaN(nx_bins, 1);
for i=1:nx_bins
    if any(i==xbin)
        currentbin = i==xbin;
        vmap_time(i, 1) = sum(dt(currentbin));
    end
end

