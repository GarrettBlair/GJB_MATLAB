function [summed_map, vmap_counts, xbin, ybin] = make_summap_2D(pos_x, pos_y, signal, bins_x, bins_y)
%%
%%%% INPUT:
% pos_x, pos_y   = x,y values of position samples, should be speed filtered already
% signal = values to be binned, can be time, spikes, etc.
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
if isrow(pos_y)
    pos_y = pos_y';
end
if isrow(signal)
    signal = signal';
end


nx_bins = length(bins_x)-1;
ny_bins = length(bins_y)-1;
[vmap_counts, ~, ~, xbin, ybin] = histcounts2(pos_x, pos_y, bins_x, bins_y);
summed_map = NaN(ny_bins, nx_bins);
for i=1:nx_bins
    for j = 1:ny_bins
        if any(j==ybin & i==xbin)
            currentbin = i==xbin & j==ybin;
            summed_map(j,i) = sum(signal(currentbin));
        end
    end
end


