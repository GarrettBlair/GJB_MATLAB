function [sfep_map, map_counts, xbin, ybin] = sfep_prob_map_2D(pos_x, pos_y, sfep_signal, bins_x, bins_y)
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
if isrow(sfep_signal)
    sfep_signal = sfep_signal';
end


nx_bins = length(bins_x)-1;
ny_bins = length(bins_y)-1;
[map_counts, ~, ~, ybin, xbin] = histcounts2(pos_y, pos_x, bins_y, bins_x);

sfep_map = NaN(ny_bins, nx_bins);

uy_bins = unique(ybin(ybin>0));
ux_bins = unique(xbin(xbin>0));

xbin(xbin==0) = interp1(find(xbin~=0), xbin(xbin~=0), find(xbin==0), 'nearest', 'extrap');
ybin(ybin==0) = interp1(find(ybin~=0), ybin(ybin~=0), find(ybin==0), 'nearest', 'extrap');

for i=1:length(uy_bins)
    for j = 1:length(ux_bins)
        if any(uy_bins(i)==ybin & ux_bins(j)==xbin)
            currentbin =  uy_bins(i)==ybin & ux_bins(j)==xbin;
            a = sum(sfep_signal(currentbin)>0);
            b = sum(sfep_signal(currentbin)<0);
            prob_room = a / (a+b);
            sfep_map(uy_bins(i),ux_bins(j)) = prob_room;
        end
    end
end

