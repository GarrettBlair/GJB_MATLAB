function [vmap_time, vmap_counts, xbin, ybin] = make_occupancymap_2D(pos_x, pos_y, signal, bins_x, bins_y)
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
[vmap_counts, ~, ~, ybin, xbin] = histcounts2(pos_y, pos_x, bins_y, bins_x);

vmap_time = NaN(ny_bins, nx_bins);

uy_bins = unique(ybin(ybin>0));
ux_bins = unique(xbin(xbin>0));

xbin(xbin==0) = interp1(find(xbin~=0), xbin(xbin~=0), find(xbin==0), 'nearest', 'extrap');
ybin(ybin==0) = interp1(find(ybin~=0), ybin(ybin~=0), find(ybin==0), 'nearest', 'extrap');

for i=1:length(uy_bins)%1:ny_bins
%     if any(i==ybin) % & j==xbin)
        for j = 1:length(ux_bins)%1:nx_bins
            if any(uy_bins(i)==ybin & ux_bins(j)==xbin)
                currentbin =  uy_bins(i)==ybin & ux_bins(j)==xbin;
                vmap_time(uy_bins(i),ux_bins(j)) = sum(signal(currentbin));
            end
        end
%     end
end

