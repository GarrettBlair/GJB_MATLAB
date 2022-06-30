function [p] = place_cell_stats(spks, pfields, spike_maps, occupancy_map)
%%
[nsegs, nframes] = size(spks);
if nsegs ~= size(pfields,1)
    error('Dim mismatch!')
end
map_dims = ndims(spike_maps);
if map_dims ~= 2 && map_dims ~= 3
    error('Unkown map type!')
end
infoPerSpike    = NaN(nsegs, 1);
infoRate        = NaN(nsegs, 1);
spkRate         = NaN(nsegs, 1);
peakRate        = NaN(nsegs, 1);
for i = 1:nsegs
    switch map_dims
        case 2
            smap = squeeze(spike_maps(i,:));
        case 3
            smap = squeeze(spike_maps(i,:,:));
    end
    [infoPerSpike(i), infoRate(i)] = infoSkaggs(smap, occupancy_map);
    spkRate(i)  = sum(spks(i,:))/nframes;
    peakRate(i) = max(smap(:)./occupancy_map(:));
end
sparsity = pfield_sparsity(pfields);

p.infoPerSpike = infoPerSpike;
p.infoRate = infoRate;
p.spkRate = spkRate;
p.peakRate = peakRate;
p.sparsity = sparsity;
end