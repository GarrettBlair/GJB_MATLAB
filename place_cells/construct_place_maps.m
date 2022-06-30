function [struct_in] = construct_place_maps(struct_in, x, y, dt, spks, params)
% input data should already be speed thresholded
pos_bins = params.pos_bins;
ksize = params.pfield_kernel_radius;
occupancy_thresh = params.occupancy_thresh;


nbins = length(pos_bins)-1;
[vmap, ~]   = make_occupancymap_2D(x, y, dt, pos_bins, pos_bins);
% vmap(vmap<occupancy_thresh) = NaN;
nsegs = size(spks, 1);
spkmap = NaN(nsegs, nbins, nbins);
pfields = NaN(nsegs, nbins, nbins);
pfields_smooth = NaN(nsegs, nbins, nbins);

% smoothing_kern = ones(ksize,ksize); smoothing_kern = smoothing_kern./(sum(smoothing_kern(:)));
smoothing_kern = gausswin(ksize*2+1);
% smoothing_kern = gausswin(ksize+1);
smoothing_kern = smoothing_kern*smoothing_kern';
smoothing_kern = smoothing_kern./(sum(smoothing_kern(:)));

for i = 1:nsegs
    %%
    % generate the spk map
    [smap, ~] = make_occupancymap_2D(x,  y,  spks(i, :), pos_bins, pos_bins);
    % threshold based on occupancy
    smap(vmap<occupancy_thresh) = NaN;
    spkmap(i,:,:) = smap;
    % divide by time spent for weighted place field
    pfields(i,:,:) = smap./vmap;
    % smoothing place map
    smap(isnan(smap)) = 0;
    sm_smap = conv2(smap, smoothing_kern, 'same');
    % remove the unvisited areas
    sm_smap(isnan(vmap)) = NaN;
    pfields_smooth(i,:,:) = sm_smap./vmap;
    if false % for plotting and checking
        figure(99); clf
        subplot(2,2,1)
        imagesc(vmap);
        subplot(2,2,2)
        imagesc(smap);
        subplot(2,2,3)
        imagesc(sm_smap);
        subplot(2,2,4)
        imagesc(sm_smap./vmap);
        drawnow
    end
end
struct_in.vmap = vmap;
struct_in.spkmap = spkmap;
struct_in.pfields = pfields;
struct_in.pfields_smooth = pfields_smooth;
end