function [struct_in] = construct_place_maps_2D(struct_in, x, y, dt, spks, bins, params)
% input data should already be speed thresholded
% pos_bins = params.pos_bins;
ksize = params.pfield_kernel_radius;
occupancy_thresh = params.occupancy_thresh;


nbins = length(bins)-1;
[vmap, ~]   = make_occupancymap_2D(x, y, dt, bins, bins);
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

vmap(vmap<occupancy_thresh) = NaN;

for i = 1:nsegs
    %%
    % generate the spk map
    [smap, ~] = make_occupancymap_2D(x,  y,  spks(i, :), bins, bins);
%     [binmap, ~, xbin, ybin] = make_occupancymap_2D(x,  y,  spks(i, :)>0, bins, bins);
%     spkprob = sum(spks(i,:)>0)./length(spks);
%     probmap = binmap./sum(spks(i,:)>0);
%     conditional_spkprob = zeros(length(x), 1);
%     for j = 1:length(x)
%         p_ix = probmap(ybin(j), xbin(j));
%         conditional_spkprob(j) = p_ix * log2( p_ix/spkprob  );
%     end
    % threshold based on occupancy
%     smap(vmap<occupancy_thresh) = NaN;
    spkmap(i,:,:) = smap;
    % divide by time spent for weighted place field
    p = smap./vmap;
    pfields(i,:,:) = p;%smap./vmap;
    % smoothing place map
%     smap(isnan(smap)) = 0;
%     sm_smap = conv2(smap, smoothing_kern, 'same');
    p(isnan(p)) = 0;
    psmooth = conv2(p, smoothing_kern, 'same');
    % remove the unvisited areas
    psmooth(isnan(vmap)) = NaN;
    pfields_smooth(i,:,:) = psmooth;%sm_smap./vmap;
%     pfields_smooth(i,:,:) = conv2(smap./vmap, smoothing_kern, 'same');
    if false % for plotting and checking
        figure(99); clf
        subplot(2,2,1)
        imagesc(vmap);
        title('vmap')
        subplot(2,2,2)
        imagesc(smap);
        title('spkmap')
        subplot(2,2,3)
%         imagesc(sm_smap);
        imagesc(p);
        title('smooth spkmap')
        subplot(2,2,4)
        imagesc(psmooth);
        title('smooth pfield')
        drawnow
    end
end
struct_in.vmap = vmap;
struct_in.spkmap = spkmap;
struct_in.pfields = pfields;
struct_in.pfields_smooth = pfields_smooth;
end