function [struct_in] = construct_place_maps_2Dpolar(struct_in, x, y, dt, spks, xbins, ybins, params)
% input data should already be speed thresholded
% pos_bins = params.pos_bins;
ksize = params.pfield_kernel_radius;
occupancy_thresh = params.occupancy_thresh;

n_xbins = length(xbins)-1;
n_ybins = length(ybins)-1;
[vmap, vmap_counts, xbin, ybin]    = make_occupancymap_2D(x, y, dt, xbins, ybins);

% firsthalf = cumsum(dt)<=sum(dt)/2;
nsplits = 2;

var_binned = sub2ind([size(vmap_counts)], ybin, xbin);
counts = vmap_counts(:)';
if isfield(struct_in, 'splits')
    splits = struct_in.splits;
else
    splits = splits_equal_occupancy(var_binned, counts, nsplits);
end
% splits = uint8(mod(cumsum(dt)./ceil(sum(dt)/(nsplits/2)), 1)<.5);
% splits(splits==0) = 2;
% splits = uint8(splits); 

[vmap1, ~]   = make_occupancymap_2D(x(splits==1), y(splits==1), dt(splits==1), xbins, ybins);
[vmap2, ~]   = make_occupancymap_2D(x(splits==2), y(splits==2), dt(splits==2), xbins, ybins);
% vmap(vmap<occupancy_thresh) = NaN;
nsegs = size(spks, 1);
spkmap = NaN(nsegs, n_ybins, n_xbins);
spkmap_smooth = NaN(nsegs, n_ybins, n_xbins);
pfields = NaN(nsegs, n_ybins, n_xbins);
pfields_smooth = NaN(nsegs, n_ybins, n_xbins);
pfields_split1 = NaN(nsegs, n_ybins, n_xbins);
pfields_smooth_split1 = NaN(nsegs, n_ybins, n_xbins);
pfields_split2 = NaN(nsegs, n_ybins, n_xbins);
pfields_smooth_split2 = NaN(nsegs, n_ybins, n_xbins);
split_corr = NaN(nsegs, 1);
split_p = NaN(nsegs, 1);

% smoothing_kern = ones(ksize,ksize); smoothing_kern = smoothing_kern./(sum(smoothing_kern(:)));
smoothing_kern = gausswin(ksize*2+1);
% smoothing_kern = gausswin(ksize+1);
smoothing_kern = smoothing_kern*smoothing_kern';
smoothing_kern = smoothing_kern./(sum(smoothing_kern(:)));

vmap(vmap<occupancy_thresh) = NaN;
vmap1(vmap1<occupancy_thresh/2) = NaN;
vmap2(vmap2<occupancy_thresh/2) = NaN;
plotting = false;
for i = 1:nsegs
    %%
    % generate the spk map
    [smap, ~] = make_occupancymap_2D(x,  y,  spks(i, :), xbins, ybins);
    % threshold based on occupancy
    spkmap(i,:,:) = smap;
    % divide by time spent for weighted place field
    p = smap./vmap;
    pfields(i,:,:) = p;%smap./vmap;
    % smoothing place map
    bad_p = isnan(p);
    p(isnan(p)) = 0;
    psmooth = conv2(p, smoothing_kern, 'same');
    % remove the unvisited areas
    psmooth(bad_p) = NaN;
    pfields_smooth(i,:,:) = psmooth;%sm_smap./vmap;
    spkmap_smooth(i,:,:)  = psmooth.*vmap;
    
    
    % doing split half correlation metric
    [smap1, ~] = make_occupancymap_2D(x(splits==1),  y(splits==1),  spks(i, (splits==1)), xbins, ybins);
    [smap2, ~] = make_occupancymap_2D(x(splits==2),  y(splits==2),  spks(i, (splits==2)), xbins, ybins);
    p1 = smap1./vmap1;
    pfields_split1(i,:,:) = p1;
    p2 = smap2./vmap2;
    pfields_split2(i,:,:) = p2;
    
    % smoothing
    bad_p = isnan(p1);
    p1(bad_p) = 0;
    p3 = cat(2, p1, p1, p1);
    p1smooth = conv2(p3, smoothing_kern, 'same');
    p1smooth = p1smooth(:, length(p1)+1:2*length(p1));
%     p1smooth = conv2(p1, smoothing_kern, 'same');
    p1smooth(bad_p) = NaN;
    bad_p = isnan(p2);
    p2(bad_p) = 0;
    p3 = cat(2, p2, p2, p2);
    p2smooth = conv2(p3, smoothing_kern, 'same');
    p2smooth = p2smooth(:, length(p2)+1:2*length(p2));
%     p2smooth = conv2(p2, smoothing_kern, 'same');
    p2smooth(bad_p) = NaN;
    
    bads = isnan(p1smooth.*p2smooth);% | (p1+p2==0);
    if any(~bads(:))
        [split_corr(i), split_p(i)] = corr(p1smooth(~bads), p2smooth(~bads));
    end
    pfields_smooth_split1(i,:,:) = p1smooth;
    pfields_smooth_split2(i,:,:) = p2smooth;
    %     pfields_smooth(i,:,:) = conv2(smap./vmap, smoothing_kern, 'same');
    if plotting == true % for plotting and checking
        %%
        ps = squeeze(pfields_smooth(i,:,:));
        if any(~bads(:))
            max_val = max([p1smooth(:); p2smooth(:); ps(:)]);
        else
            max_val = 1;
        end
        figure(99); clf
        subplot(2,3,1)
        imagesc(vmap);
        title('vmap')
        subplot(2,3,2)
        imagesc(smap); hold on
        title('spkmap')
        subplot(2,3,3)
%         imagesc(sm_smap);
        imagesc(squeeze(pfields(i,:,:)));
        title('raw pfield')
        subplot(2,3,6)
        imagesc(squeeze(pfields_smooth(i,:,:)), [-max_val max_val]);
        title('smooth pfield')
        subplot(2,3,4)
        imagesc(p1smooth, [-max_val max_val]);
        title(sprintf('1st split c=%.2f', split_corr(i)))
        subplot(2,3,5)
        imagesc(p2smooth, [-max_val max_val]);
        title(sprintf('2nd split p=%0.4f', split_p(i)))
        drawnow
%         pause(.5)
    end
end
struct_in.vmap           = vmap;
struct_in.vmap_counts    = vmap_counts;
struct_in.vmap_split1    = vmap1;
struct_in.vmap_split2    = vmap2;
struct_in.spkmap         = spkmap;
struct_in.spkmap_smooth  = spkmap_smooth;
struct_in.pfields        = pfields;
struct_in.pfields_smooth = pfields_smooth;
struct_in.pfield_split_vec      = splits;
struct_in.pfields_split1        = pfields_split1;
struct_in.pfields_smooth_split1 = pfields_smooth_split1;
struct_in.pfields_split2        = pfields_split2;
struct_in.pfields_smooth_split2 = pfields_smooth_split2;
struct_in.split_corr     = split_corr;
struct_in.split_p        = split_p;
end