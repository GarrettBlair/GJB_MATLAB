function [momentary_pos_info] = ipos_test(ms, integration_time, params)
% integration_time = .1; % time in seconds to bin activity
% spks_all = normalize_rows(ms.neuron.S_matw);
spks_bin = bin_spks_time(ms.neuron.S_matw>0, integration_time, ms.timestamps./1000, false);
max_n_spks = ceil(size(ms.neuron.S_matw,2)/size(spks_bin,2));

% [~, speed_epochs] = get_speed_epochs(ms.arena.speed_smooth, params);


params.pos_bins = [-45, -36:4:36, 45]; % in cm, x and y

bins = params.pos_bins; 
ksize = params.pfield_kernel_radius;
occupancy_thresh = params.occupancy_thresh;

xo = ms.room.x; 
yo = ms.room.y;
dt = ms.dt;
to = ms.timestamps;

x_av = average_spks_time(xo', integration_time, ms.timestamps./1000, false, 'mean');
y_av = average_spks_time(yo', integration_time, ms.timestamps./1000, false, 'mean');
t_av = average_spks_time(to', integration_time, ms.timestamps./1000, false, 'mean');


nbins = length(bins)-1;
[vmap, ~, xbin, ybin]   = make_occupancymap_2D(x_av, y_av, dt, bins, bins);
nsegs = size(spks_bin, 1);

vmap(vmap<occupancy_thresh) = NaN;
%%
momentary_pos_info = NaN(size(spks_bin));
plotting = true;
if plotting
figure(99); clf; 
set(gcf, 'Position', [200   359   900   341])
end
for i = 1:nsegs
    %%
    % generate the conditional probability maps for each n spikes
    smap_prob = NaN(max_n_spks+1, nbins, nbins);
    smb_prob = NaN(max_n_spks+1, 1);
    for spikeval = 0:max_n_spks
        sm = spks_bin(i, :);
        smb = sm==spikeval;%(sm~=ii) = NaN;
        smb_prob(spikeval+1) = sum(smb)/length(smb);
        [smap] = make_probabilitymap_2D(x_av,  y_av,  smb, bins, bins);
        smap_prob(spikeval+1,:,:) = smap;
        if plotting
        figure(99);
        subplot(2,max_n_spks+1,spikeval+1); cla
        imagesc(smap)
        end
    end
    for j = 1:length(xbin)
        nspks_ind = spks_bin(i,j) + 1;
        p_ix = smap_prob(nspks_ind, ybin(j), xbin(j));
        momentary_pos_info(i, j) = p_ix * log2( p_ix/smb_prob(nspks_ind)  );
    end
    if plotting % for plotting and checking
        subplot(2,1,2); cla
        yyaxis('left'); cla
        plot(t_av./1000, spks_bin(i,:))
        ylim([-1 max_n_spks+1])
        hold on
        yyaxis('right'); cla
        plot(t_av./1000, momentary_pos_info(i,:), 'LineWidth', 2)
        ylim([-2 11])
        xlim([-50 t_av(end)./1000+50])
        drawnow
    end
end






