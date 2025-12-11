momentary_pos_info = NaN(size(spks_bin));
spks_bin = ceil(normalize_rows(spks_bin)*10);
spks_bin(isnan(spks_bin)) = 0;
max_n_spks = max(spks_bin(:));
[nsegs, nframes ] = size(spks_bin);
plotting=false;%true
if plotting==true
    figure(99); clf;
    set(gcf, 'Position', [200   359   900   341])
end
for i = 1:nsegs
    %%
    % generate the conditional probability maps for each n spikes
    if plotting==true
        figure(99); clf
    end
    sm = spks_bin(i, :);
    for spikeval = 0:max_n_spks
        smb = sm==spikeval;%(sm~=ii) = NaN;
        smb_idx = find(smb);
        p_i  = mean(smb); % overall prob
        for j = 1:length(smb_idx)
            pos_idx = find(ybin == ybin(smb_idx(j)) & xbin == xbin(smb_idx(j)));
            smb_pos = sm(pos_idx) == spikeval;
            p_ix = mean(smb_pos); % prob at position
            momentary_pos_info(i, smb_idx(j)) = p_ix * log2( p_ix / p_i  );
        end
    end
    if any(isnan(momentary_pos_info(i,:)))
        figure;
    end
    if plotting==true % for plotting and checking
        subplot(2,1,2); cla
        yyaxis('left'); cla
        plot(t_av, spks_bin(i,:), 'LineWidth', 2)
        ylim([-1 max_n_spks+1])
        hold on
        yyaxis('right'); cla
        plot(t_av, momentary_pos_info(i,:))
%         ylim([-1 1.2+max(momentary_pos_info(i,:))])
        xlim([-5 t_av(end)+5])
        drawnow
%         input('')
    end
end