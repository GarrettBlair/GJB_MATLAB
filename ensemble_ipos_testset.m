
% testset = load("D:\APA recordings\Acc20832\processed_files\2023_01_13_H18_48_32_TR6_@placecells.mat");
testset = load("D:\APA recordings\Acc20832\processed_files\2023_01_28_H11_01_51_NEW6_@placecells.mat");
% spks_bin = bin_spks_time(ms.neuron.S_matw>0, integration_time, ms.timestamps./1000, false);
%%
% [ms] = downsample_ms_struct(testset.ms, [1 20000]);
ms = testset.ms;
[nsegs, nsamples] = size(ms.neuron.S_matw);
x = cos(linspace(1, 48*pi, nsamples))*40;
y = sin(linspace(1, 48*pi, nsamples))*40;
ms.room.x                = x';
ms.room.y                = y';
% diff_prop = [0:.05:1];
diff_prop = [ .05 .1:.1:.9 .95];
% noise_diff = [-.09, -.05, 0:.1:.8];
% noise_diff = [-.09, -.05, 0:.05:.8, .85, .89];
spk_thresh = single(.9);
noise_diff = single([.99 .975:-.025:.8]);% .9 .85 .8];
% diff_prop = .7;
% noise_diff = .2;
rng(1)
numCells = 300;
int_time = .25; % seconds
[qqqq, qqq, qq] = Fenton_ipos(ms, int_time, 'room', []);

ensemble_out = NaN(length(diff_prop), length(noise_diff), length(qq));
ipos_out = NaN(length(diff_prop), length(noise_diff), length(qq));
spks_av = NaN(length(diff_prop), length(noise_diff), length(qq));
split = NaN(length(diff_prop), length(noise_diff));
rand_spks = rand(numCells, nsamples);
for i = 8:length(diff_prop)
    for j = 1:length(noise_diff)
        ms.neuron.S_matw   = rand_spks;           
        cut1 = floor(nsamples*diff_prop(i))+1;
        if cut1 > nsamples
            cut1 = nsamples;
        end
        if cut1 == nsamples
            ms.neuron.S_matw(:, 1:cut1)     = ms.neuron.S_matw(:, 1:cut1)>spk_thresh;
        else
            ms.neuron.S_matw(:, 1:cut1)     = ms.neuron.S_matw(:, 1:cut1)>spk_thresh;
            ms.neuron.S_matw(:, cut1+1:end) = ms.neuron.S_matw(:, cut1+1:end)> noise_diff(j);%(spk_thresh-noise_diff(j));
        end
%         [a_theta, a_rho] = cart2pol(sub_room.x, sub_room.y);
%         [a_dist] = abs(angular_distance(a_theta, pi/2));
%         ms.neuron.S_matw(:, cut1+1:end) = ms.neuron.S_matw(:, cut1+1:end)> (spk_thresh-noise_diff(j));
        
        [ipos, subroom2, ensemble_out(i,j,:)] = Fenton_ipos(ms, int_time, 'room', []);
        spks_av(i,j,:) = nanmean(subroom2.spks_bin, 1);
        ipos_out(i,j,:) = nanmean(ipos,1);
        split(i,j) = subroom2.spks_bin_group(cut1);
    end
    i
end

ensemble_scale = ensemble_out./max(max(max(ensemble_out(:))));
ipos_scale = ipos_out./max(max(max(ipos_out(:))));
spk_scale = spks_av./max(max(max(spks_av(:))));
% ensemble_out = ensemble_out./max(ensemble_out(:));
%%
e_pwr = NaN(length(diff_prop), length(noise_diff));
i_pwr = NaN(length(diff_prop), length(noise_diff));
clrs = jet(length(noise_diff));
% figure(54); clf
figure(53); clf
rng(1)
con_j = find(single(noise_diff)==single(spk_thresh));
for i = 1:length(diff_prop)
    for j = 1:length(noise_diff)
        cut1 = floor(nsamples*diff_prop(i));
        sp = squeeze(squeeze(spk_scale(i,j,:)));
        ind = sub2ind([length(noise_diff) length(diff_prop) ], j, i);
%         figure(54); subplot_tight(length(diff_prop),length(noise_diff), ind); hold on
%         axis off tight
%         ylim([-.5 1.5])
        figure(53); subplot_tight(length(diff_prop),length(noise_diff), ind, [0 0]);
        e  = (squeeze(squeeze(ensemble_scale(i,j,:))));
        ip = (squeeze(squeeze(ipos_scale(i,j,:))));
%         s2 = interp1(1/length(s):1/length(s):1, s, 1/length(e):1/length(e):1);
        plot(sp-1, 'Color', clrs(j,:)/3); hold on
        plot(e, 'Color', clrs(j,:)/2, 'LineWidth', 2);
        plot(ip+1, '-','Color', clrs(j,:)/2, 'LineWidth', 2);
%         figure(53+ind); clf
%         plot(s2-1, 'Color', clrs(j,:)/3); hold on
%         plot(e, 'Color', clrs(j,:)/2, 'LineWidth', 2);
%         plot(ip+1, '.-', 'Color', clrs(j,:)/2, 'LineWidth', 2);
        ylim([-1 3])
        axis off
        drawnow
        cut1 = split(i,j);
        e  = (squeeze(squeeze(ensemble_out(i,j,:))));
        e_control = (squeeze(squeeze(ensemble_out(i,con_j,:))));
        ip = (squeeze(squeeze(ipos_out(i,j,:))));
        ip_control = (squeeze(squeeze(ipos_out(i,con_j,:))));
%         e_pwr(i,j) = nansum(e(cut1+1:end))/nansum(e);
%         i_pwr(i,j) = nansum(ip(cut1+1:end))/nansum(ip);
%         e_pwr(i,j) = nansum(e(cut1+1:end))/nansum(e_control(cut1+1:end));
%         i_pwr(i,j) = nansum(ip(cut1+1:end))/nansum(ip_control(cut1+1:end));
        e_pwr(i,j) = nansum(e-e_control);
        i_pwr(i,j) = nansum(ip-ip_control);
    end
end

figure(55); clf
subplot(1,2,1);
imagesc(e_pwr)%, [0 1])
title('Ensemble info')
set(gca, 'XTickLabel', noise_diff-spk_thresh, 'XTick', 1:length(noise_diff), 'XTickLabelRotation', -90,...
    'YTickLabel', 1-diff_prop, 'YTick', 1:length(diff_prop))
xlabel('Firing threshold difference')
ylabel('Proportion of time modulated')
subplot(1,2,2); 
imagesc(i_pwr)%, [0 1])
title('Positional info I_{pos}')
set(gca, 'XTickLabel', noise_diff-spk_thresh, 'XTick', 1:length(noise_diff), 'XTickLabelRotation', -90,...
    'YTickLabel', diff_prop, 'YTick', 1:length(diff_prop))
xlabel('Firing threshold difference')
ylabel('Proportion of time modulated')
colormap viridis















