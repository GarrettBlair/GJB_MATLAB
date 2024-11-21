cd('E:\RecordingData\GarrettBlair\PKCZ_imaging\figs\andre grant oct 2024')
fs      = ["E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\processed_files\2023_07_25_H17_11_23_TR7_@placecells_HPC_miniscope1.mat",...
    "E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\processed_files\2023_08_12_H17_39_31_RET8_@placecells_HPC_miniscope1.mat"];
fs(2,:) = ["E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\processed_files\2023_07_25_H17_44_28_TR7_@placecells_HPC_miniscope1.mat",...
    "E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\processed_files\2023_08_12_H16_55_19_RET8_@placecells_HPC_miniscope1.mat"];

fs(3,:) = ["E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC24457\processed_files\2024_06_06_H16_48_40_TR19_@placecells_HPC_miniscope1.mat",...
    "E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC24457\processed_files\2024_06_27_H16_15_20_RET20_@placecells_HPC_miniscope1.mat"];
fs(4,:) = ["E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC24458\processed_files\2024_05_30_H13_29_16_TR12_@placecells_HPC_miniscope1.mat",...
    "E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC24458\processed_files\2024_06_19_H13_05_43_TR13_@placecells_HPC_miniscope1.mat"];
fs(5,:) = ["E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC24459\processed_files\2024_05_30_H14_20_04_TR13_@placecells_HPC_miniscope1.mat",...
    "E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC24459\processed_files\2024_06_19_H13_39_28_RET14_@placecells_HPC_miniscope1.mat"];
conts = ["E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\matching_contours\manual_alignment_HPC\2023_07_25_H17_11_23_TR7_@placecells_HPC_miniscope1.mat",...
    "E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23459\matching_contours\manual_alignment_HPC\2023_08_12_H17_39_31_RET8_@placecells_HPC_miniscope1.mat"];
conts(2,:) = ["E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\matching_contours\manual_alignment_HPC\2023_07_25_H17_44_28_TR7_@placecells_HPC_miniscope1.mat",...
    "E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\matching_contours\manual_alignment_HPC\2023_08_12_H16_55_19_RET8_@placecells_HPC_miniscope1.mat"];
numan = 5;
data = cell(numan,2);
contours = cell(numan,2);
for an = 1:numan
figure(10+an); clf;
subplot_tight(2,2,1);
title('Pre'); axis off
subplot_tight(2,2,2);
title('Post'); axis off
for i = 1:2
    %%
    temp = load(fs{an,i});
    data{an,i} = temp;
    if an < 3
        temp2 = load(conts{an,i});
        im = temp.ms.neuron.maxFrame;
        im = im ./ max(im(:));
        c = temp2.contours_shifted;
        c = normalize_matrix(c);
        c(c<=.6) = 0;
        c = normalize_matrix(c);
        contours{an,i} = c;
        
        [ns, ~, ~] = size(c);
        cm = jet(ns);
        cm = cm(randperm(ns), :);
        [rgb_im] = color_contours_im(c, cm);%, im, 2);
        subplot_tight(2,2,i);
        imagesc(im, [.1 1.1]); colormap gray; axis image off
        subplot_tight(2,2,i+2);
        image(rgb_im*2); axis image off
        title(sprintf('# cells = %d', ns));
    end
end
end
%%
inds = 1:5000;
close all
% cells = [1:10:300];
for an = 3 % 1:numan
for i = 1:2
    temp = data{an,i};
    spks = temp.ms.spks;
    craw = temp.ms.neuron.YrA + temp.ms.neuron.C;
    figure(an)
    imagesc(normalize_rows(spks>0))
    
    craw = normalize_rows(craw);
    spks2 = craw;
    spks2(spks==0) = NaN;
    figure(10+an); clf
    stacked_traces(craw(:, inds), 1, {'k'})
    stacked_traces(spks2(:, inds), 1, {'r.'})
end
end
%%
si_bins = linspace(0, 1.5, 50);
si_bins_cent = si_bins(1:end-1) + mean(abs(diff(si_bins)));
kern = gausswin(11); kern = kern*kern';
kern = kern./(sum(kern(:)));
percent_place = NaN(numan,2);
ind = 0;
%
for an = 1:5%1:numan
histfig = figure(102+an*10); clf
for i = 1:2
    %%
    temp = data{an,i};
    if an < 3
        temp2 = contours{an,i};
    else
        temp2=[];
    end
    %%
    r = temp.ms.room;
    c = r.split_corr;
    si = r.pcell_stats.infoPerSpike;
%     si = r.pcell_stats.infoProb;
    percent_place(an,i) = mean(r.pcell_stats.infoProb<=.05);
%     h = histcounts(si, si_bins, 'Normalization', 'probability');
%     figure(histfig); hold on; plot(si_bins_cent, h)
    good = r.pcell_stats.infoProb <= .05 & r.split_p <=.05;
%     good = r.pcell_stats.infoProb <= .01 & r.pcell_stats.sparsity >.5;
    nspks = sum(temp.ms.spks,2);
    nspks = r.pcell_stats.infoPerSpike.*nspks;
    nspks(good==0) = 0;
    [~, ord] = sort(nspks, 'descend');
%     si(good==0) = 0;
%     [~, ord] = sort(si, 'descend');
    %
    figure(i+1000+an*10); clf
    set(gcf, 'Position', [50+ind*50   359   405   557])
    xbins = -24:2:24;
    ybins = -24:2:24;
    for j = 1:4
        cellnum = ord(j*2);
       subplot(4,3,(j-1)*3 + 1); hold on
       if ~isempty(temp2)
       cell_inset(temp2, cellnum);
       end
       subplot(4,3,(j-1)*3 + 2); hold on
       plot(r.x, r.y, 'k-')
       s = temp.ms.spks(cellnum,:);% > 0;
       scatter(r.x(s>0), r.y(s>0), s(s>0)*200, 'r.')
       axis image off
       [vmap, vmap_counts, xbin, ybin]    = make_occupancymap_2D(r.x, r.y, r.x*0 +1, xbins, ybins);
       [smap, ~, ~, ~]    = make_occupancymap_2D(r.x, r.y, s, xbins, ybins);
       vmap(vmap<5) = NaN;
       p = smap./vmap;
%         vmap = r.vmap_counts;
%         p = squeeze(r.pfields(cellnum,:,:));
       alpha = vmap;
       
       p(isnan(vmap)) = 0;
       p2 = conv2(p, kern, 'same');
       p2(isnan(vmap)) = 0;
       
       
       subplot(4,3,(j-1)*3 + 3); hold on
       imagesc(p2, 'AlphaData', ~isnan(alpha)); axis image off
       title(sprintf('%2.2f bits/s', si(cellnum)))
       %        p = squeeze(r.pfields_smooth(ord(j), :, :));
%        imagesc(p);
    end
    colormap magma
    drawnow()
    ind = ind+1;
end
end
%%
figure(1023);
clf; hold on
bar(mean(100*percent_place,1))
plot(100.*percent_place', 'k-')
scatter(ones(numan,1)*1, 100*percent_place(:,1), 'r')
scatter(ones(numan,1)*2, 100*percent_place(:,2), 'r')
axis([ .5 2.5 -2.5 52.5])
set(gca, 'XTick', [1 2], 'XTickLabel', {'PRE' 'POST'}, 'YTick', 100*[0:.1:.5])
ylabel('% place cells')
[h,p,stats] = ttest2(percent_place(:,1), percent_place(:,2))
%%
function cell_inset(temp2, cellnum)
[ns, h, w] = size(temp2);
win = 20;
c = squeeze(temp2(cellnum,:,:));
w1 = find(any(c,1), 1, 'first') - win;
w2 = find(any(c,1), 1, 'last')  + win;
h1 = find(any(c,2), 1, 'first') - win;
h2 = find(any(c,2), 1, 'last')  + win;

% c2 = c(h1:h2, w1:w2);

w1 = max(1, w1);
w2 = min(w, w2);
h1 = max(1, h1);
h2 = min(h, h2);

cm = .3.*ones(ns,3);
cm(cellnum, :) = [.3, 1, 1];
[rgb_im] = color_contours_im(temp2, cm);%, im, 2);
image(rgb_im)
axis([w1 w2 h1 h2]); axis square off
end