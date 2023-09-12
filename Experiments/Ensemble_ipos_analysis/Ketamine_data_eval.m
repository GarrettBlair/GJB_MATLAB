% fnames = {'2022-07-27_07-41-00_M015_SAL_PFC_HPC_0_0_0mpk_HPCprobe_extracted.mat',...
%     '2022-07-28_13-19-00_M016_SAL_PFC_HPC_0_0_0mpk_HPCprobe_extracted.mat',...
%     '2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted.mat',...
%     '2022-08-03_12-45-00_M015_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted.mat',...
%     '2022-08-04_03-45-00_M016_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted.mat',...
%     '2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted.mat',...
%     '2022-08-08-04-05-00_M017_SAL_mPFC_HPC_0_0_0mpk_HPCprobe_extracted.mat',...
%     '2022-08-11-01-55-00_M018_SAL_mPFC_HPC_0_0_0mpk_HPCprobe_extracted.mat',...
%     '2022-08-12-02-35-00_M017_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted.mat',...
%     '2022-08-13-03-57-00_M018_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted.mat',...
%     '2022-08-15-02-25-00_M018_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted.mat'};
fnames = {'2022-07-27_07-41-00_M015_SAL_PFC_HPC_0_0_0mpk_HPCprobe_extracted_recency.mat',...
    '2022-07-28_13-19-00_M016_SAL_PFC_HPC_0_0_0mpk_HPCprobe_extracted_recency.mat',...
    '2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted_recency.mat',...
    '2022-08-03_12-45-00_M015_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted_recency.mat',...
    '2022-08-04_03-45-00_M016_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted_recency.mat',...
    '2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted_recency.mat',...
    '2022-08-08-04-05-00_M017_SAL_mPFC_HPC_0_0_0mpk_HPCprobe_extracted_recency.mat',...
    '2022-08-11-01-55-00_M018_SAL_mPFC_HPC_0_0_0mpk_HPCprobe_extracted_recency.mat',...
    '2022-08-12-02-35-00_M017_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted_recency.mat',...
    '2022-08-13-03-57-00_M018_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted_recency.mat',...
    '2022-08-15-02-25-00_M018_RSK_mPFC_HPC_3_10_30mpk_HPCprobe_extracted_recency.mat'};


nf = length(fnames);
distr_dist_hpc = NaN(nf, 4, 4);
kl_div_hpc = NaN(nf, 4, 4);
distr_dist_pfc = NaN(nf, 4, 4);
kl_div_pfc = NaN(nf, 4, 4);
rsk_label = false(nf,1);
for i = 1%:nf
    %%
    temp = load(fnames{i});
    h_sn = double(temp.HPC_spks_norm);
    p_sn = double(temp.PFC_spks_norm);
    h_s = double(temp.HPC_spks);
    p_s = double(temp.PFC_spks);
    hsamp = size(h_s,2);
    psamp = size(p_s,2);
    min_samples = min(hsamp, psamp);
    h_s = h_s(:,1:min_samples);
    p_s = p_s(:,1:min_samples);
    h_sn = h_sn(:,1:min_samples);
    p_sn = p_sn(:,1:min_samples);
    timing = ones(10, size(h_s,2));
    
    h_en = temp.HPC_ensemble_prob;
    p_en = temp.PFC_ensemble_prob;
%     h_en = temp.HPC_ensemble_dist;
%     p_en = temp.PFC_ensemble_dist;
    
    en_hists = cell(2,4);
    hpc_vals = cell(1,4);
    pfc_vals = cell(1,4);
    probbins = [0:.05:1];
    for j = 1:3
        inds = 1800*(j-1)+1:1800*j;
        timing(:, inds) = j;
        en_hists{1,j} = histcounts(h_en(inds), probbins, 'Normalization', 'probability');
        en_hists{2,j} = histcounts(p_en(inds), probbins, 'Normalization', 'probability');
        hpc_vals{1,j} = h_en(inds);
        pfc_vals{1,j} = p_en(inds);
    end
    j = 4; inds = 1800*(j-1)+1:min_samples;
    timing(:, inds) = j;
    en_hists{1,j} = histcounts(h_en(inds), probbins, 'Normalization', 'probability');
    en_hists{2,j} = histcounts(p_en(inds), probbins, 'Normalization', 'probability');
    hpc_vals{1,j} = h_en(inds);
    pfc_vals{1,j} = p_en(inds);
    
    for j = 1:3
        for k = j+1:4
            p1 = en_hists{1, j};
            p1(p1==0) = eps;
            p2 = en_hists{1, k};
            p2(p2==0) = eps;
            p1 = p1./sum(p1);
            p2 = p2./sum(p2);
            distr_dist_hpc(i, j, k) = (nanmedian(hpc_vals{k}) - nanmedian(hpc_vals{j}));
            kl_div_hpc(i, j, k) = kldiv(probbins(1:end-1), p1, p2, 'sym');
            kl_div_hpc(i, k, j) = kl_div_hpc(i, j, k); % kldiv(probbins(1:end-1), p2, p1);
            distr_dist_hpc(i, k, j) = distr_dist_hpc(i, j, k); % kldiv(probbins(1:end-1), p2, p1);
            
            p1 = en_hists{2, j};
            p1(p1==0) = eps;
            p2 = en_hists{2, k};
            p2(p2==0) = eps;
            p1 = p1./sum(p1);
            p2 = p2./sum(p2);
            distr_dist_pfc(i, j, k) = (nanmedian(pfc_vals{k}) - nanmedian(pfc_vals{j})); % abs(nanmedian(p1) - nanmedian(p2));
            kl_div_pfc(i, j, k) = kldiv(probbins(1:end-1), p1, p2, 'sym');
            kl_div_pfc(i, k, j) = kl_div_pfc(i, j, k); % kldiv(probbins(1:end-1), p2, p1);
            distr_dist_pfc(i, k, j) = distr_dist_pfc(i, j, k); % kldiv(probbins(1:end-1), p2, p1);
        end
    end
    
    all_s = [h_s; timing; p_s];
    all_sn = [h_sn; timing; p_sn];
    
    h_i = double(temp.HPC_isoMap');
    p_i = double(temp.PFC_isoMap');
    hsamp = size(h_i,2);
    psamp = size(p_i,2);
    h_i = normalize_rows(h_i(:,1:min(hsamp, psamp)));
    p_i = normalize_rows(p_i(:,1:min(hsamp, psamp)));

%     figure(100*i); clf
%     set(gcf, 'Position', [680   313   759   665])
%     subplot_tight(2,1,1)
%     imagesc(normalize_rows(all_s)); 
%     subplot_tight(2,1,2)
%     imagesc(normalize_rows(all_sn)); 

    all_iso = ([h_i; timing; p_i]);
    if true
    figure(i); clf
    set(gcf, 'Position', [680   313   759   665])
    subplot_tight(3,2,1)
    imagesc(normalize_rows(all_s)); 
    if contains(temp.fname_HPC(end-60:end), '_RSK_')
        colormap(magma)
        cm = magma(4);
        rsk_label(i)=1;
    else
        colormap(viridis)
        cm = viridis(4);
        rsk_label(i)=0;
    end
    title(['...' temp.fname_HPC(end-60:end)], 'interpreter', 'none')
    
    subplot_tight(3,2,3); hold on
    plot(temp.HPC_ensemble_prob, 'b')
    plot(temp.PFC_ensemble_prob*-1, 'c')
    axis([0 length(temp.PFC_ensemble_prob) -1.1 1.1])

    t = squeeze(timing(1,:));
    subplot_tight(3,2,5); hold on
    scatter3(h_i(1,:),h_i(2,:),h_i(3,:), 5, 'k', 'Marker', 'o', 'MarkerEdgeAlpha', .2, 'MarkerFaceColor', 'none')
    scatter3(h_i(1,:),h_i(2,:),h_i(3,:), 15, cm(int16(t),:), 'Marker', '.', 'MarkerEdgeAlpha', .75)
    axis([0 1 0 1 0 1])
    set(gca, 'View', [36 15], 'Box', 'on')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot_tight(3,2,6); hold on
    scatter3(p_i(1,:),p_i(2,:),p_i(3,:), 5, 'k', 'Marker', 'o', 'MarkerEdgeAlpha', .2, 'MarkerFaceColor', 'none')
    scatter3(p_i(1,:),p_i(2,:),p_i(3,:), 15, cm(int16(t),:), 'Marker', '.', 'MarkerEdgeAlpha', .75)
    axis([0 1 0 1 0 1])
    set(gca, 'View', [36 15], 'Box', 'on')
    
    subplot_tight(3,2,4); hold on
    for j = 1:4
%     subplot_tight(3,8,12+j); hold on
        plot(probbins(1:end-1), en_hists{1,j}./max(en_hists{1,j}), 'Color', cm(j,:)/1.2, 'LineWidth', 2)
        plot([nanmedian(hpc_vals{j}) nanmedian(hpc_vals{j})], [0 1.2], 'Color', cm(j,:)/1, 'LineWidth', 2)
        plot(probbins(1:end-1), en_hists{2,j}./max(en_hists{2,j}) -2, 'Color', cm(j,:)/1.2, 'LineWidth', 2)
        plot([nanmedian(pfc_vals{j}) nanmedian(pfc_vals{j})], [0 1.2]-2, 'Color', cm(j,:)/1, 'LineWidth', 2)
    end
    axis([-.1 1.1 -2.5 1.5])
    subplot_tight(3,4,3); hold on
    ch = corr([h_s]);
    cp = corr([p_s]);
    imagesc(ch, [0 1])
    axis image
    subplot_tight(3,4,4); hold on
    imagesc(cp, [0 1])
    axis image
    if contains(temp.fname_HPC(end-60:end), '_RSK_')
        colormap(magma)
        cm = magma(4);
    else
        colormap(viridis)
        cm = viridis(4);
    end

    end
end
%%
varrs = {'kl_div_hpc' 'kl_div_pfc' 'distr_dist_hpc' 'distr_dist_pfc'};
figure(1003); clf
for v = 1:length(varrs)
eval(sprintf('a = %s(rsk_label==true,:,:);', varrs{v}));
eval(sprintf('b = %s(rsk_label==false,:,:);', varrs{v}));
figure(1002); 
subplot(2,2,v)
imagesc(squeeze(nanmean(a,1)))
if any(v==[1,2])
    caxis([0 40])
else
    caxis([-.5 1])
end
figure(1002);
subplot(2,2,v)
imagesc(squeeze(nanmean(b,1)))
if any(v==[1,2])
    caxis([0 40])
else
    caxis([-.5 1])
end
figure(1003); 
subplot(2,2,v)
hold on;
for i = 1:size(b,1)
    plot([1,2,3]-.1, [b(i, 1, 2), b(i, 2, 3), b(i, 3, 4)], 'k')
    plot([1,2,3]-.1, [b(i, 1, 2), b(i, 2, 3), b(i, 3, 4)], 'k.')
%     plot([2,3], [b(i, 2, 3), b(i, 3, 4)], 'k')
%     scatter([1],b(i, 2, 3) -  b(i, 1, 2) , 'ko')
%     scatter([2],b(i, 3, 4) -  b(i, 2, 3) , 'ko')
end
for i = 1:size(a,1)
    plot([1,2,3]+.1, [a(i, 1, 2), a(i, 2, 3), a(i, 3, 4)], 'r')
    plot([1,2,3]+.1, [a(i, 1, 2), a(i, 2, 3), a(i, 3, 4)], 'r.')
%     plot([2,3], [a(i, 2, 3), a(i, 3, 4)], 'r')
%     scatter([1],a(i, 2, 3) -  a(i, 1, 2) , 'ro')
%     scatter([2],a(i, 3, 4) -  a(i, 2, 3) , 'ro')
end
title(varrs{v}, 'Interpreter', 'none')
if any(v==[1,2])
    axis([0 4 -9 40])
else
    axis([0 4 -.6 1.2])
end
mv = ceil(10*max([max(a(:)), abs(min(a(:))), max(b(:)), abs(min(b(:)))]))/10;
end



















