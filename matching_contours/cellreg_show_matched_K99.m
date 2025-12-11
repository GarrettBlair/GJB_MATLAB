clear 

% ddir = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_HPC\cellreg_subIL\";
% crfile = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_HPC\cellreg_subIL\cellRegistered_20251014_165122.mat";
% ridx = 1
ddir = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_ACC\cellreg_subIL\";
crfile = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_ACC\cellreg_subIL\cellRegistered_20251014_165409.mat";
ridx = 2


% ddir = "D:\GarrettBlair\APA\HPCACC24500\matching_contours\manual_alignment_ACC\cellreg\";
% crfile = "D:\GarrettBlair\APA\HPCACC24500\matching_contours\manual_alignment_ACC\cellreg\cellRegistered_20231116_151132.mat";

% ddir = "D:\GarrettBlair\APA\HPCACC24500\matching_contours\manual_alignment_ACC\cellreg\";
% crfile = "D:\GarrettBlair\APA\HPCACC24500\matching_contours\manual_alignment_HPC\cellreg\cellRegistered_20231116_152054.mat";

temp = load(crfile);
cmap = temp.cell_registered_struct.cell_to_index_map;
ns = size(cmap, 2);
%%

figure; 
subplot(1,3,1)
[o, l] = pop_overlap(cmap);
imagesc(cmap)
subplot(1,3,2)
imagesc(o, [0, 1])
colorbar
subplot(1,3,3)
shadedErrorBar(1:ns, nanmean(l,2), nanstd(l,[],2))
ylim([0 1])
colormap magma
%%
[~,ny,nx] = size(temp.cell_registered_struct.spatial_footprints_corrected{1});
close all
contour_thresh = .65;
for i = 2% 1:ns
%     shared = cmap(:,j)>0 & cmap(:,i)>0;
    c_im = zeros(ny,nx,3);
    idx=0;
    for j = [i-1:i+1]
        idx=idx+1;
        if j>0 && j<=ns
            if j == i
                if j==1
                    shared = (cmap(:,i+1)>0 & cmap(:,i)>0);
                elseif j==ns
                    shared = (cmap(:,i-1)>0 & cmap(:,i)>0);
                else
                    shared = (cmap(:,i+1)>0 & cmap(:,i)>0) | (cmap(:,i-1)>0 & cmap(:,i)>0);
                end
            else
                shared = cmap(:,j)>0 & cmap(:,i)>0;
            end
            segs = cmap(shared,j);
            cs = temp.cell_registered_struct.spatial_footprints_corrected{j}(segs,:,:);
            for k = 1:size(cs,1); 
                ccs = squeeze(cs(k,:,:));
                ccs = ccs-min(ccs(:));
                ccs = ccs./max(ccs(:));
                ccs(ccs<=contour_thresh) = 0;
                cs(k,:,:) = ccs;
            end
                
%             cs(cs<=.4)=0;
            c_im(:,:, mod(j,3)+1) = sum(cs,1).*1;
        end
    end    
    figure(i); 
    image(c_im)
    title(i)
end
%%
sessns = 1:3;
contour_thresh = .5;
ns=length(sessns);

sharedvals = [.4 .4 .4];
othervals  = [.4 .4 .4];
figure(1001); clf
% set(gcf, 'Position', [200, 200, 1000, 450])
set(gcf, 'Position', [367   170   677   662])

for ridx = 1:2
    if ridx ==1
        ddir = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_HPC\cellreg_subIL\";
        crfile = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_HPC\cellreg_subIL\cellRegistered_20251014_165122.mat";
%         ridx = 1
    else
        ddir = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_ACC\cellreg_subIL\";
        crfile = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_ACC\cellreg_subIL\cellRegistered_20251014_165409.mat";
%         ridx = 2
        
    end
    
    temp = load(crfile);
    cmap = temp.cell_registered_struct.cell_to_index_map;
%     ns = size(cmap, 2);
    
[~,ny,nx] = size(temp.cell_registered_struct.spatial_footprints_corrected{1});
shared = sum(cmap(:,sessns)>0, 2) == ns; %  & cmap(:,i)>0
    
    % figure(1001); clf
    % set(gcf, 'Position', [200, 200, 1000, 450])
    set(gcf, 'Position', [367   170   677   662])
    
    for i = 1:ns
        %     shared = cmap(:,j)>0 & cmap(:,i)>0;
        c_im = zeros(ny,nx,3);
        segs1 = cmap(shared,i);
        segs2 = cmap(shared==false & cmap(:,i)>0, i);
        im1 = seg_im(temp.cell_registered_struct.spatial_footprints_corrected{i}(segs1,:,:), contour_thresh);
        
        im2 = seg_im(temp.cell_registered_struct.spatial_footprints_corrected{i}(segs2,:,:), contour_thresh);
        
        c_im(:,:, i) = im1*.95;
        for ii = 1:3; c_im(:,:, ii) = squeeze(c_im(:,:, ii)) + squeeze(im2*.4); end
        for ii = 1:3
            if i~=ii
                c_im(:,:, ii) = squeeze(c_im(:,:, ii)) + squeeze(im1*.2);
            end
        end
        subplot_tight(2,ns,ns*(ridx-1)+i, [.1,.01])
        image(c_im)
        axis image off
        title(sprintf('sess: %d,  %d of %d', sessns(i), sum(shared), sum(cmap(:,i)>0)))
        
    end
    disp(sum(shared))
end
%%
% close all
use_simple = false;

matchdirs = [   "D:\GarrettBlair\APA\HPCACC34990\processed_files\subIL_cellmatching_HPC.mat",
                "D:\GarrettBlair\APA\HPCACC34990\processed_files\subIL_cellmatching_ACC.mat"];
% matchdirs = [   "D:\GarrettBlair\APA\HPCACC34990\processed_files\RAR_cellmatching_ACC.mat",
%                 "D:\GarrettBlair\APA\HPCACC34990\processed_files\RAR_cellmatching_HPC.mat"];
if use_simple == false
    datadir = "D:\GarrettBlair\APA\HPCACC34990\processed_files\";
else
    datadir = "D:\GarrettBlair\APA\HPCACC34990\simple_files\";
end
manstr = 'manual_alignment_';
sessns = [1,2,3];
matdat_all = cell(2,3);
ms_images = cell(2,3);
cmap_all = cell(2,1);
nsess = length(sessns);
data_fn = cell(2, nsess);
for ii = 1:2
    temp = load(matchdirs(ii));
    cmap_all{ii} = temp.cellmap(:,sessns);
    for jj = 1:nsess
        fn = temp.file_names{sessns(jj)};
        datastring = fn(strfind(fn, manstr)+length(manstr)+4:length(fn));
        data_fn{ii,jj} = sprintf('%s%s', datadir, datastring);
        disp(data_fn{ii,jj})
        matdat_all{ii,jj} = load(data_fn{ii,jj});
        if use_simple == false
            ms_images{ii,jj} = matdat_all{ii,jj}.ms.neuron.meanFrame;
        end
    end
end
for ii = 1:2
    for jj = 1:nsess
        if use_simple == false
            ms_images{ii,jj} = matdat_all{ii,jj}.ms.neuron.maxFrame;
        end
    end
end
figure; set(gcf, 'Position', [367   170   677   662])
subplot_tight(2,1,1,[.1, .01])
imagesc(ms_images{2,1}, [50 200]); colormap gray; axis image off
subplot_tight(2,1,2,[.1, .01])
imagesc(ms_images{1,1}, [50 200]); colormap gray; axis image off


% prob_corr = nan(2,nsess);
% prob_left = nan(2,nsess);
% for ii = 1:2
%     for jj = 1:3 % 1:nsess
%         matdat = matdat_all{ii,jj}; % load(data_fn);
%         cmap_sub = cmap(shared, jj);
% %         prob_corr(ii,jj)= matdat.ms.zonestruct.metrics.p_correct;
% %         prob_left(ii,jj)= matdat.ms.zonestruct.metrics.p_left;
%         prob_corr(ii,jj)= matdat.p_correct;
%         prob_left(ii,jj)= matdat.p_correct;
%     end
% end

t_max = nan(2,3);
for ii = 1:2
    for jj = 1:nsess
        matdat = matdat_all{ii,jj}; % load(data_fn);
        if use_simple == false
            t    = matdat.ms.timestamps;
            t_dt    = matdat.ms.dt';
        else
            t    = matdat.time_ms;
            t_dt    = matdat.dt';
        end
        t_max(ii,jj) = t(end)-t(1);
    end
end
for jj = 1:nsess
    t_max(:,jj) = min(t_max(:,jj));
end
%%
% close all
dt=100;
movmean_k = 5;
offset=0;
bins = linspace(-1,1,41);
n = floor((length(bins)-1)/2);
pcf = figure(10001); clf
lc = nan(2,n,3);
le = nan(2,n,3);
fieldsdist = zeros(n,2);
all_fields = cell(2, nsess);
region = {'HPC', 'ACC'};
fstart = 'D:\GarrettBlair\APA\HPCACC34990\k99_data\SMART_IL_4-8_';
% fname = 'D:\GarrettBlair\APA\HPCACC34990\k99_data\SMART_IL_4-8_HPC.mat'
figure(12345); clf
clims = [0, 1];
for ii = 1:2
    cmap = cmap_all{ii}; % temp.cellmap(:,sessns);
    shared = sum(cmap>0, 2) == nsess; %  & cmap(:,i)>0
    spks_all = [];
    time_all = [];
    corr_all = [];
    sessnum = [];
    
    for jj = 1:nsess
        matdat = matdat_all{ii,jj}; % load(data_fn);
        cmap_sub = cmap(shared, jj);
        if use_simple == false
            s    = matdat.ms.spks(cmap_sub, :);
            craw = matdat.ms.neuron.C(cmap_sub, :) + matdat.ms.neuron.YrA(cmap_sub, :);
            t    = matdat.ms.timestamps - matdat.ms.timestamps(1);
            t_dt    = matdat.ms.dt';
            cm   = matdat.ms.zonestruct.choices.paths.linear_vec;
            e    = matdat.ms.zonestruct.choices.paths.error==true;
            cm(e) = cm(e)*-1;
            lm   = matdat.ms.zonestruct.choices.paths.linear_vec;
            e    = matdat.ms.zonestruct.choices.paths.right==true;
            lm(e) = lm(e)*-1;
        else
            s    = matdat.spks(cmap_sub, :);
            craw = matdat.craw(cmap_sub, :);
            t    = matdat.time_ms - matdat.time_ms(1);
            t_dt    = matdat.dt';
            cm   = matdat.correct_linvec;
            lm   = matdat.left_linvec;
        end
        %         s = matdat.ms.neuron.S_matw(cmap_sub, :);
        %         sm = movmean(s, movmean_k, 2);
        %         sm = diff(cat(2, sm(:,1), sm), [], 2);
        %         cm = movmean(cm, movmean_k, 2);
%         sm = s;
        if dt>0
%             [sm, ~,~]      = average_spks_time(sm, dt, t, false, 'sum');
%             [cm, ~,timeout] = average_spks_time(cm, dt, t, false, 'mean');
%             [lm, ~,~] = average_spks_time(lm, dt, t, false, 'mean');

            if movmean_k>0; s = movmean(s, movmean_k, 2); end
            if movmean_k>0; craw = movmean(craw, movmean_k, 2); end
            t_ds = 0:dt:t_max(ii,jj);
            sm = zeros(size(s,1), length(t_ds));
            craw_ds = zeros(size(s,1), length(t_ds));
            for ss = 1:size(sm,1)
                sm(ss,:) = interp1(t, s(ss,:), t_ds, 'linear');
                craw_ds(ss,:) = interp1(t, craw(ss,:), t_ds, 'linear');
            end
                cm = interp1(t, cm, t_ds, 'linear', 'extrap');
                lm = interp1(t, lm, t_ds, 'linear', 'extrap');
%             t_ds = 0:dt:dt*(size(sm,2)-1);
            t_dt = t_ds*0 + dt;
        else
            cm   = cm';
            lm   = lm';
            t_ds = t' - t(1);
%             t_dt = t_ds*0 + dt;
        end
        if jj == 1
            spks_all = sm;
            craw_all = craw_ds;
            time_all = t_ds;
            dt_all   = t_dt;
            corr_all = cm;
            left_all = lm;
            sessnum  = cm*0 + jj;
        else
            spks_all = cat(2, spks_all, sm);
            craw_all = cat(2, craw_all, craw_ds);
            time_all = cat(2, time_all, t_ds);
            dt_all   = cat(2, dt_all,   t_dt);
            corr_all = cat(2, corr_all, cm);
            left_all = cat(2, left_all, lm);
            sessnum  = cat(2, sessnum, cm*0 + jj);
        end
    end
    if true
    spks_all(isnan(spks_all))=0;
%     spks_all = zscore(spks_all, 0, 2);
    
    timefield = zeros(size(spks_all,1), length(bins)-1);
    fields = zeros(size(spks_all,1), length(bins)-1);
    for ss = 1:size(spks_all,1)
        s = double(spks_all(ss,:));
        timefield(ss,:) = binned_statistic1d(corr_all, dt_all, bins, 'sum');
        fields(ss,:) = binned_statistic1d(corr_all, s, bins, 'sum');
%         timefield(ss,:) = binned_statistic1d(left_all, dt_all, bins, 'sum');
%         fields(ss,:) = binned_statistic1d(left_all, s, bins, 'sum');
    end
    fields = fields./timefield;
    all_fields{ii, jj} = fields;
    fields = normalize_rows(fields);
    [~, m] = max(fields, [], 2);
    [~, ord] = sort(m, 'descend');
%     fields = normalize_rows(fields);
    fields = fields(ord,:);
    fields_correct = fields(:,n+1+offset:end);
        
    fields_incorrect = fields(:,n:-1:1+offset);
    for kk = 1:size(lc, 2)
%         fieldsdist(kk,ii) = mean( pdist2(fields_incorrect(:,kk)', fields_correct(:,kk)', 'euclidean') );
        fieldsdist(kk,ii) = corr(fields_incorrect(:,kk), fields_correct(:,kk));
    end
    
%     latent_err = sc(1:3,:)';
    figure(12345)
    subplot(3,2,2*(ii-1)+1)
    imagesc(fields_correct, clims)
    nnn = (floor(sum(shared)/100));
    set(gca, 'XTick', [0, n/2, n-1]+.5, 'XTickLabel', [0,50,100],...
            'YTick', [linspace(1, 100*nnn, nnn+1)]-.5,... 
            'YTickLabel', [floor(linspace(1, 100*nnn, nnn+1))])
    subplot(3,2,2*(ii-1)+2)
    imagesc(fields_incorrect, clims)
    set(gca, 'XTick', [0, n/2, n-1]+.5, 'XTickLabel', [0,50,100],...
            'YTick', [linspace(1, 100*nnn, nnn+1)]-.5,... 
            'YTickLabel', [])
    colormap viridis
    if false
        [coeff, sc, latent] = pca(fields_correct');
        latent_correct = sc(:, 1:3);
        %     latent_correct = sc(1:3,:)';
        [coeff, sc, latent] = pca( fields_incorrect' );
        latent_err = sc(:, 1:3);
        subplot(2,2,3)
        plot(latent_correct)
        subplot(2,2,4)
        plot(latent_err)
        lc(ii,1+offset:end,:) = latent_correct;
        le(ii,1+offset:end,:) = latent_err;
        %
        figure(pcf)
        if true
            hold on
            scatter3(latent_correct(1,1), latent_correct(1,2), latent_correct(1,3), 'kx')
            plot3(latent_correct(:,1), latent_correct(:,2), latent_correct(:,3), 'g', 'LineWidth', 2)
            scatter3(latent_err(1,1), latent_err(1,2), latent_err(1,3), 'kx')
            plot3(latent_err(:,1), latent_err(:,2), latent_err(:,3), 'r', 'LineWidth', 2)
        else
            subplot(121)
            hold on
            scatter(latent_correct(1,1), latent_correct(1,2), 'gx')
            plot(latent_correct(:,1), latent_correct(:,2), 'g', 'LineWidth', 2)
            scatter(latent_err(1,1), latent_err(1,2), 'rx')
            plot(latent_err(:,1), latent_err(:,2), 'r', 'LineWidth', 2)
            
            subplot(122)
            hold on
            scatter(latent_correct(1,3), latent_correct(1,2), 'gx')
            plot(latent_correct(:,3), latent_correct(:,2), 'g', 'LineWidth', 2)
            scatter(latent_err(1,3), latent_err(1,2), 'rx')
            plot(latent_err(:,3), latent_err(:,2), 'r', 'LineWidth', 2)
            
        end
    end
    end
    %%
    if false
        v2 = Ziv_LaplacianEigenVectors(spks_all, true, [.075/15, .075]);
        v = v2{3}(:,2:4);
        figure()
        idx = corr_all<0;
        v = v2{3}(:,2:4)-.01;
        scatter3(v(idx,1), v(idx,2), v(idx,3), 5, -1*corr_all(idx), 'Marker', '.') 
        hold on
        idx = corr_all>0;
        v = v2{3}(:,2:4)+.01;
        scatter3(v(idx,1), v(idx,2), v(idx,3), 5, corr_all(idx), 'Marker', '.') 
        drawnow()
        disp(' ')
    end
    fname = sprintf('%s%s.mat', fstart, region{ii});
    vars2save = {'spks_all', 'craw_all', 'time_all', 'corr_all', 'left_all', 'sessnum', 'cmap'};
    for jj = 1:nsess
        eval(sprintf('sessname%d = data_fn{ii,jj};', jj));
        vars2save{end+1} = sprintf('sessname%d', jj);
    end
    if false; save(fname, vars2save{:}); end
end
figure(12345)
subplot(3,1,3); cla; hold on
plot(fieldsdist(:,1), 'b')
plot(fieldsdist(:,2), 'm')
set(gca, 'XTick', [0, n/2, n], 'XTickLabel', [0,50,100], 'YTick', [0, .5, 1], 'YLim', [-.1 1.1])

if false
%%% data set for HARSHA
h = load('SMART_IL_4-8_HPC.mat');
a = load('SMART_IL_4-8_ACC.mat');
hpc_spks    = h.spks_all;
acc_spks    = a.spks_all;
hpc_dff     = h.craw_all;
acc_dff     = a.craw_all;
sess_num    = a.sessnum;
time_ms     = a.time_all;
time_full   = 0:dt:dt*(length(time_ms)-1);
is_correct  = a.corr_all;
is_left     = a.left_all;

% hpc_spks  - filtered calcium signal (low passed, 1.5 z-score threshold)
% acc_spks  - filtered calcium signal (low passed, 1.5 z-score threshold)
% hpc_dff   - raw, denoised calcium fluorescence from hippocampus
% acc_dff   - raw, denoised calcium fluorescence from hippocampus
% sess_num  - (scalar 1-4) which session this data was from before concatenation
% time_ms   - (milliseconds) within session time; resets to 0 between sessions
% time_full - (milliseconds) within session time; resets to 0 between sessions   = 0:dt:dt*(length(time_ms)-1);
% is_correct - (0-1 scalar) marking the start to end (1) of a trial/lap. >0 is correct/rewarded, <0 incorrect/error
% is_left    - (0-1 scalar) marking the start to end (1) of a trial/lap. >0 is a left trajectory, <0 is right
vars2save = {'hpc_spks', 'acc_spks', 'hpc_dff', 'acc_dff', 'sess_num', 'time_ms', 'time_full', 'is_correct', 'is_left'};
fname = sprintf('%s%s.mat', fstart, 'HPC_ACC_aligned_nosmooth');
save(fname, vars2save{:})
end
%%
%
d= nan(n,2);
for ii = 1:n
    d(ii,1) = pdist2(squeeze(lc(1,ii,:))',  squeeze(lc(2,ii,:))', 'euclidean');
    d(ii,2) = pdist2(squeeze(le(1,ii,:))',  squeeze(le(2,ii,:))', 'euclidean');
%     d(ii,1) = sum(sqrt( (squeeze(lc(1,ii,:)) - squeeze(lc(2,ii,:))).^2 ));
%     d(ii,2) = sum(sqrt( (squeeze(le(1,ii,:)) - squeeze(le(2,ii,:))).^2 ));
end
figure()
plot(d)
%% K99 Bayesian decdoing example
correct_only = true;
spd = is_correct*0;
spd(2:end) = abs(diff(is_correct));
spd(1) = spd(2);
if correct_only
%     close all
    valid = is_correct>=0.1;
    pos  = is_correct(valid);
    spd = spd(valid)
    valid_bins = bins(2:end)>0;
else
    % valid = abs(is_correct)>0;
    % pos  = abs(is_correct(valid));
    valid = is_correct<=-0.1;
    pos  = -1*is_correct(valid);
    valid_bins = bins(2:end)<0;
end
fields_acc = normalize_rows(all_fields{2, 3}(:,valid_bins));
% fields_acc(isnan(fields_acc)) = 0;
[~, m] = max(fields_acc, [], 2);
[~, ord] = sort(m, 'descend');
acc_s = movmean(acc_dff(ord,:), 3, 2);
% acc_s = movmean(acc_spks(ord,:), 5, 2);
b_acc = fields_acc(ord,:)'*acc_s;
fields_hpc = normalize_rows(all_fields{1, 3}(:,valid_bins));
% fields_hpc(isnan(fields_hpc)) = 0;
[~, m] = max(fields_hpc, [], 2);
[~, ord] = sort(m, 'descend');
hpc_s = movmean(hpc_dff(ord,:), 3, 2);
% hpc_s = movmean(hpc_spks(ord,:), 5, 2);
b_hpc = fields_hpc(ord,:)'*hpc_s;

figure(101); clf; 
a = normalize_cols(b_acc(:,valid));
h = normalize_cols(b_hpc(:,valid));
[~, mh] = max(h, [], 1);
[~, ma] = max(a, [], 1);

pf = cat(1,a,h);
subplot_tight(2,1,1, [.05 .05]);
imagesc(h, [0, 1.1])
hold on
plot(pos*20, 'w.-', 'LineWidth', 1)
plot(mh, 'r.', 'LineWidth', 1)
xlim([0 565])
colorbar
% xlim([16960 17529])

subplot_tight(2,1,2, [.05 .05]);
imagesc(a, [0, 1.1])
hold on
plot(pos*20, 'w.-', 'LineWidth', 1)
plot(ma, 'r.', 'LineWidth', 1)
xlim([0 565])
% xlim([16960 17529])
colormap viridis

% valid = 
c = NaN(length(h),1);
for jj = 1:length(h)
    v = (isnan(h(:, jj)) == 0) & (isnan(a(:, jj)) == 0);
    if sum(v)>5
    c(jj) = corr(h(v, jj), a(v, jj));
    end
end

nb = 11;
bb = linspace(0.1,1,nb)
b = binned_statistic1d(pos, c, bb, 'nanmedian');

if correct_only
    figure(102); clf
    subplot_tight(2,1,2, [.05 .05]);
    subplot_tight(2,1,1, [.05 .05]);
    e1 = binned_statistic1d(pos, abs(pos-ma/size(a,1)), bb, 'nanmedian');
    e2 = binned_statistic1d(pos, abs(pos-mh/size(h,1)), bb, 'nanmedian');
    p = binned_statistic1d(pos, spd, bb, 'nanmedian');
else
    figure(102); 
    e1 = binned_statistic1d(pos, abs((pos)-ma/size(a,1)), bb, 'nanmedian');
    e2 = binned_statistic1d(pos, abs((1*pos)-mh/size(h,1)), bb, 'nanmedian');
end
subplot(211)
hold on; 
plot(e1, 'm')
plot(e2, 'b')
axis([ 0 11 0 .35])
subplot(212)
hold on; 
plot(b)
axis([ 0 11 .20 .8])
% yyaxis('right')
% plot(spd)
% corr
%%





function im = seg_im(cs, contour_thresh)

for k = 1:size(cs,1)
    ccs = squeeze(cs(k,:,:));
    ccs = ccs-min(ccs(:));
    ccs = ccs./max(ccs(:));
    ccs(ccs<=contour_thresh) = 0;
    cs(k,:,:) = ccs;
end
im = sum(cs,1).*1;

end