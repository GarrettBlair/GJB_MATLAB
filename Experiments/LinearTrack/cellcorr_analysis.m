%%
topdir = 'D:\Sample Data\doi_10.5068_D1ZT2S__v1\Blair_et_al._2022_repo\tadblair-main\Blair_et_al_DATA\'; % 'D:\MiniScopeData\Blair_et_al\'; %top level folder where data files are found
shockdir=[topdir 'shocktimes\']; %folder where the goodies are stored
cellmapdir=[topdir 'cellmaps\']; %folder where the goodies are stored
sessiondir=[topdir 'sessiondata\'];
% datadir=[topdir 'pretrn\'];

% cellmat_shockcols = [...
%     %    3 6 9;...  %Hipp6 3 day
%     5 6 8;...  %Hipp6 ext
%     6 7 10;... %Hipp7
%     3 4 7;...  %Hipp8
%     4 5 8;...  %Hipp9
%     8 9 12;... %Hipp15
%     5 6 9;...  %Hipp31
%     
%     5 6 9;... %Hipp12 (barrier)
%     5 6 9;... %Hipp13 (barrier)
%     6 7 8;... %Hipp18 (barrier) ext
%     7 8 11;... %Hipp35 (barrier)
%     5 6 7;... %Hipp30 (barrier) ext
%     7 8 11;... %Hipp34 (barrier)
%     
%     9 10 13;... %Hipp12 (shock)
%     9 10 13;... %Hipp13 (shock)
%     12 13 14;... %Hipp18 (shock) ext
%     11 12 15;...   %Hipp35 (shock)
%     
%     9 10 13;... %Hipp30 (shock+scop)
%     10 11 14;... %Hipp32 (shock+scop)
%     11 12 15;... %Hipp34 (shock+scop)
%     11 12 15;...   %Hipp36 (shock+scop)
%     
%     22 23 26;... %Hipp12 (shock+scop)
%     22 23 26;... %Hipp13 (shock+scop)
%     22 23 26;... %Hipp15 (shock+scop)
%     22 23 26;...   %Hipp18 (shock+scop)
%     19 20 21;...   %Hipp31 (shock+scop)
%     
%     13 14 15;... %Hipp30 (drug free shock2)
%     19 20 24;... %Hipp32 (drug free shock2)
%     
%     14 15 16;... %Hipp31 (scopolamine alone)
%     5 6 7;... %Hipp32 (scopolamine alone)
%     6 7 8]; %Hipp36 (scopolamine alone)
% sesslabel = [1 1 1 1 1 1 3 3 3 3 3 3 1 1 1 1 0 0 0 0 0 0 0 0 0 2 2 4 4 4];
% % 0==scopo+shock, 1==shock, 2==2nd shock, 3==barrier, 4==scopo alone
% cellmat_Anum = [6 7 8 9 15 31 12 13 18 35 30 34 12 13 18 35 30 32 34 36 12 13 15 18 31 30 32 31 32 36];
% is_shock  = [1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0];

% from fig4A_analysis.m
cellmat_shockcols = [...
    %    3 6 9;...  %Hipp6 3 day
    5 6 8;...  %r=1 Hipp6 ext
    6 7 11;... %r=2 Hipp7
    3 4 7;...  %r=3 Hipp8
    4 5 8;...  %r=4 Hipp9
    8 9 12;... %r=5 Hipp15
    5 6 9;...  %r=6 Hipp31

    5 6 9;... %r=7 Hipp12 (barrier)
    5 6 9;... %r=8 Hipp13 (barrier)
    6 7 9;... %r=9 Hipp18 (barrier) ext
    7 8 11;... %r=10 Hipp35 (barrier)
    5 6 9;... %r=11 Hipp30 (barrier) ext
    7 8 11;... %r=12 Hipp34 (barrier)

    9 10 13;... %r=13 Hipp12 (shock)
    9 10 13;... %r=14 Hipp13 (shock)
    12 13 16;... %r=15 Hipp18 (shock) ext
    11 12 15;...   %r=16 Hipp35 (shock)
    
    9 10 13;... %r=17 Hipp30 (shock+scop)
    10 11 14;... %r=18 Hipp32 (shock+scop)
    11 12 15;... %r=19 Hipp34 (shock+scop)
    11 12 15;...   %r=20 Hipp36 (shock+scop)
    
    22 23 26;... %r=21 Hipp12 (shock+scop)
    22 23 26;... %r=22 Hipp13 (shock+scop)
    22 22 26;... %r=23 Hipp15 (shock+scop)    %%%% 23 trainig sess --> 22
    22 23 26;...   %r=24 Hipp18 (shock+scop)
    19 20 23;...   %r=25 Hipp31 (shock+scop)
    
    13 14 17;... %r=26 Hipp30 (drug free shock2)
    19 20 21;... %r=27 Hipp32 (drug free shock2)
    
    14 15 16;... %r=28 Hipp31 (scopolamine alone)
    5 6 7;... %r=29 Hipp32 (scopolamine alone)
    6 7 8]; %r=30 Hipp36 (scopolamine alone)

sesslabel = [1 1 1 1 1 1 3 3 3 3 3 3 1 1 1 1 0 0 0 0 0 0 0 0 0 2 2 4 4 4];
% 0==scopo+shock, 1==shock, 2==2nd shock, 3==barrier, 4==scopo alone
is_shock  = [1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0];
cellmat_Anum = [6 7 8 9 15 31 12 13 18 35 30 34 12 13 18 35 30 32 34 36 12 13 15 18 31 30 32 31 32 36];

rats_to_analyze=[1:4 6 13:17 20:27]; %dfdf rats
% rats_to_analyze=[1:4 6:27];
cellmat_shockcols = cellmat_shockcols(rats_to_analyze, :);
sesslabel = sesslabel(rats_to_analyze);
is_shock = is_shock(rats_to_analyze);
cellmat_Anum = cellmat_Anum(rats_to_analyze);

load([shockdir 'shocktimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
load([shockdir 'shocktimes2']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
load([shockdir 'scopshocktimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
load([shockdir 'barriertimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
shocktime_anum = [6 7 8 9 12 13 15 18 30 31 32 34 35 36 37];


%% load matching data
corr_pre = NaN(length(cellmat_Anum),1); corr_post = corr_pre;
overlap_pre = corr_pre; overlap_post = corr_pre;
deltaT = 1; % second integration for spike binning
smooth_binning = false;
ndims = 10;
run_LapEig = false;
eig_angle_corr = NaN(ndims, 4, length(cellmat_Anum));
% ANIMALNUM == 15 recieved an immediate shock in the scop+shock training sess (23), so exclude since theres no pre shock data
for aLoop = 1:length(cellmat_Anum)
    %%
    clearvars sessionNums frame*
    ANIMALNUM = cellmat_Anum(aLoop);%6;
    s1 = cellmat_shockcols(aLoop, 1);
    s2 = cellmat_shockcols(aLoop, 2);
    s3 = cellmat_shockcols(aLoop, 3);
    switch sesslabel(aLoop)
        case 0
            load([cellmapdir '\Hipp' num2str(ANIMALNUM) '_scopshock_cmap.mat'])
        case 1
            load([cellmapdir '\Hipp' num2str(ANIMALNUM) '_shock_cmap.mat'])
        case 2
            load([cellmapdir '\Hipp' num2str(ANIMALNUM) '_shock2_cmap.mat'])
        case 3
            load([cellmapdir '\Hipp' num2str(ANIMALNUM) '_barrier_cmap.mat'])
        case 4
            load([cellmapdir '\Hipp' num2str(ANIMALNUM) '_scopo_cmap.mat'])
    end
    % load session data
    f1 = [sessiondir 'Hipp' num2str(ANIMALNUM) '_linear' num2str(s1) '_sess.mat'];
    f2 = [sessiondir 'Hipp' num2str(ANIMALNUM) '_linear' num2str(s2) '_sess.mat'];
    f3 = [sessiondir 'Hipp' num2str(ANIMALNUM) '_linear' num2str(s3) '_sess.mat'];
    
    plotting_sessions = (ANIMALNUM == 7 && sesslabel(aLoop)==1) || (ANIMALNUM == 36 && sesslabel(aLoop)==0);
    if all(isfile({f1, f2, f3}))
        %
    load(f1)
    load(f2)
    load(f3)
    frame_pre = eval(sprintf('frame%d', s1));
    frame_shock = eval(sprintf('frame%d', s2));
    frame_post = eval(sprintf('frame%d', s3));
    %
    sess1 = (sessionNums==s1);
    sess2 = (sessionNums==s2);
    sess3 = (sessionNums==s3);
    
    time1 = frame_pre.time./1000;
    time2 = frame_shock.time./1000;
    time3 = frame_post.time./1000;
    if ANIMALNUM == 15 && sesslabel(aLoop)==0
        % use the prveious session since shock was given too early
        shocktimes(shocktime_anum == ANIMALNUM, 1) = 600;
    end
    preshock_inds1 = time1<=shocktimes(shocktime_anum == ANIMALNUM, 1);
    preshock_inds2 = time2<=shocktimes(shocktime_anum == ANIMALNUM, 1);
    preshock_inds3 = time3<=shocktimes(shocktime_anum == ANIMALNUM, 1);
    time1 = time1(preshock_inds1);
    time2 = time2(preshock_inds2);
    time3 = time3(preshock_inds3);
%     matched = cmap(:, [sess1 | sess2 | sess3]);
    matched = cmap(:, [sess1 | sess2]);
    segs = sum(matched>0, 2) == size(matched,2); %matched(:,1)>0 & matched(:,2)>0;
    spks1 = frame_pre.deconv(cmap(segs, sess1), time1<=shocktimes(shocktime_anum == ANIMALNUM, 1));
    spks2 = frame_shock.deconv(cmap(segs, sess2), time2<=shocktimes(shocktime_anum == ANIMALNUM, 1));
    bin1 = bin_spks_time(spks1, deltaT, time1, smooth_binning);
    bin2 = bin_spks_time(spks2, deltaT, time2, smooth_binning);
    [corrpre1, p1] = corr(bin1');
    [corrpost1, p2] = corr(bin2');
    
    [~, ord1] = sort(corrpre1(1,:), 'descend');
    bs11 = [bin1(ord1,:)];
    bs12 = [bin2(ord1,:)];
    corrpre1 = corrpre1(ord1, :);
    corrpre1 = corrpre1(:,ord1);
    corrpost1 = corrpost1(ord1, :);
    corrpost1 = corrpost1(:,ord1);
%     figure(942); clf;
%     subplot(4,3,1:2); stacked_traces(normalize_rows(bs11(1:20,1:300)), .9);
%     subplot(4,3,4:5); stacked_traces(normalize_rows(bs12(1:20,1:300)), .9);
%     subplot(4,3,3); imagesc(corrpre1(1:20,1:20), [0 1]); axis square
%     subplot(4,3,6); imagesc(corrpost1(1:20,1:20), [0 1]); colormap hot; axis square
%     subplot(2,2,3); imagesc(corrpre1, [0 1]); axis square
%     subplot(2,2,4); imagesc(corrpost1, [0 1]); colormap hot; axis square
    if plotting_sessions==true
        clim = [-1 1];
    figure(900+aLoop); clf;
    set(gcf, 'Position', [124   161   745   730])
    subplot(4,3,1:2); stacked_traces(normalize_rows(bs11(1:20,1:600)), .9);
    title('Binned ca activity from PRE, sorted by corr with cell #1 in PRE session'); ylabel('Cell ID')
    subplot(4,3,4:5); stacked_traces(normalize_rows(bs12(1:20,1:600)), .9);
    title('Binned ca activity from TRAIN, sorted by corr with cell #1 in PRE'); xlabel('Time bin (sec)'); ylabel('Cell ID')
    subplot(4,3,3); imagesc(corrpre1(1:20,1:20), clim); axis square
    title('Corr mat of [left]');
    subplot(4,3,6); imagesc(corrpost1(1:20,1:20), clim); axis square
    title('Corr mat of [left]');
    subplot(2,2,3); imagesc(corrpre1, clim); axis square; colorbar
    title(sprintf('Full corr mat of PRE sorted by PRE'));
    subplot(2,2,4); imagesc(corrpost1, clim); colormap viridis; axis square; colorbar
    end
    
    prop_overlap = sum(sum(p1<=.05 & p2<=.05)) / sum(sum(p1<=.05 | p2<=.05));
        
%     bin2 = bin_spks_time(spks2, deltaT, time2, true);
%     bin1 = bin_spks_time(spks1, deltaT, time1, true);
% %     [popcorr_pre,  ~] = corr([bin1, bin2], 'type', 'Pearson');
% %     popcorr_pre((eye(size(popcorr_pre))==1)) = NaN;
%     popcorr_prob((eye(size(popcorr))==1)) = NaN;

    corrpre1(find(eye(size(corrpre1,1)))) = NaN;
    corrpost1(find(eye(size(corrpost1,1)))) = NaN;
    valid = ~isnan(corrpre1.*corrpost1);
    sig = p1<=.05 | p2<=.05;
%     diffmat = corr(corrpre1(valid & sig), corrpost1(valid & sig));
    diffmat = corr(corrpre1(valid), corrpost1(valid));
    if plotting_sessions==true
    title(sprintf('Full corr mat of TRAIN sorted by PRE\n PRE-TRAIN corr = %0.3f', diffmat));
    end
    if ANIMALNUM == 15 && sesslabel(aLoop)==0
        corr_pre(aLoop) = NaN;
        overlap_pre(aLoop) = NaN;
        shock_ang = NaN*time2;
    else
        corr_pre(aLoop) = diffmat;
        overlap_pre(aLoop) = prop_overlap;
    end
    if run_LapEig == true
        [pre_ang, ~]    = cart2pol(frame_pre.x(preshock_inds1),     frame_pre.y(preshock_inds1) - 25);
        if any(isnan(pre_ang))
            pre_ang(isnan(pre_ang)) = interp1(find(~isnan(pre_ang)), pre_ang(~isnan(pre_ang)), find(isnan(pre_ang)), 'nearest', 'extrap');
        end
        [shock_ang, ~]  = cart2pol(frame_shock.x(preshock_inds2),   frame_shock.y(preshock_inds2) - 25);
        if any(isnan(shock_ang))
            shock_ang(isnan(shock_ang)) = interp1(find(~isnan(shock_ang)), shock_ang(~isnan(shock_ang)), find(isnan(shock_ang)), 'nearest', 'extrap');
        end
        [post_ang, ~]   = cart2pol(frame_post.x(preshock_inds3),    frame_post.y(preshock_inds3) - 25);
        if any(isnan(post_ang))
            post_ang(isnan(post_ang)) = interp1(find(~isnan(post_ang)), post_ang(~isnan(post_ang)), find(isnan(post_ang)), 'nearest', 'extrap');
        end
        preds = [];
        preds.bins = [-pi:pi/24:pi]; % range fo binning angular position
        preds.bin_center = preds.bins(2:end) - abs(diff(preds.bins))/2;
        [preds.counts, ~,  preds.ang_binned] = histcounts(pre_ang, preds.bins);
        % setup the weights
        ang_weights = 1 - (preds.counts./sum(preds.counts));
        wmat = ones(length(pre_ang),1) * ang_weights;
        preds.weights = NaN(length(pre_ang), 1);
        for i = 1:length(preds.counts)
            ind = preds.ang_binned==i;
            preds.weights(ind) = wmat(ind, i);
        end
        
        bin1 = bin_spks_time(spks1, deltaT, time1, true);
        bin2 = bin_spks_time(spks2, deltaT, time2, true);
        v1_iinds = [zeros(size(bin1,2), 1 ); ones(size(bin2,2), 1)];
        [v1, v1_valid] = Ziv_LaplacianEigenVectors([bin1, bin2], false);
        v1a = v1{3}(v1_iinds(v1_valid)==0, :);
        v1b = v1{3}(v1_iinds(v1_valid)==1, :);
        
%         iii = preds.ang_binned ~= 1 & preds.ang_binned ~= 23;
%         Mdl = fitcecoc(v1a, preds.ang_binned',...
%             'Weights', preds.weights);
%         preds.predicted     = predict(Mdl, v1b);
%         preds.predicted     = predict(Mdl, v1b);        
    end
    
%     matched = cmap(:, [sess2 | sess3]);
    matched = cmap(:, [sess1 | sess3]);
    segs = sum(matched>0, 2) == size(matched,2); %matched(:,1)>0 & matched(:,2)>0;
%     spks2 = frame_shock.deconv(cmap(segs, sess2), time2<=shocktimes(shocktime_anum == ANIMALNUM, 1));
    %%%%%%%%% using pre shock session to avoid scopo contamination
    spks2 = frame_pre.deconv(cmap(segs, sess1), time1<=shocktimes(shocktime_anum == ANIMALNUM, 1));
    %%%%%%%%%
    spks3 = frame_post.deconv(cmap(segs, sess3), time3<=shocktimes(shocktime_anum == ANIMALNUM, 1));    

%     popcorr_prob((eye(size(popcorr))==1)) = NaN;
    %     pos3 = average_spks_time(frame_post.posbin', 1, time3, false, 'median');
    %     v = Ziv_LaplacianEigenVectors([spks2, spks3], true);
%     temp = Ziv_TopoClustering( v{3}, frame_post.posbin);
    bin2 = bin_spks_time(spks2, deltaT, time1, smooth_binning);
    bin3 = bin_spks_time(spks3, deltaT, time3, smooth_binning);
    [corrpre2, p2] = corr(bin2');
    [corrpost2, p3] = corr(bin3');
    prop_overlap = sum(sum(p2<=.05 & p3<=.05)) / sum(sum(p2<=.05 | p3<=.05));
    [~, ord2] = sort(corrpre2(1,:), 'descend');
    bs21 = [bin2(ord2,:)];
    bs22 = [bin3(ord2,:)];
    corrpre2 = corrpre2(ord2, :);
    corrpre2 = corrpre2(:,ord2);
    corrpost2 = corrpost2(ord2, :);
    corrpost2 = corrpost2(:,ord2);
    if plotting_sessions==true
    figure(901+aLoop); clf;
    set(gcf, 'Position', [924   161   745   730])
    subplot(4,3,1:2); stacked_traces(normalize_rows(bs21(1:20,1:600)), .9);
    title('Binned ca activity from PRE, sorted by corr with cell #1 in PRE session'); ylabel('Cell ID')
    subplot(4,3,4:5); stacked_traces(normalize_rows(bs22(1:20,1:600)), .9);
    title('Binned ca activity from POST, sorted by corr with cell #1 in PRE'); xlabel('Time bin (sec)'); ylabel('Cell ID')
    subplot(4,3,3); imagesc(corrpre2(1:20,1:20), clim); axis square
    title('Corr mat of [left]');
    subplot(4,3,6); imagesc(corrpost2(1:20,1:20), clim); axis square
    title('Corr mat of [left]');
    subplot(2,2,3); imagesc(corrpre2, clim); axis square; colorbar
    title(sprintf('Full corr mat of PRE sorted by PRE'));
    subplot(2,2,4); imagesc(corrpost2, clim); colormap viridis; axis square; colorbar
    end
    %     subplot(2,3,6); imagesc(corrpre2-corrpost2, [0 1]); axis square
%     bin2 = bin_spks_time(spks2, deltaT, time2, true);
%     bin3 = bin_spks_time(spks3, deltaT, time3, true);
    % next calc the correlation between time bins i and j
% %     [popcorr_post,  ~] = corr([bin2, bin3], 'type', 'Pearson');
% %     popcorr_post((eye(size(popcorr_post))==1)) = NaN;
% %     figure(aLoop); clf; 
% %     subplot(121); imagesc(popcorr_pre, [0 .6]); subplot(122); imagesc(popcorr_post, [0 .6]); colormap hot
    corrpre2(find(eye(size(corrpre2,1)))) = NaN;
    corrpost2(find(eye(size(corrpost2,1)))) = NaN;
    valid = ~isnan(corrpre2.*corrpost2);
    sig = p2<=.05 | p3<=.05;
%     diffmat = nansum(nansum((corrpre1 - corrpost1).^2))/sum(segs);
%     diffmat = corr(corrpre2(valid & sig), corrpost2(valid & sig));
%     diffmat = corr(corrpre2(valid & sig), corrpost2(valid & sig));
    diffmat = corr(corrpre2(valid), corrpost2(valid));
%     dd = pdist2(corrpre2(:), corrpost2(:), 'mahalanobis');
    if plotting_sessions==true % Hipp7 shock data
    title(sprintf('Full corr mat of POST sorted by PRE\n PRE-POST corr = %0.3f', diffmat));
    end
    %     dpost(aLoop) = nanmean(diffmat(:));
    corr_post(aLoop) = diffmat;%o;%nanmean(diffmat(:));
    overlap_post(aLoop) = prop_overlap;%o;%nanmean(diffmat(:));

    if run_LapEig == true
    bin2 = bin_spks_time(spks2, deltaT, time1, true);
    bin3 = bin_spks_time(spks3, deltaT, time3, true);
    v2_iinds = [zeros(size(bin2,2), 1 ); ones(size(bin3,2), 1)];
    [v2, v2_valid] = Ziv_LaplacianEigenVectors([bin2, bin3], false);
    v2a = v2{3}(v2_iinds(v2_valid)==0, :);
    v2b = v2{3}(v2_iinds(v2_valid)==1, :);

%     v2a = v2{3}(1:size(bin2,2), 2:4);
%     v2b = v2{3}(size(bin2,2)+1:end, 2:4);
%     c1a = corr(frame_pre.x(preshock_inds1), v1a)';
%     c1b = corr(frame_shock.x(preshock_inds1), v1b)';
%     c2a = corr(frame_pre.x(preshock_inds1), v2a)';
%     c2b = corr(frame_post.x(preshock_inds3), v2b)';
    c1a = corr(pre_ang(v1_valid(v1_iinds==0)), v1a)';
    c1b = corr(shock_ang(v1_valid(v1_iinds==1)), v1b)';
    c2a = corr(pre_ang(v2_valid(v2_iinds==0)), v2a)';
    c2b = corr(post_ang(v2_valid(v2_iinds==1)), v2b)';
    eig_angle_corr(:,:,aLoop) = [c1a c1b c2a c2b];
%     figure(aLoop); clf; 
%     subplot(211); cla; hold on
%     plot(v1{3}(:,2), 'k')
%     yyaxis('right')
%     plot(find(v1_iinds==0), pre_ang, 'b-')
%     plot(find(v1_iinds==1), shock_ang, 'r-')
% %     plot3(v1a(:,1), v1a(:,2), v1a(:,3), 'k.')
% %     plot3(v1b(:,1), v1b(:,2), v1b(:,3), 'r.')
% 
%     subplot(212); cla; hold on
%     plot(v2{3}(:,2), 'k')
%     yyaxis('right')
%     plot(find(v2_iinds==0), pre_ang, 'b-')
%     plot(find(v2_iinds==1), post_ang, 'r-')
% %     plot3(v2a(:,1), v2a(:,2), v2a(:,3), 'k.')
% %     plot3(v2b(:,1), v2b(:,2), v2b(:,3), 'r.')
    end

    else
        if ~isfile(f1); fprintf('%s\n', f1); end
        if ~isfile(f2); fprintf('%s\n', f2); end
        if ~isfile(f3); fprintf('%s\n', f3); end
    end
end
% corr_pre(14) = NaN; % repeat session for pre
% overlap_pre(14) = NaN; % repeat session for pre

%

Blair_et_al_cellcorr_plotting




