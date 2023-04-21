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
nantemplate = NaN(length(cellmat_Anum), 1);
% pred_pretr_eigs = nantemplate;
% pred_prepost_eigs = nantemplate;
% pred_pretr_spks = nantemplate;
% pred_prepost_spks = nantemplate;
vol_change_train = nantemplate;
vol_change_post = nantemplate;
train_in_pre = nantemplate;
post_in_pre = nantemplate;

x_fold_training = 2;
deltaT = .25; % second integration for spike binning
smooth_binning = false;
ndims = 10;
run_LapEig = true;
run_SVM = false;
eig_angle_corr = NaN(ndims, 4, length(cellmat_Anum));
%%
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
        spks_PRE_TR = frame_pre.deconv(cmap(segs, sess1), time1<=shocktimes(shocktime_anum == ANIMALNUM, 1))>0;
        spks_TR_PRE = frame_shock.deconv(cmap(segs, sess2), time2<=shocktimes(shocktime_anum == ANIMALNUM, 1))>0;
        bin_PRE_TR = bin_spks_time(spks_PRE_TR, deltaT, time1, smooth_binning);
        bin_TR_PRE = bin_spks_time(spks_TR_PRE, deltaT, time2, smooth_binning);
        if run_LapEig == true
            v1_iinds = [zeros(size(bin_PRE_TR,2), 1 ); ones(size(bin_TR_PRE,2), 1)];
            [v1, v1_valid] = Ziv_LaplacianEigenVectors([bin_PRE_TR, bin_TR_PRE], false);
            v1a = v1{3}(v1_iinds(v1_valid)==0, :);
            v1b = v1{3}(v1_iinds(v1_valid)==1, :);

        end
        if run_SVM == true
            sess_pre = zeros(size(v1a,1),1)+1;
            sess_tr = zeros(size(v1b,1),1)+2;
            all_spks = cat(1, v1a, v1b); all_label = cat(1, sess_pre, sess_tr);
            testlabel = mod(0:length(all_label)-1, x_fold_training) + 1;
            pred_label = all_label*NaN;
            for xx = 1:x_fold_training
                    rx = all_spks(testlabel~=xx,:); 
                    ry = all_label(testlabel~=xx);
%                     Mdl = fitcecoc(rx, ry);
                    Mdl = fitcsvm(rx, ry);
                    px = all_spks(testlabel==xx,:);
                    pred_label(testlabel==xx)     = predict(Mdl, px);
            end
            pred_pretr_eigs(aLoop) = sum(pred_label==all_label)/length(pred_label);
            %%%%%%%%%%%%%%%%%%%
            
            sess_pre = zeros(size(bin_PRE_TR,2),1)+1;
            sess_tr = zeros(size(bin_TR_PRE,2),1)+2;
            all_spks = cat(2, bin_PRE_TR, bin_TR_PRE)'; all_label = cat(1, sess_pre, sess_tr);
            testlabel = mod(0:length(all_label)-1, x_fold_training) + 1;
            pred_label = all_label*NaN;
            for xx = 1:x_fold_training
                    rx = all_spks(testlabel~=xx,:); 
                    ry = all_label(testlabel~=xx);
                    Mdl = fitcsvm(rx, ry);
                    px = all_spks(testlabel==xx,:);
                    pred_label(testlabel==xx)     = predict(Mdl, px);
            end
            pred_pretr_spks(aLoop) = sum(pred_label==all_label)/length(pred_label);
        end
        %             figure(24); clf; hold on; plot(all_label,'k'); plot(pred_label,'r.'); ylim([0 3])
%             title(sprintf('BINSPKS hit rate: %1.3f', pred_pretr_spks(aLoop)))
            drawnow
        
        %     matched = cmap(:, [sess2 | sess3]);
        matched = cmap(:, [sess1 | sess3]);
        segs = sum(matched>0, 2) == size(matched,2); %matched(:,1)>0 & matched(:,2)>0;
        %     spks_PRE_POST = frame_shock.deconv(cmap(segs, sess2), time2<=shocktimes(shocktime_anum == ANIMALNUM, 1));
        %%%%%%%%% using pre shock session to avoid scopo contamination
        spks_PRE_POST = frame_pre.deconv(cmap(segs, sess1), time1<=shocktimes(shocktime_anum == ANIMALNUM, 1))>0;
        %%%%%%%%%
        spks_POST_PRE = frame_post.deconv(cmap(segs, sess3), time3<=shocktimes(shocktime_anum == ANIMALNUM, 1))>0;
        
        %     popcorr_prob((eye(size(popcorr))==1)) = NaN;
        %     pos3 = average_spks_time(frame_post.posbin', 1, time3, false, 'median');
        %     v = Ziv_LaplacianEigenVectors([spks_PRE_POST, spks_POST_PRE], true);
        %     temp = Ziv_TopoClustering( v{3}, frame_post.posbin);
        bin_PRE_POST = bin_spks_time(spks_PRE_POST, deltaT, time1, smooth_binning);
        bin_POST_PRE = bin_spks_time(spks_POST_PRE, deltaT, time3, smooth_binning);
                
        if run_LapEig == true
            v2_iinds = [zeros(size(bin_PRE_POST,2), 1 ); ones(size(bin_POST_PRE,2), 1)];
            [v2, v2_valid] = Ziv_LaplacianEigenVectors([bin_PRE_POST, bin_POST_PRE], false);
            v2a = v2{3}(v2_iinds(v2_valid)==0, :);
            v2b = v2{3}(v2_iinds(v2_valid)==1, :);
            
            aabb = alphaShape(v1{3}(:,2), v1{3}(:,3), v1{3}(:,4));
            hhgg = alphaShape(v2{3}(:,2), v2{3}(:,3), v2{3}(:,4));
%             aabb = alphaShape(v1a(:,2), v1a(:,3), v1a(:,4));
%             hhgg = alphaShape(v2a(:,2), v2a(:,3), v2a(:,4));
            a = criticalAlpha(aabb,'one-region');
            h = criticalAlpha(hhgg,'one-region');
            aa = alphaShape(v1a(:,2), v1a(:,3), v1a(:,4), a);
%             a = criticalAlpha(aa,'one-region');
% %             a = .005;
%             aa = alphaShape(v1a(:,2), v1a(:,3), v1a(:,4), a);
            
            %             bb = alphaShape(v1b(:,2), v1b(:,3), v1b(:,4));
            %             a = criticalAlpha(bb,'one-region');
%             bb = alphaShape(v1b(:,2), v1b(:,3), v1b(:,4), a);
%             bb = alphaShape(v1b(:,2), v1b(:,3), v1b(:,4), 'RegionThreshold', .0002);
            bb = alphaShape(v1b(:,2), v1b(:,3), v1b(:,4), a);

            
            hh = alphaShape(v2a(:,2), v2a(:,3), v2a(:,4), h);
%             a = criticalAlpha(hh,'one-region');
%             hh = alphaShape(v2a(:,2), v2a(:,3), v2a(:,4), a);
%             gg = alphaShape(v2b(:,2), v2b(:,3), v2b(:,4), a);        
%             gg = alphaShape(v2b(:,2), v2b(:,3), v2b(:,4), a);
            gg = alphaShape(v2b(:,2), v2b(:,3), v2b(:,4), h);%,h);
            
            maxx = max([max(aa.Points(:)) max(bb.Points(:)) max(hh.Points(:)) max(gg.Points(:))]);
            minn = min([min(aa.Points(:)) min(bb.Points(:)) min(hh.Points(:)) min(gg.Points(:))]);
            emax = max(abs([maxx minn]));
            
            maxx = [max(aa.Points(:)) max(bb.Points(:)) max(hh.Points(:)) max(gg.Points(:))];
            minn = [min(aa.Points(:)) min(bb.Points(:)) min(hh.Points(:)) min(gg.Points(:))];

            vol_change_train(aLoop) = volume(bb)/volume(aa);
            vol_change_post(aLoop) = volume(gg)/volume(hh);

            qp1 = rand(50000,3);
            for q_loop = 1:3
                qmin = min([v1a(:,q_loop+1); v1b(:,q_loop+1)]);
                qmax = max([v1a(:,q_loop+1); v1b(:,q_loop+1)]);
                qp1(:,q_loop) = qmin + (qmax-qmin).*qp1(:,q_loop);
            end
            isinshape1 = inShape(aa, qp1);
            isinshape2 = inShape(bb, qp1);
            train_in_pre(aLoop) = sum(isinshape1==1 & isinshape2==1)/sum(isinshape1==1 | isinshape2==1);
            
            qp2 = rand(50000,3);
            for q_loop = 1:3
                qmin = min([v2a(:,q_loop+1); v2b(:,q_loop+1)]);
                qmax = max([v2a(:,q_loop+1); v2b(:,q_loop+1)]);
                qp2(:,q_loop) = qmin + (qmax-qmin).*qp2(:,q_loop);
            end
            isinshape1 = inShape(hh, qp2);
            isinshape2 = inShape(gg, qp2);
            post_in_pre(aLoop)= sum(isinshape1==1 & isinshape2==1)/sum(isinshape1==1 | isinshape2==1);

            n = ['Hipp' num2str(ANIMALNUM) '_linear' num2str(s1) '   sesstype_' num2str(sesslabel(aLoop))];
            figure(100+aLoop); clf; 
            set(gcf, 'Position', [114   161   645   596], 'Name', n);
            subplot_tight(2,3,1);  hold on; 
            plot3(v1a(:,2), v1a(:,3), v1a(:,4), 'r.'); 
            aa.plot('FaceColor', 'r', 'FaceAlpha', .2, 'LineStyle', 'none')
            axis([-emax emax -emax emax -emax emax]); view([31.3 53])
            
            subplot_tight(2,3,2);  hold on;
            plot3(v1b(:,2), v1b(:,3), v1b(:,4), 'b.')
            bb.plot('FaceColor', 'b', 'FaceAlpha', .2, 'LineStyle', 'none')
            axis([-emax emax -emax emax -emax emax]); view([31.3 53])

            subplot_tight(2,3,3);  hold on;
            aa.plot('FaceColor', 'r', 'FaceAlpha', .2, 'LineStyle', 'none')
            bb.plot('FaceColor', 'b', 'FaceAlpha', .2, 'LineStyle', 'none')
%             plot3(qp1(:,1), qp1(:,2), qp1(:,3), 'k.')
            axis(1.*[-emax emax -emax emax -emax emax]); view([31.3 53])
            title(sprintf('%0.2f', train_in_pre(aLoop)))
            
            subplot_tight(2,3,4);  hold on; 
            plot3(v2a(:,2), v2a(:,3), v2a(:,4), 'r.'); hold on
            hh.plot('FaceColor', 'r', 'FaceAlpha', .2, 'LineStyle', 'none')
            axis([-emax emax -emax emax -emax emax]); view([31.3 53])
            
            subplot_tight(2,3,5);  hold on; 
            plot3(v2b(:,2), v2b(:,3), v2b(:,4), 'b.');
            gg.plot('FaceColor', 'b', 'FaceAlpha', .2, 'LineStyle', 'none')
            axis([-emax emax -emax emax -emax emax]); view([31.3 53])
            
            subplot_tight(2,3,6);  hold on; 
            hh.plot('FaceColor', 'r', 'FaceAlpha', .2, 'LineStyle', 'none')
            gg.plot('FaceColor', 'b', 'FaceAlpha', .2, 'LineStyle', 'none')
%             plot3(qp2(:,1), qp2(:,2), qp2(:,3), 'k.')
            axis([-emax emax -emax emax -emax emax]); view([31.3 53])
            title(sprintf('%0.2f', post_in_pre(aLoop)))
            drawnow

            %%%%%%%%%%%%%%%%%%%
        end
%             figure(29); clf; hold on;
%             hh.plot('FaceColor', 'r', 'FaceAlpha', .2, 'LineStyle', 'none')
%             gg.plot('FaceColor', 'b', 'FaceAlpha', .2, 'LineStyle', 'none')
        if run_SVM == true
            sess_pre = zeros(size(v2a,1),1)+1;
            sess_post = zeros(size(v2b,1),1)+2;
            all_spks = cat(1, v2a, v2b); all_label = cat(1, sess_pre, sess_post);
            testlabel = mod(0:length(all_label)-1, x_fold_training) + 1;
            pred_label = all_label*NaN;
            for xx = 1:x_fold_training
                    rx = all_spks(testlabel~=xx,:); 
                    ry = all_label(testlabel~=xx);
%                     Mdl = fitcecoc(rx, ry);
                    Mdl = fitcsvm(rx, ry);
                    px = all_spks(testlabel==xx,:);
                    pred_label(testlabel==xx)     = predict(Mdl, px);
            end
            pred_prepost_eigs(aLoop) = sum(pred_label==all_label)/length(pred_label);
%             figure(25); clf; hold on; plot(all_label,'k'); plot(pred_label,'r.'); ylim([0 3])
%             title(sprintf('EIGS hit rate: %1.3f', pred_prepost_eigs(aLoop)))
            
            sess_pre = zeros(size(bin_PRE_POST,2),1)+1;
            sess_tr = zeros(size(bin_POST_PRE,2),1)+2;
            all_spks = cat(2, bin_PRE_POST, bin_POST_PRE)'; all_label = cat(1, sess_pre, sess_tr);
            testlabel = mod(0:length(all_label)-1, x_fold_training) + 1;
            pred_label = all_label*NaN;
            for xx = 1:x_fold_training
                    rx = all_spks(testlabel~=xx,:); 
                    ry = all_label(testlabel~=xx);
%                     Mdl = fitcecoc(rx, ry);
                    Mdl = fitcsvm(rx, ry);
                    px = all_spks(testlabel==xx,:);
                    pred_label(testlabel==xx)     = predict(Mdl, px);
            end
            pred_prepost_spks(aLoop) = sum(pred_label==all_label)/length(pred_label);
        end
        
    else
        if ~isfile(f1); fprintf('%s\n', f1); end
        if ~isfile(f2); fprintf('%s\n', f2); end
        if ~isfile(f3); fprintf('%s\n', f3); end
    end
end
% corr_pre(14) = NaN; % repeat session for pre
% overlap_pre(14) = NaN; % repeat session for pre

%%
figure(122); clf;
subplot(121)
hold on; 
plot([vol_change_train(sesslabel==1) vol_change_post(sesslabel==1)]', 'r')
plot([vol_change_train(sesslabel==0) vol_change_post(sesslabel==0)]', 'm')


% figure(124); clf;
subplot(122)
hold on; 
plot([train_in_pre(sesslabel==1) post_in_pre(sesslabel==1)]', 'r')
plot([train_in_pre(sesslabel==0) post_in_pre(sesslabel==0)]', 'm')
% Blair_et_al_cellcorr_plotting


% pred = [pred_pretr_spks pred_prepost_spks];
% figure(125); clf;
% plot(pred(sesslabel==1,:)', 'r')
% hold on; 
% plot(pred(sesslabel==0,:)', 'm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% ALL shared %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load matching data
nantemplate = NaN(length(cellmat_Anum), 1);
pred_pretr_eigs = nantemplate;
pred_prepost_eigs = nantemplate;
pred_pretr_spks = nantemplate;
pred_prepost_spks = nantemplate;
vol_change_train = nantemplate;
vol_change_post = nantemplate;
train_in_pre = nantemplate;
post_in_pre = nantemplate;
train_pre_dist = nantemplate;
train_post_dist = nantemplate;

x_fold_training = 5;
deltaT = .25; % second integration for spike binning
smooth_binning = false;
ndims = 10;
run_LapEig = true;
run_SVM = true;
eig_angle_corr = NaN(ndims, 4, length(cellmat_Anum));
%
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
        
        matched = cmap(:, [sess1 | sess2 | sess3]);
        segs = sum(matched>0, 2) == size(matched,2); %matched(:,1)>0 & matched(:,2)>0;
%         spks_PRE_TR = frame_pre.deconv(cmap(segs, sess1), time1<=shocktimes(shocktime_anum == ANIMALNUM, 1));
%         spks_TR_PRE = frame_shock.deconv(cmap(segs, sess2), time2<=shocktimes(shocktime_anum == ANIMALNUM, 1));
%         spks_POST_PRE = frame_post.deconv(cmap(segs, sess3), time3<=shocktimes(shocktime_anum == ANIMALNUM, 1));
        spks_PRE_TR = frame_pre.deconv(cmap(segs, sess1), time1<=shocktimes(shocktime_anum == ANIMALNUM, 1))>0;
        spks_TR_PRE = frame_shock.deconv(cmap(segs, sess2), time2<=shocktimes(shocktime_anum == ANIMALNUM, 1))>0;
        spks_PRE_POST = frame_pre.deconv(cmap(segs, sess1), time1<=shocktimes(shocktime_anum == ANIMALNUM, 1))>0;
        spks_POST_PRE = frame_post.deconv(cmap(segs, sess3), time3<=shocktimes(shocktime_anum == ANIMALNUM, 1))>0;
        bin_PRE_TR = bin_spks_time(spks_PRE_TR, deltaT, time1, smooth_binning);
        bin_TR_PRE = bin_spks_time(spks_TR_PRE, deltaT, time2, smooth_binning);
        bin_PRE_POST = bin_spks_time(spks_PRE_POST, deltaT, time1, smooth_binning);
        bin_POST_PRE = bin_spks_time(spks_POST_PRE, deltaT, time3, smooth_binning);

        v1_iinds = [zeros(size(bin_PRE_TR,2), 1 ); ones(size(bin_TR_PRE,2), 1); 1+ones(size(bin_POST_PRE,2), 1)];
        [v1, v1_valid] = Ziv_LaplacianEigenVectors([bin_PRE_TR, bin_TR_PRE, bin_POST_PRE], false);
        v1a = v1{3}(v1_iinds(v1_valid)==0, :);
        v1b = v1{3}(v1_iinds(v1_valid)==1, :);
        v1c = v1{3}(v1_iinds(v1_valid)==2, :);
        
        if run_SVM == true
            sess_pre = zeros(size(v1a,1),1)+1;
            sess_tr = zeros(size(v1b,1),1)+2;
            all_spks = cat(1, v1a, v1b); all_label = cat(1, sess_pre, sess_tr);
            testlabel = mod(0:length(all_label)-1, x_fold_training) + 1;
            pred_label = all_label*NaN;
            for xx = 1:x_fold_training
                    rx = all_spks(testlabel~=xx,:); 
                    ry = all_label(testlabel~=xx);
%                     Mdl = fitcecoc(rx, ry);
                    Mdl = fitcsvm(rx, ry);
                    px = all_spks(testlabel==xx,:);
                    pred_label(testlabel==xx)     = predict(Mdl, px);
            end
            pred_pretr_eigs(aLoop) = sum(pred_label==all_label)/length(pred_label);
            %%%%%%%%%%%%%%%%%%%
            
            sess_pre = zeros(size(bin_PRE_TR,2),1)+1;
            sess_tr = zeros(size(bin_TR_PRE,2),1)+2;
            all_spks = cat(2, bin_PRE_TR, bin_TR_PRE)'; all_label = cat(1, sess_pre, sess_tr);
            testlabel = mod(0:length(all_label)-1, x_fold_training) + 1;
            pred_label = all_label*NaN;
            for xx = 1:x_fold_training
                    rx = all_spks(testlabel~=xx,:); 
                    ry = all_label(testlabel~=xx);
%                     Mdl = fitcecoc(rx, ry);
                    Mdl = fitcsvm(rx, ry);
                    px = all_spks(testlabel==xx,:);
                    pred_label(testlabel==xx)     = predict(Mdl, px);
            end
            pred_pretr_spks(aLoop) = sum(pred_label==all_label)/length(pred_label);
        end
        %             figure(24); clf; hold on; plot(all_label,'k'); plot(pred_label,'r.'); ylim([0 3])
%             title(sprintf('BINSPKS hit rate: %1.3f', pred_pretr_spks(aLoop)))
            drawnow
%                     preds.predicted     = predict(Mdl, v1b);
        
        %     matched = cmap(:, [sess2 | sess3]);
%         matched = cmap(:, [sess1 | sess3]);
                
        if run_LapEig == true
%             v2_iinds = [zeros(size(bin_PRE_POST,2), 1 ); ones(size(bin_POST_PRE,2), 1)];
%             [v2, v2_valid] = Ziv_LaplacianEigenVectors([bin_PRE_POST, bin_POST_PRE], false);
%             v2a = v2{3}(v2_iinds(v2_valid)==0, :);
%             v2b = v2{3}(v2_iinds(v2_valid)==1, :);
            v2a = v1{3}(v1_iinds(v1_valid)==0, :);
            v2b = v1{3}(v1_iinds(v1_valid)==2, :);
            
% %             [v2hh, v2hh_valid] = Ziv_LaplacianEigenVectors([bin_PRE_POST], false);
% %             [v2gg, v2gg_valid] = Ziv_LaplacianEigenVectors([bin_POST_PRE], false);
% %             v2 = [];
% %             v2{3} = cat(1, v2hh{3}, v2gg{3});
% %             v2a = v2hh{3};
% %             v2b = v2gg{3};
            %%
%             a = .01;
%             aabb = alphaShape(v1{3}(:,2), v1{3}(:,3), v1{3}(:,4));
%             hhgg = alphaShape(v2{3}(:,2), v2{3}(:,3), v2{3}(:,4));
% % %             aabb = alphaShape(v1a(:,2), v1a(:,3), v1a(:,4));
% % %             hhgg = alphaShape(v2a(:,2), v2a(:,3), v2a(:,4));
%             a = criticalAlpha(aabb,'one-region');
%             h = criticalAlpha(hhgg,'one-region');
%             aa = alphaShape(v1a(:,2), v1a(:,3), v1a(:,4));
%             a = criticalAlpha(aa,'one-region');
% %             a = .005;
            aa = alphaShape(v1a(:,2), v1a(:,3), v1a(:,4));
            
            %             bb = alphaShape(v1b(:,2), v1b(:,3), v1b(:,4));
            %             a = criticalAlpha(bb,'one-region');
%             bb = alphaShape(v1b(:,2), v1b(:,3), v1b(:,4), a);
%             bb = alphaShape(v1b(:,2), v1b(:,3), v1b(:,4), 'RegionThreshold', .0002);
            bb = alphaShape(v1b(:,2), v1b(:,3), v1b(:,4));

            
            hh = alphaShape(v2a(:,2), v2a(:,3), v2a(:,4));
%             a = criticalAlpha(hh,'one-region');
%             hh = alphaShape(v2a(:,2), v2a(:,3), v2a(:,4), a);
%             gg = alphaShape(v2b(:,2), v2b(:,3), v2b(:,4), a);        
%             gg = alphaShape(v2b(:,2), v2b(:,3), v2b(:,4), a);
            gg = alphaShape(v2b(:,2), v2b(:,3), v2b(:,4));%,h);

            aam = median(aa.Points, 1);
            bbm = median(bb.Points, 1);
            hhm = median(hh.Points, 1);
            ggm = median(gg.Points, 1);
            train_pre_dist(aLoop) = sqrt(sum((aam - bbm).^2));
            train_post_dist(aLoop) = sqrt(sum((hhm - ggm).^2));
            
            maxx = max([max(aa.Points(:)) max(bb.Points(:)) max(hh.Points(:)) max(gg.Points(:))]);
            minn = min([min(aa.Points(:)) min(bb.Points(:)) min(hh.Points(:)) min(gg.Points(:))]);
            emax = max(abs([maxx minn]));
            
            maxx = [max(aa.Points(:)) max(bb.Points(:)) max(hh.Points(:)) max(gg.Points(:))];
            minn = [min(aa.Points(:)) min(bb.Points(:)) min(hh.Points(:)) min(gg.Points(:))];
%             gg = alphaShape(v2b(:,2), v2b(:,3), v2b(:,4));
%             a = criticalAlpha(gg,'one-region');
            vol_change_train(aLoop) = volume(bb)/volume(aa);
            vol_change_post(aLoop) = volume(gg)/volume(hh);
%             qp = 2*emax.*rand(50000,3) - emax;'
            qp1 = rand(50000,3);
            for q_loop = 1:3
                qmin = min([v1a(:,q_loop+1); v1b(:,q_loop+1)]);
                qmax = max([v1a(:,q_loop+1); v1b(:,q_loop+1)]);
                qp1(:,q_loop) = qmin + (qmax-qmin).*qp1(:,q_loop);
            end
            isinshape1 = inShape(aa, qp1);
            isinshape2 = inShape(bb, qp1);
            train_in_pre(aLoop) = sum(isinshape1==1 & isinshape2==1)/sum(isinshape1==1 | isinshape2==1);
            
            qp2 = rand(50000,3);
            for q_loop = 1:3
                qmin = min([v2a(:,q_loop+1); v2b(:,q_loop+1)]);
                qmax = max([v2a(:,q_loop+1); v2b(:,q_loop+1)]);
                qp2(:,q_loop) = qmin + (qmax-qmin).*qp2(:,q_loop);
            end
            isinshape1 = inShape(hh, qp2);
            isinshape2 = inShape(gg, qp2);
            post_in_pre(aLoop)= sum(isinshape1==1 & isinshape2==1)/sum(isinshape1==1 | isinshape2==1);

            
            n = ['Hipp' num2str(ANIMALNUM) '_linear' num2str(s1) '   sesstype_' num2str(sesslabel(aLoop))];
            figure(100+aLoop); clf; 
            set(gcf, 'Position', [114   161   645   596], 'Name', n);
            subplot_tight(2,3,1);  hold on; 
            scatter3(v1a(:,2), v1a(:,3), v1a(:,4), 'r.', 'MarkerEdgeAlpha', .2); 
            scatter3(aam(1), aam(2), aam(3), 300, 'k.'); 
            scatter3(aam(1), aam(2), aam(3), 100, 'ko'); 
            aa.plot('FaceColor', 'r', 'FaceAlpha', .2, 'LineStyle', 'none')
            axis([-emax emax -emax emax -emax emax]); view([31.3 53])
            
            subplot_tight(2,3,2);  hold on;
            scatter3(v1b(:,2), v1b(:,3), v1b(:,4), 'b.', 'MarkerEdgeAlpha', .2)
            scatter3(bbm(1), bbm(2), bbm(3), 300, 'k.'); 
            scatter3(bbm(1), bbm(2), bbm(3), 100, 'ko'); 
            bb.plot('FaceColor', 'b', 'FaceAlpha', .2, 'LineStyle', 'none')
            axis([-emax emax -emax emax -emax emax]); view([31.3 53])

            subplot_tight(2,3,3);  hold on;
            aa.plot('FaceColor', 'r', 'FaceAlpha', .2, 'LineStyle', 'none')
            bb.plot('FaceColor', 'b', 'FaceAlpha', .2, 'LineStyle', 'none')
            scatter3(aam(1), aam(2), aam(3), 300, 'r.'); 
            scatter3(bbm(1), bbm(2), bbm(3), 300, 'b.'); 
%             plot3(qp1(:,1), qp1(:,2), qp1(:,3), 'k.')
            axis(1.*[-emax emax -emax emax -emax emax]); view([31.3 53])
            title(sprintf('%0.2f', train_in_pre(aLoop)))
            
            subplot_tight(2,3,4);  hold on; 
            scatter3(v2a(:,2), v2a(:,3), v2a(:,4), 'r.', 'MarkerEdgeAlpha', .2); hold on
            scatter3(hhm(1), hhm(2), hhm(3), 300, 'k.'); 
            scatter3(hhm(1), hhm(2), hhm(3), 100, 'ko'); 
            hh.plot('FaceColor', 'r', 'FaceAlpha', .2, 'LineStyle', 'none')
            axis([-emax emax -emax emax -emax emax]); view([31.3 53])
            
            subplot_tight(2,3,5);  hold on; 
            scatter3(v2b(:,2), v2b(:,3), v2b(:,4), 'b.', 'MarkerEdgeAlpha', .2);
            scatter3(ggm(1), ggm(2), ggm(3), 300, 'k.'); 
            scatter3(ggm(1), ggm(2), ggm(3), 100, 'ko'); 
            gg.plot('FaceColor', 'b', 'FaceAlpha', .2, 'LineStyle', 'none')
            axis([-emax emax -emax emax -emax emax]); view([31.3 53])
            
            subplot_tight(2,3,6);  hold on; 
            hh.plot('FaceColor', 'r', 'FaceAlpha', .2, 'LineStyle', 'none')
            gg.plot('FaceColor', 'b', 'FaceAlpha', .2, 'LineStyle', 'none')
            scatter3(hhm(1), hhm(2), hhm(3), 300, 'r.'); 
            scatter3(ggm(1), ggm(2), ggm(3), 300, 'b.'); 
%             plot3(qp2(:,1), qp2(:,2), qp2(:,3), 'k.')
            axis([-emax emax -emax emax -emax emax]); view([31.3 53])
            title(sprintf('%0.2f', post_in_pre(aLoop)))
            drawnow

            %%%%%%%%%%%%%%%%%%%
        end
%             figure(29); clf; hold on;
%             hh.plot('FaceColor', 'r', 'FaceAlpha', .2, 'LineStyle', 'none')
%             gg.plot('FaceColor', 'b', 'FaceAlpha', .2, 'LineStyle', 'none')
        if run_SVM == true
            sess_pre = zeros(size(v2a,1),1)+1;
            sess_post = zeros(size(v2b,1),1)+2;
            all_spks = cat(1, v2a, v2b); all_label = cat(1, sess_pre, sess_post);
            testlabel = mod(0:length(all_label)-1, x_fold_training) + 1;
            pred_label = all_label*NaN;
            for xx = 1:x_fold_training
                    rx = all_spks(testlabel~=xx,:); 
                    ry = all_label(testlabel~=xx);
%                     Mdl = fitcecoc(rx, ry);
                    Mdl = fitcsvm(rx, ry);
                    px = all_spks(testlabel==xx,:);
                    pred_label(testlabel==xx)     = predict(Mdl, px);
            end
            pred_prepost_eigs(aLoop) = sum(pred_label==all_label)/length(pred_label);
%             figure(25); clf; hold on; plot(all_label,'k'); plot(pred_label,'r.'); ylim([0 3])
%             title(sprintf('EIGS hit rate: %1.3f', pred_prepost_eigs(aLoop)))
            
            sess_pre = zeros(size(bin_PRE_POST,2),1)+1;
            sess_tr = zeros(size(bin_POST_PRE,2),1)+2;
            all_spks = cat(2, bin_PRE_POST, bin_POST_PRE)'; all_label = cat(1, sess_pre, sess_tr);
            testlabel = mod(0:length(all_label)-1, x_fold_training) + 1;
            pred_label = all_label*NaN;
            for xx = 1:x_fold_training
                    rx = all_spks(testlabel~=xx,:); 
                    ry = all_label(testlabel~=xx);
%                     Mdl = fitcecoc(rx, ry);
                    Mdl = fitcsvm(rx, ry);
                    px = all_spks(testlabel==xx,:);
                    pred_label(testlabel==xx)     = predict(Mdl, px);
            end
            pred_prepost_spks(aLoop) = sum(pred_label==all_label)/length(pred_label);
        end
        
    else
        if ~isfile(f1); fprintf('%s\n', f1); end
        if ~isfile(f2); fprintf('%s\n', f2); end
        if ~isfile(f3); fprintf('%s\n', f3); end
    end
end
% ANIMALNUM == 15 recieved an immediate shock in the scop+shock training sess (23), so exclude since theres no pre shock data
% aLoop = 14
train_pre_dist(14) = NaN; % repeat session for pre
train_in_pre(14) = NaN; % repeat session for pre

svm_LapEig_decoding_plotting;
%
% figure(123); clf;
% subplot(121)
% hold on; 
% % plot([vol_change_train(sesslabel>=1) vol_change_post(sesslabel>=1)]', 'r')
% % plot([vol_change_train(sesslabel==0) vol_change_post(sesslabel==0)]', 'm')
% plot([train_pre_dist(sesslabel>=1) train_post_dist(sesslabel>=1)]', 'r')
% plot([train_pre_dist(sesslabel==0) train_post_dist(sesslabel==0)]', 'm')
% 
% 
% % figure(124); clf;
% subplot(122)
% hold on; 
% plot([train_in_pre(sesslabel>=1) post_in_pre(sesslabel>=1)]', 'r')
% plot([train_in_pre(sesslabel==0) post_in_pre(sesslabel==0)]', 'm')















