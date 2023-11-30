%% matching across sessions
dataDir = 'D:\GarrettBlair\APA\';

animals = {'HPCACC24500', 'HPCACC24502', 'HPCACC24503', 'Acc20832', 'Acc19947', 'Hipp18240'};% animals = {'Acc20832', 'Acc19947'};

matching_dir = 'matching_contours\manual_alignment\cellreg\';
hpc_matching_dir = 'matching_contours\manual_alignment_HPC\cellreg\';
acc_matching_dir = 'matching_contours\manual_alignment_ACC\cellreg\';

%%
disp(' ')
disp('Finding files...')
rewrite_matchmatrix = false;
for aLoop = 1:length(animals)
    %%
    ddir = [dataDir animals{aLoop} '\' matching_dir];
    processedfile_dir = [dataDir animals{aLoop} '\processed_files'];
    cmapdir = dir([ddir 'cellregistered*']);
    cmapfile = [cmapdir.folder '\' cmapdir.name];
    if isfile(cmapfile)
        disp(cmapfile)
        if contains(animals{aLoop}, 'Acc')==1 || contains(animals{aLoop}, 'ACC')==1
            region = 'ACC';
        elseif contains(animals{aLoop}, 'Hipp')==1 || contains(animals{aLoop}, 'HPC')==1
            region = 'HPC';
        end
        cmap_out = [processedfile_dir '\cellmatching_' region '.mat'];
        Setup_matching_marix([cmapdir.folder '\'], cmap_out, rewrite_matchmatrix)
    end
    ddir = [dataDir animals{aLoop} '\' hpc_matching_dir];
    cmapdir = dir([ddir 'cellregistered*']);
    cmapfile = [cmapdir.folder '\' cmapdir.name];
    if isfile(cmapfile)
        disp(cmapfile)
        region = 'HPC';
        cmap_out = [processedfile_dir '\cellmatching_' region '.mat'];
        Setup_matching_marix([cmapdir.folder '\'], cmap_out, rewrite_matchmatrix)
    end
    ddir = [dataDir animals{aLoop} '\' acc_matching_dir];
    cmapdir = dir([ddir 'cellregistered*']);
    cmapfile = [cmapdir.folder '\' cmapdir.name];
    if isfile(cmapfile)
        disp(cmapfile)
        region = 'ACC';
        cmap_out = [processedfile_dir '\cellmatching_' region '.mat'];
        Setup_matching_marix([cmapdir.folder '\'], cmap_out, rewrite_matchmatrix)
    end
end

%%

region = 'HPC';
numA = length(animals);
minsharedcells = 1;
max_day_sep = 11;
max_reps = 50;

HPC_pop_overlap = NaN(numA, 100, max_reps);
HPC_sess_types = cell(numA, 100, 100);
HPC_sep_idx_count = zeros(numA, 100);
HPC_cofire_change_time = NaN(numA, max_day_sep, max_reps);
HPC_cofiresig_change_time = NaN(numA, max_day_sep, max_reps);
HPC_pfieldcorr_room_time = NaN(numA, max_day_sep, max_reps);
HPC_pfieldcorr_arena_time = NaN(numA, max_day_sep, max_reps);
HPC_SFEPchange_time = NaN(numA, max_day_sep, max_reps);
%
for aLoop = 1:length(animals)
    %%
    ddir = [dataDir animals{aLoop} '\' matching_dir];
    processedfile_dir = [dataDir animals{aLoop} '\processed_files'];
    cmap_out = [processedfile_dir '\cellmatching_' region '.mat'];
    if isfile(cmap_out)
        cm = load(cmap_out);
        cellmap = cm.cellmap;
        [nsegs, nsess] = size(cm.cellmap);
        day_sep = NaN(nsess);
        %         HPC_cofire_change_an{aLoop} = NaN(nsess);
        for i = 1:nsess-1
            i_name = [processedfile_dir '\' cm.session_names{i}];
            [sessDate1, trialname1, trial_type1, trial_num1] = recdata_from_filename(i_name);
            for j = i+1:nsess
                %%
                j_name = [processedfile_dir '\' cm.session_names{j}];
                [sessDate2, trialname2, trial_type2, trial_num2] = recdata_from_filename(j_name);
                if strcmp(trial_type1, trial_type2)==true
                    if isempty(sess_i)
                        disp([' > ' i_name])
                        sess_i = load(i_name);
                    end
                    if isempty(sess_j)
                        disp([' ---> ' j_name])
                        sess_j = load(j_name);
                    end
                    disp('')
                    day_sep(i,j) = abs(round(days(sessDate2 - sessDate1)));
                    HPC_sess_types{aLoop, i, j} = [trial_type1 '-' trial_type2];
                    daysep_idx = day_sep(i,j) + 1;
                    HPC_sep_idx_count(aLoop, daysep_idx) = HPC_sep_idx_count(aLoop, daysep_idx)+1;
                    rep = HPC_sep_idx_count(aLoop, daysep_idx);
                    
                    sharedcells = cellmap(:,i)>0 & cellmap(:,j)>0;
                    totalcells  = cellmap(:,i)>0 | cellmap(:,j)>0;
                    HPC_pop_overlap(aLoop, daysep_idx, rep) = sum(sharedcells>0)/sum(totalcells>0);
                    if sum(sharedcells>0)>=minsharedcells && daysep_idx<=max_day_sep
                        
                        ns = sum(sharedcells>0);
                        segsi = cellmap(sharedcells,i);
                        segsj = cellmap(sharedcells,j);
                        
                        spksi = sess_i.ms.room.svm_decoding.spks_bin(segsi, :);
                        spksj = sess_j.ms.room.svm_decoding.spks_bin(segsj, :);
                        [cc1, pp1] = nancorr(spksi');
                        [cc2, pp2] = nancorr(spksj');
                        sigs = pp1<=.05 & pp2 <=.05;
                        mask = triu(ones(ns), 1)==1;
                        ccdiff = cc1-cc2;
                        ccdiff(~mask) = NaN;
                        HPC_cofire_change_time(aLoop, daysep_idx, rep) = nanmedian(abs(ccdiff(:)));
                        ccdiff(~sigs) = NaN;
                        HPC_cofiresig_change_time(aLoop, daysep_idx, rep) = nanmedian(abs(ccdiff(:)));
                        %                     HPC_cofire_change_an{aLoop}(i,j) = nanmedian(abs(ccdiff(:)));
                        
                        p1 = sess_i.ms.room.pfields_smooth(segsi, :,:);
                        p2 = sess_j.ms.room.pfields_smooth(segsj, :,:);
                        [d1, d2, d3] = size(p1);
                        p1 = reshape(p1, [d1, d2*d3]);
                        p2 = reshape(p2, [d1, d2*d3]);
                        c = NaN(d1,1); p = NaN(d1,1);
                        for ii = 1:d1
                            [c(ii), p(ii)] = nancorr(p1(ii,:),p2(ii,:));
                        end
                        HPC_pfieldcorr_room_time(aLoop, daysep_idx, rep) = nanmedian(c(:));
                        
                        p1 = sess_i.ms.arena.pfields_smooth(segsi, :,:);
                        p2 = sess_j.ms.arena.pfields_smooth(segsj, :,:);
                        [d1, d2, d3] = size(p1);
                        p1 = reshape(p1, [d1, d2*d3]);
                        p2 = reshape(p2, [d1, d2*d3]);
                        c = NaN(d1,1); p = NaN(d1,1);
                        for ii = 1:d1
                            [c(ii), p(ii)] = nancorr(p1(ii,:),p2(ii,:));
                        end
                        HPC_pfieldcorr_arena_time(aLoop, daysep_idx, rep) = nanmedian(c(:));
                        
                        iposi = abs(sess_i.ms.room.momentary_pos_info(segsi, :)) - abs(sess_i.ms.arena.momentary_pos_info(segsi, :));
                        iposj = abs(sess_j.ms.room.momentary_pos_info(segsj, :)) - abs(sess_j.ms.arena.momentary_pos_info(segsj, :));
                        cellSFEP_i = sum(iposi>0, 2)./sum(iposi>0 | iposi<0, 2);
                        cellSFEP_j = sum(iposj>0, 2)./sum(iposj>0 | iposj<0, 2);
                        HPC_SFEPchange_time(aLoop, daysep_idx, rep) = nanmedian(cellSFEP_i-cellSFEP_j);
                    end
                end % TR TR only
                sess_j = [];
            end
            sess_i = [];
        end
    end
end



%%
region = 'ACC';
numA = length(animals);
ACC_pop_overlap = NaN(numA, 100, max_reps);
ACC_sess_types = cell(numA, 100, 100);
ACC_sep_idx_count = zeros(numA, 100);
ACC_cofire_change_time = NaN(numA, max_day_sep, max_reps);
ACC_cofiresig_change_time = NaN(numA, max_day_sep, max_reps);
ACC_pfieldcorr_room_time = NaN(numA, max_day_sep, max_reps);
ACC_pfieldcorr_arena_time = NaN(numA, max_day_sep, max_reps);
ACC_SFEPchange_time = NaN(numA, max_day_sep, max_reps);
% ACC_cofire_change_an = cell(numA,1);
%
for aLoop = 1:length(animals)
    %%
    ddir = [dataDir animals{aLoop} '\' matching_dir];
    processedfile_dir = [dataDir animals{aLoop} '\processed_files'];
    cmap_out = [processedfile_dir '\cellmatching_' region '.mat'];
    if isfile(cmap_out)
        cm = load(cmap_out);
        cellmap = cm.cellmap;
        [nsegs, nsess] = size(cm.cellmap);
        day_sep = NaN(nsess);
        %         ACC_cofire_change_an{aLoop} = NaN(nsess);
        for i = 1:nsess-1
            i_name = [processedfile_dir '\' cm.session_names{i}];
            [sessDate1, trialname1, trial_type1, trial_num1] = recdata_from_filename(i_name);
            for j = i+1:nsess
                %%
                j_name = [processedfile_dir '\' cm.session_names{j}];
                [sessDate2, trialname2, trial_type2, trial_num2] = recdata_from_filename(j_name);
                if strcmp(trial_type1, trial_type2)==true
                    if isempty(sess_i)
                        disp([' > ' i_name])
                        sess_i = load(i_name);
                    end
                    if isempty(sess_j)
                        disp([' ---> ' j_name])
                        sess_j = load(j_name);
                    end
                    disp('')
                    day_sep(i,j) = abs(round(days(sessDate2 - sessDate1)));
                    ACC_sess_types{aLoop, i, j} = [trial_type1 '-' trial_type2];
                    daysep_idx = day_sep(i,j) + 1;
                    ACC_sep_idx_count(aLoop, daysep_idx) = ACC_sep_idx_count(aLoop, daysep_idx)+1;
                    rep = ACC_sep_idx_count(aLoop, daysep_idx);
                    
                    sharedcells = cellmap(:,i)>0 & cellmap(:,j)>0;
                    totalcells  = cellmap(:,i)>0 | cellmap(:,j)>0;
                    ACC_pop_overlap(aLoop, daysep_idx, rep) = sum(sharedcells>0)/sum(totalcells>0);
                    if sum(sharedcells>0)>=minsharedcells && daysep_idx<=max_day_sep
                        
                        ns = sum(sharedcells>0);
                        segsi = cellmap(sharedcells,i);
                        segsj = cellmap(sharedcells,j);
                        
                        spksi = sess_i.ms.room.svm_decoding.spks_bin(segsi, :);
                        spksj = sess_j.ms.room.svm_decoding.spks_bin(segsj, :);
                        [cc1, pp1] = nancorr(spksi');
                        [cc2, pp2] = nancorr(spksj');
                        sigs = pp1<=.05 & pp2 <=.05;
                        mask = triu(ones(ns), 1)==1;
                        ccdiff = cc1-cc2;
                        ccdiff(~mask) = NaN;
                        ACC_cofire_change_time(aLoop, daysep_idx, rep) = nanmedian(abs(ccdiff(:)));
                        ccdiff(~sigs) = NaN;
                        ACC_cofiresig_change_time(aLoop, daysep_idx, rep) = nanmedian(abs(ccdiff(:)));
                        %                     ACC_cofire_change_an{aLoop}(i,j) = nanmedian(abs(ccdiff(:)));
                        
                        p1 = sess_i.ms.room.pfields_smooth(segsi, :,:);
                        p2 = sess_j.ms.room.pfields_smooth(segsj, :,:);
                        [d1, d2, d3] = size(p1);
                        p1 = reshape(p1, [d1, d2*d3]);
                        p2 = reshape(p2, [d1, d2*d3]);
                        c = NaN(d1,1); p = NaN(d1,1);
                        for ii = 1:d1
                            [c(ii), p(ii)] = nancorr(p1(ii,:),p2(ii,:));
                        end
                        ACC_pfieldcorr_room_time(aLoop, daysep_idx, rep) = nanmedian(c(:));
                        
                        p1 = sess_i.ms.arena.pfields_smooth(segsi, :,:);
                        p2 = sess_j.ms.arena.pfields_smooth(segsj, :,:);
                        [d1, d2, d3] = size(p1);
                        p1 = reshape(p1, [d1, d2*d3]);
                        p2 = reshape(p2, [d1, d2*d3]);
                        c = NaN(d1,1); p = NaN(d1,1);
                        for ii = 1:d1
                            [c(ii), p(ii)] = nancorr(p1(ii,:),p2(ii,:));
                        end
                        ACC_pfieldcorr_arena_time(aLoop, daysep_idx, rep) = nanmedian(c(:));
                        
                        iposi = abs(sess_i.ms.room.momentary_pos_info(segsi, :)) - abs(sess_i.ms.arena.momentary_pos_info(segsi, :));
                        iposj = abs(sess_j.ms.room.momentary_pos_info(segsj, :)) - abs(sess_j.ms.arena.momentary_pos_info(segsj, :));
                        cellSFEP_i = sum(iposi>0, 2)./sum(iposi>0 | iposi<0, 2);
                        cellSFEP_j = sum(iposj>0, 2)./sum(iposj>0 | iposj<0, 2);
                        ACC_SFEPchange_time(aLoop, daysep_idx, rep) = nanmedian(cellSFEP_i-cellSFEP_j);
                    end
                end % TR TR only
                sess_j = [];
            end
            sess_i = [];
        end
    end
end
%%
startfig = 2000;
figure(startfig+1); clf; 
subplot(2,2,1)
imagesc(squeeze(nansum(HPC_pfieldcorr_room_time>0,1)))
subplot(2,2,2)
imagesc(squeeze(nansum(HPC_pfieldcorr_arena_time>0,1)))
subplot(2,2,3)
imagesc(squeeze(nansum(ACC_pfieldcorr_room_time>0,1)))
subplot(2,2,4)
imagesc(squeeze(nansum(ACC_pfieldcorr_arena_time>0,1)))

figure(startfig+2); clf; 
subplot(2,2,1)
imagesc(squeeze(nanmean(HPC_pfieldcorr_room_time,1)))
subplot(2,2,2)
imagesc(squeeze(nanmean(HPC_pfieldcorr_arena_time,1)))
subplot(2,2,3)
imagesc(squeeze(nanmean(ACC_pfieldcorr_room_time,1)))
subplot(2,2,4)
imagesc(squeeze(nanmean(ACC_pfieldcorr_arena_time,1)))

%%
figure(startfig+3); clf; 
hold on

aclr  = jet(numA);
x = 0:max_day_sep-1;
for aLoop = 1:numA
    y = nanmean(HPC_pop_overlap(aLoop,1:max_day_sep,:), 3);
    subplot(2,1,1); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
    
    y = nanmean(ACC_pop_overlap(aLoop,1:max_day_sep,:), 3);
    subplot(2,1,2); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
end

y = nanmean(nanmean(HPC_pop_overlap(:,1:max_day_sep,:),3), 1);
subplot(2,1,1); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep 0 1.1])
xlabel('days apart')
y = nanmean(nanmean(ACC_pop_overlap(:,1:max_day_sep,:),3), 1);
subplot(2,1,2); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep 0 1.1])
xlabel('days apart')
%%
figure(startfig+4); clf; 
hold on

aclr  = jet(numA);
x = 0:max_day_sep-1;
for aLoop = 1:numA
    y = nanmean(HPC_SFEPchange_time(aLoop,1:max_day_sep,:), 3);
    subplot(2,1,1); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
    
    y = nanmean(ACC_SFEPchange_time(aLoop,1:max_day_sep,:), 3);
    subplot(2,1,2); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
end

y = nanmean(nanmean(HPC_SFEPchange_time(:,1:max_day_sep,:),3), 1);
subplot(2,1,1); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep -.25 .25])

y = nanmean(nanmean(ACC_SFEPchange_time(:,1:max_day_sep,:),3), 1);
subplot(2,1,2); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep -.25 .25])

%%
figure(startfig+5); clf; 
hold on

aclr  = jet(numA);
x = 0:max_day_sep-1;
for aLoop = 1:numA
    y = nanmean(HPC_cofire_change_time(aLoop,1:max_day_sep,:), 3);
    subplot(2,2,1); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
    
    y = nanmean(HPC_cofiresig_change_time(aLoop,1:max_day_sep,:), 3);
    subplot(2,2,2); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
    y = nanmean(ACC_cofire_change_time(aLoop,1:max_day_sep,:), 3);
    subplot(2,2,3); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
    y = nanmean(ACC_cofiresig_change_time(aLoop,1:max_day_sep,:), 3);
    subplot(2,2,4); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
end

y = nanmean(nanmean(HPC_cofire_change_time(:,1:max_day_sep,:),3), 1);
subplot(2,2,1); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep 0 .03])

y = nanmean(nanmean(HPC_cofiresig_change_time(:,1:max_day_sep,:),3), 1);
subplot(2,2,2); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep 0 .1])

y = nanmean(nanmean(ACC_cofire_change_time(:,1:max_day_sep,:),3), 1);
subplot(2,2,3); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep 0 .03])

y = nanmean(nanmean(ACC_cofiresig_change_time(:,1:max_day_sep,:),3), 1);
subplot(2,2,4); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep 0 .1])


%%
figure(startfig+6); clf; 
hold on

aclr  = jet(numA);
x = 0:max_day_sep-1;
for aLoop = 1:numA
    y = nanmean(HPC_pfieldcorr_room_time(aLoop,:,:), 3);
    subplot(2,2,1); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
    y = nanmean(HPC_pfieldcorr_arena_time(aLoop,:,:), 3);
    subplot(2,2,2); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
    y = nanmean(ACC_pfieldcorr_room_time(aLoop,:,:), 3);
    subplot(2,2,3); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
    y = nanmean(ACC_pfieldcorr_arena_time(aLoop,:,:), 3);
    subplot(2,2,4); hold on
    plot(x,y, 'Color', aclr(aLoop,:))
    
end

y = nanmean(nanmean(HPC_pfieldcorr_room_time,3), 1);
subplot(2,2,1); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep -.2 .6])

y = nanmean(nanmean(HPC_pfieldcorr_arena_time,3), 1);
subplot(2,2,2); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep -.2 .6])

y = nanmean(nanmean(ACC_pfieldcorr_room_time,3), 1);
subplot(2,2,3); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep -.2 .6])

y = nanmean(nanmean(ACC_pfieldcorr_arena_time,3), 1);
subplot(2,2,4); hold on
plot(x,y, 'Color', 'k')
axis([-1 max_day_sep -.2 .6])
%%





















function [stability, stabprob, sfeppref] = sub_analyze(stri, segsi, strj, segsj)
    spksi = stri.svm_decoding.spks_bin(segsi, :);
    spksj = strj.svm_decoding.spks_bin(segsj, :);
    [cc1, pp1] = nancorr(spksi');
    [cc2, pp2] = nancorr(spksj');
%     sigs = pp1<=.05 & pp2 <=.05;
    mask = triu(ones(d1), 1)==1;
    ccdiff = cc1-cc2; ccdiff(~mask) = NaN;

    p1 = stri.pfields_smooth(segsi, :,:);
    p2 = strj.pfields_smooth(segsj, :,:);
    [d1, d2, d3] = size(p1);
    p1 = reshape(p1, [d1, d2*d3]);
    p2 = reshape(p2, [d1, d2*d3]);
    stability = NaN(d1); stabprob = NaN(d1);
    for ii = 1:d1
        [stability(ii), stabprob(ii)] = nancorr(p1(ii,:), p2(ii,:), 'Kendall');
    end
end





