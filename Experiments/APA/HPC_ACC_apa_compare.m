animals = {'HPCACC24500', 'HPCACC24502', 'HPCACC24503', 'Acc20832', 'Acc19947', 'Hipp18240'};% animals = {'Acc20832', 'Acc19947'};
% experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/';
experiment_folder = [];
experiment_folder{1} = 'D:/GarrettBlair/APA/';
experiment_folder{2} = 'D:/APA recordings/';


vars = {'pfield_decode.decode_dist', 'pfield_decode.decode_dist_shuffle_median',...
    'split_corr', 'pcell_stats.infoPerSpike'};

var_outname = {'BayDec_mean', 'BayDec_mean_rand', 'pfield_stability', 'pfield_info'};

top_struct = {'room', 'arena'};
analysis_version = 'v1.42';
this_ver = str2double(analysis_version(2:end));

time_origin = datetime(2020, 1, 1, 0, 0, 0); %

numAnimals = length(animals);


%% CAUSALITY
nan_templ = NaN(40, 2);
acc_GC_F = nan_templ;
acc_GC_pval = nan_templ;
acc_GC_shuffprob = nan_templ;

hpc_GC_F = nan_templ;
hpc_GC_pval = nan_templ;
hpc_GC_shuffprob = nan_templ;
run_granger = true;


acc_ac = nan_templ;
hpc_ac = nan_templ;
hpc_and_acc = nan_templ;
hpc_and_acc_rand = nan_templ;
hpc_ac_whenacc = nan_templ;
rand_ac_whenacc = nan_templ;
acc_ac_whenhpc = nan_templ;
rand_ac_whenhpc = nan_templ;


animal_idx = 0;

decode_mat = cell(numAnimals,1);
cind = 0;
for animalLoop = [1,3] % HPCACC24500 & HPCACC24503
    animal_name = animals{animalLoop};
    processed_dir = sprintf('%s%s/processed_files/', experiment_folder{1}, animal_name);
    fprintf('\n\nREPROCESS files for %s in folder  \n%s...\n', animal_name, processed_dir)
    temp = dir([processed_dir '*@placecells*']);
    nfiles = length(temp);
    this_idx = 0;
    animal_idx = animal_idx+1;
    decode_mat{animal_idx} = NaN(40, 4,4);
    for sessionLoop = 1:2:nfiles-1
        %%
        fname_acc = temp(sessionLoop).name;
        processedFile_acc   = [temp(sessionLoop).folder '/' temp(sessionLoop).name];
        processedFile_hpc   = [temp(sessionLoop+1).folder '/' temp(sessionLoop+1).name];
        hpc_acc_exist = isfile(processedFile_acc) & isfile(processedFile_hpc);
        prev_version    = load(processedFile_acc, 'analysis_version');
        prev_ver = str2double(prev_version.analysis_version(2:end));
        if prev_ver >= this_ver && hpc_acc_exist
            
            fprintf('~~~ACC file: %s...\n', temp(sessionLoop).name)
            fprintf('~~~HPC file: %s...\n', temp(sessionLoop+1).name)
            fprintf('\t\tanalysis_version : %1.2f\n', prev_ver)
            clearvars A H
            A = load(processedFile_acc, 'ms', 'params');
            H = load(processedFile_hpc, 'ms', 'params');
            [sessDate1, trialname1, trial_type1, trial_num1] = recdata_from_parentdir(A.ms.parentDir);
            [sessDate2, trialname2, trial_type2, trial_num2] = recdata_from_parentdir(H.ms.parentDir);
            if (sessDate1 - sessDate2) ~= 0
                warning('sesssion mismatch')
%                 break
            else
                a_ipos = nanmean( abs(A.ms.room.momentary_pos_info) - abs(A.ms.arena.momentary_pos_info), 1);
                h_ipos = nanmean( abs(H.ms.room.momentary_pos_info) - abs(H.ms.arena.momentary_pos_info), 1);
                %                 as = nanmean( A.ms.spks, 1);
                %                 hs = nanmean( H.ms.spks, 1);
                %                 a_ipos = nanmean( bin_spks_time(as, A.params.ipos_int_time, A.ms.timestamps./1000, false), 1);
                %                 h_ipos = nanmean( bin_spks_time(hs, H.params.ipos_int_time, H.ms.timestamps./1000, false), 1);
                if length(a_ipos)<length(h_ipos) % interp to the smaller vector
                    a_t = average_spks_time((A.ms.timestamps./1000)', A.params.ipos_int_time, A.ms.timestamps./1000, false, 'mean');
                    h_t = average_spks_time((H.ms.timestamps./1000)', H.params.ipos_int_time, H.ms.timestamps./1000, false, 'mean');
                    h_ipos = interp1(h_t, h_ipos, a_t, 'linear');
                elseif length(a_ipos)>length(h_ipos)
                    a_t = average_spks_time((A.ms.timestamps./1000)', A.params.ipos_int_time, A.ms.timestamps./1000, false, 'mean');
                    h_t = average_spks_time((H.ms.timestamps./1000)', H.params.ipos_int_time, H.ms.timestamps./1000, false, 'mean');
                    a_ipos = interp1(a_t, a_ipos, h_t, 'linear');
                end
                goods = ~isnan(h_ipos) & ~isnan(a_ipos);
                h_ipos = h_ipos(goods);
                a_ipos = a_ipos(goods);
                this_idx = this_idx+1;
                if run_granger == true
%                 [outstr, randstr] = granger_causality_Seth2015([a_ipos; h_ipos], .25, {'ACC' 'HPC'}, [], 10, false, false);
                [outstr, randstr] = APA_GrangerCausality_HPCACC_epochs(processedFile_hpc, processedFile_acc);
                acc_GC_F(this_idx, animal_idx) = outstr.F(2,1);%
                hpc_GC_F(this_idx, animal_idx) = outstr.F(1,2);%
                acc_GC_pval(this_idx, animal_idx) = outstr.pval(2,1);%
                hpc_GC_pval(this_idx, animal_idx) = outstr.pval(1,2);%
                acc_GC_shuffprob(this_idx, animal_idx) = outstr.prob(2,1);%
                hpc_GC_shuffprob(this_idx, animal_idx) = outstr.prob(1,2);%
                end
                %%                
                acc_ac(this_idx, animal_idx) = nanmean(A.ms.room.svm_decoding.var_binned==A.ms.room.svm_decoding.pred_real);
                hpc_ac(this_idx, animal_idx) = nanmean(H.ms.room.svm_decoding.var_binned==H.ms.room.svm_decoding.pred_real);
                
                q1 = A.ms.room.svm_decoding.pred_real;
                q2 = H.ms.room.svm_decoding.pred_real;
                if length(q1)<length(q2) % interp to the smaller vector
                    q2 = interp1(H.ms.room.svm_decoding.t, q2, A.ms.room.svm_decoding.t, 'nearest');
                elseif length(q1)>length(q2)
                    q1 = interp1(A.ms.room.svm_decoding.t, q1, H.ms.room.svm_decoding.t, 'nearest');
                end
                
                hpc_and_acc(this_idx, animal_idx) = nanmean(q1==q2);

                q1 = A.ms.room.svm_decoding.pred_rand;
                q2 = H.ms.room.svm_decoding.pred_rand;
                if length(q1)<length(q2) % interp to the smaller vector
                    q2 = interp1(H.ms.room.svm_decoding.t, q2, A.ms.room.svm_decoding.t, 'nearest');
                elseif length(q1)>length(q2)
                    q1 = interp1(A.ms.room.svm_decoding.t, q1, H.ms.room.svm_decoding.t, 'nearest');
                end
                hpc_and_acc_rand(this_idx, animal_idx) = nanmean(q1==q2);
                
                goods = A.ms.room.svm_decoding.var_binned==A.ms.room.svm_decoding.pred_real;
                hgoods = interp1(A.ms.room.svm_decoding.t,...
                    goods*1, H.ms.room.svm_decoding.t, 'nearest') == 1;
                hpc_ac_whenacc(this_idx, animal_idx)  = nanmean(H.ms.room.svm_decoding.var_binned(hgoods)==H.ms.room.svm_decoding.pred_real(hgoods));
                rand_ac_whenacc(this_idx, animal_idx) = nanmean(H.ms.room.svm_decoding.var_binned(hgoods)==H.ms.room.svm_decoding.pred_rand(hgoods));
                
                goods = H.ms.room.svm_decoding.var_binned==H.ms.room.svm_decoding.pred_real;
                agoods = interp1(H.ms.room.svm_decoding.t,...
                    goods*1, A.ms.room.svm_decoding.t, 'nearest') == 1;
                acc_ac_whenhpc(this_idx, animal_idx)  = nanmean(A.ms.room.svm_decoding.var_binned(agoods)==A.ms.room.svm_decoding.pred_real(agoods));
                rand_ac_whenhpc(this_idx, animal_idx) = nanmean(A.ms.room.svm_decoding.var_binned(agoods)==A.ms.room.svm_decoding.pred_rand(agoods));

                if false
                qq = {'A' 'H'};
%                 qqq = {'varopp_x', 'x', 'pred_x', 'rand_x'};
%                 yname = {'varopp_y', 'y', 'pred_y', 'rand_y'};
                qqq = [];%{'varopp', 'var_binned', 'pred_real', 'pred_rand'};
                q = NaN(length(qqq));
%                 A.ms.room.svm_decoding.varopp = H.ms.room.svm_decoding.var_binned;
                A.ms.room.svm_decoding.varopp = interp1(H.ms.room.svm_decoding.t,...
                    H.ms.room.svm_decoding.var_binned, A.ms.room.svm_decoding.t, 'nearest');
                A.ms.room.svm_decoding.varopp_x = interp1(H.ms.room.svm_decoding.t,...
                    H.ms.room.svm_decoding.x, A.ms.room.svm_decoding.t, 'linear');
                A.ms.room.svm_decoding.varopp_y = interp1(H.ms.room.svm_decoding.t,...
                    H.ms.room.svm_decoding.y, A.ms.room.svm_decoding.t, 'linear');
%                 H.ms.room.svm_decoding.varopp = A.ms.room.svm_decoding.var_binned;
                H.ms.room.svm_decoding.varopp = interp1(A.ms.room.svm_decoding.t,...
                    A.ms.room.svm_decoding.var_binned, H.ms.room.svm_decoding.t, 'nearest');
                H.ms.room.svm_decoding.varopp_x = interp1(A.ms.room.svm_decoding.t,...
                    A.ms.room.svm_decoding.x, H.ms.room.svm_decoding.t, 'linear');
                H.ms.room.svm_decoding.varopp_y = interp1(A.ms.room.svm_decoding.t,...
                    A.ms.room.svm_decoding.y, H.ms.room.svm_decoding.t, 'linear');
                            cind = cind+1;
                for i = 1:length(qqq)
                    for j = 1:length(qqq)
                        if contains(qqq{i}, 'x') || contains(qqq{j}, 'x')
%                             disp([i '_' j])
                            x1 = eval(sprintf('%s.ms.room.svm_decoding.%s;', qq{1}, qqq{i}));
                            y1 = eval(sprintf('%s.ms.room.svm_decoding.%s;', qq{1}, yname{i}));
                            x2 = eval(sprintf('%s.ms.room.svm_decoding.%s;', qq{2}, qqq{j}));
                            y2 = eval(sprintf('%s.ms.room.svm_decoding.%s;', qq{2}, yname{j}));
                            if length(x1)<length(x2) % interp to the smaller vector
                                x2 = interp1(H.ms.room.svm_decoding.t, x2, A.ms.room.svm_decoding.t, 'linear');
                                y2 = interp1(H.ms.room.svm_decoding.t, y2, A.ms.room.svm_decoding.t, 'linear');
                            elseif length(x1)>length(x2)
                                x1 = interp1(A.ms.room.svm_decoding.t, x1, H.ms.room.svm_decoding.t, 'linear');
                                y1 = interp1(A.ms.room.svm_decoding.t, y1, H.ms.room.svm_decoding.t, 'linear');
                            end
                            ddd = sqrt( (x1-x2).^2 + (y1-y2).^2 );
                            q(i,j) = nanmedian( ddd );
%                             hh = histcounts(ddd, [0:2:80], 'Normalization', 'probability');
                            figure(1235); subplot(4,4, sub2ind([4,4], i, j));
                            histogram(ddd, [0:2:80], 'Normalization', 'probability');
%                             plot(hh); 
                            ylim([0 .2]) 
                            drawnow
                        else % if  strcmp(qqq{i}, 'pred_real')
                            q1 = eval(sprintf('%s.ms.room.svm_decoding.%s;', qq{1}, qqq{i}));
                            q2 = eval(sprintf('%s.ms.room.svm_decoding.%s;', qq{2}, qqq{j}));
                            if length(q1)<length(q2) % interp to the smaller vector
                                q2 = interp1(H.ms.room.svm_decoding.t, q2, A.ms.room.svm_decoding.t, 'nearest');
                            elseif length(q1)>length(q2)
                                q1 = interp1(A.ms.room.svm_decoding.t, q1, H.ms.room.svm_decoding.t, 'nearest');
                            end
                            q(i,j) = nanmean(q1==q2);
                            ub = union(unique(q1), unique(q2));
                            ub = ub(~isnan(ub));
                            hh = histcounts2(q1,q2,ub, ub);
                            hh = normalize_matrix(hh);
                            figure(1234+cind); subplot(4,4, sub2ind([4,4], i, j));
                            imagesc(hh, [0 .2]); drawnow
                            
                            
                            %                             scatter(q1,q2, 'k.'); drawnow
                        end
                        
                    end
                end
                %%
                decode_mat{animal_idx}(this_idx,:,:) = q;
                end
                
            end
        else
            disp(hpc_acc_exist)
            
        end
    end
end

% a1= [decode_mat{1}(:,3,2); decode_mat{1}(:,2,3)]; 
% % b1= [decode_mat{1}(:,3,3)]; 
% c1= [decode_mat{1}(:,3,4); decode_mat{1}(:,4,3)]; 
% a2= [decode_mat{2}(:,3,2); decode_mat{2}(:,2,3)]; 
% % b2= [decode_mat{2}(:,3,3)]; 
% c2= [decode_mat{2}(:,3,4); decode_mat{2}(:,4,3)]; 


%%


% clear
a = load("D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_20_H13_11_53_TR15_@placecells_ACC_miniscope2.mat");
h = load("D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_20_H13_11_53_TR15_@placecells_HPC_miniscope1.mat");
% a = load("D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_20_H13_11_53_TR15_@placecells_ACC_miniscope2.mat");
% h = load("D:\GarrettBlair\APA\HPCACC24502\processed_files\2023_06_26_H11_39_11_TR18_@placecells_HPC_miniscope1.mat");

% aipos = nanmean( abs(a.ms.room.momentary_pos_info) - abs(a.ms.arena.momentary_pos_info), 1);
% hipos = nanmean( abs(h.ms.room.momentary_pos_info) - abs(h.ms.arena.momentary_pos_info), 1);
%
timeint = 1;

params = [];
params.pos_bins = [-45, -36:4:36, 45]; % in cm, x and y
params.occupancy_thresh = timeint; % min in seconds
params.skip_ensemble = true;
% aspks(isnan(aspks)) = 0;

if exist('a', 'var') == true
aspks = a.ms.spks;%ceil(4*normalize_rows(a.ms.neuron.S_matw));
    
    minshift = floor(length(aspks)*.1);
    maxshift = floor(length(aspks)*.9);
    shiftval = randi([minshift maxshift], size(aspks,1), 1);
    randdir = 1*(rand(size(aspks,1), 1)>=.5);
    randdir(randdir==0) = -1;
    shiftval = shiftval.*randdir;
    randaspks = NaN(size(aspks));
    % randaspks = circshift(aspks, shiftval(1), 2);
    for i = 1:size(aspks,1)
        randaspks(i,:) = circshift(aspks(i,:), shiftval(i));
    end
    
    a.ms.spks = aspks;
    [acc_roomipos,  ~, acc_ens_room]  = Fenton_ipos(a.ms, timeint, 'room',  params);
    [acc_arenaipos, ~, acc_ens_arena]  = Fenton_ipos(a.ms, timeint, 'arena', params);
    acc_ipos = abs(acc_roomipos) - abs(acc_arenaipos);
    
%     a.ms.spks = randaspks;
%     [rand_acc_roomipos,  ~, rand_acc_ens_room]     = Fenton_ipos(a.ms, timeint, 'room',  params);
%     [rand_acc_arenaipos, ~, rand_acc_ens_arena]  = Fenton_ipos(a.ms, timeint, 'arena', params);
%     a.ms.spks = aspks;
%     rand_acc_ipos = abs(rand_acc_roomipos) - abs(rand_acc_arenaipos);
end
%
hspks = h.ms.spks;

minshift = floor(length(hspks)*.1);
maxshift = floor(length(hspks)*.9);
shiftval = randi([minshift maxshift], size(hspks,1), 1);
randdir = 1*(rand(size(hspks,1), 1)>=.5);
randdir(randdir==0) = -1;
shiftval = shiftval.*randdir;
randhspks = NaN(size(hspks));
% randhspks = circshift(hspks, shiftval(1), 2);
for i = 1:size(hspks,1)
randhspks(i,:) = circshift(hspks(i,:), shiftval(i), 2);
end

h.ms.spks = hspks;
[hpc_roomipos,   room_ds, hpc_room_ensemb]   = Fenton_ipos(h.ms, timeint, 'room',  params);
[hpc_arenaipos, arena_ds, hpc_arena_ensemb]  = Fenton_ipos(h.ms, timeint, 'arena', params);
h_ipos = abs(hpc_roomipos) - abs(hpc_arenaipos);

% h.ms.spks = randhspks;
% [rand_hpc_roomipos,  ~, rand_hpc_ens_room]   = Fenton_ipos(h.ms, timeint, 'room',  params);
% [rand_hpc_arenaipos, ~, rand_hpc_ens_arena]  = Fenton_ipos(h.ms, timeint, 'arena', params);
% h.ms.spks = hspks;
% h_ipos_rand = abs(rand_hpc_roomipos) - abs(rand_hpc_arenaipos);

% [th, rho] = cart2pol(room_ds.x, room_ds.y);
% [th, rho] = cart2pol(h.ms.room.pfield_decode.xbin, h.ms.room.pfield_decode.ybin);
[th, rho] = cart2pol(h.ms.room.svm_decoding.x, h.ms.room.svm_decoding.y);
% [hth, ~] = cart2pol(h.ms.room.pfield_decode.xdecoded, h.ms.room.pfield_decode.ydecoded);
% [ath, ~] = cart2pol(a.ms.room.pfield_decode.xdecoded, a.ms.room.pfield_decode.ydecoded);
[hth, ~]  = cart2pol(h.ms.room.svm_decoding.pred_x, h.ms.room.svm_decoding.pred_y);
[rhth, ~] = cart2pol(h.ms.room.svm_decoding.rand_x, h.ms.room.svm_decoding.rand_y);
[ath, ~]  = cart2pol(a.ms.room.svm_decoding.pred_x, a.ms.room.svm_decoding.pred_y);
[rath, ~] = cart2pol(a.ms.room.svm_decoding.rand_x, a.ms.room.svm_decoding.rand_y);
% [th2, rho] = cart2pol(dx, dy);
[ang_dist, ~]   = angular_distance(th,  pi/2, 0);
[h_ang_dist, ~] = angular_distance(hth, pi/2, 0);
[a_ang_dist, ~] = angular_distance(ath, pi/2, 0);

% bd= sqrt( (a.ms.room.pfield_decode.xdecoded - h.ms.room.pfield_decode.xdecoded).^2 + (a.ms.room.pfield_decode.ydecoded - h.ms.room.pfield_decode.ydecoded).^2 );

figure(2); clf; hold on; 
plot(ang_dist, 'k.-'); 
plot(h_ang_dist, 'b.')
plot(a_ang_dist, 'm.'); 

her = nanmedian( sqrt( (h.ms.room.svm_decoding.pred_x - h.ms.room.svm_decoding.x).^2 + (h.ms.room.svm_decoding.pred_y - h.ms.room.svm_decoding.y).^2))
aer = nanmedian( sqrt( (a.ms.room.svm_decoding.pred_x - a.ms.room.svm_decoding.x).^2 + (a.ms.room.svm_decoding.pred_y - a.ms.room.svm_decoding.y).^2))
qq = nanmedian( sqrt( (h.ms.room.svm_decoding.pred_x - a.ms.room.svm_decoding.pred_x).^2 + (h.ms.room.svm_decoding.pred_y - a.ms.room.svm_decoding.pred_y).^2))
qq = nanmedian( sqrt( (h.ms.room.svm_decoding.pred_x - h.ms.room.svm_decoding.rand_x).^2 + (h.ms.room.svm_decoding.pred_y - h.ms.room.svm_decoding.rand_y).^2))
qq = nanmedian( sqrt( (a.ms.room.svm_decoding.pred_x - a.ms.room.svm_decoding.rand_x).^2 + (a.ms.room.svm_decoding.pred_y - a.ms.room.svm_decoding.rand_y).^2))

q = NaN(4);
qq = {'a' 'h'};
a.ms.room.svm_decoding.varopp = h.ms.room.svm_decoding.var_binned;
h.ms.room.svm_decoding.varopp = a.ms.room.svm_decoding.var_binned;
qqq = {'varopp', 'var_binned', 'pred_real', 'pred_rand'};
for i = 1:4
    for j = 1:4
        q1 = eval(sprintf('%s.ms.room.svm_decoding.%s;', qq{1}, qqq{i}));
        q2 = eval(sprintf('%s.ms.room.svm_decoding.%s;', qq{2}, qqq{j}));
        q(i,j) = nanmean(q1==q2);
    end
end
h_hr  = nanmean(h.ms.room.svm_decoding.pred_real == h.ms.room.svm_decoding.pred_rand);
a_ar  = nanmean(a.ms.room.svm_decoding.pred_real == a.ms.room.svm_decoding.pred_rand);
a_h   = nanmean(a.ms.room.svm_decoding.pred_real == h.ms.room.svm_decoding.pred_real);
ar_hr = nanmean(a.ms.room.svm_decoding.pred_rand == h.ms.room.svm_decoding.pred_rand);

%%
figure(6); clf
subplot(2,2,1)
hold on
plot(room_ds.t, nansum(acc_ipos,1), 'r');
plot(room_ds.t, nansum(rand_acc_ipos,1), 'k');
title('real ''r'' - rand ''k''')

subplot(2,2,3)
hold on
plot(room_ds.t, nansum(acc_ipos,1) - nansum(rand_acc_ipos,1));
title('real - rand')

subplot(2,2,2)
hold on
plot(room_ds.t, nansum(h_ipos,1), 'r');
plot(room_ds.t, nansum(h_ipos_rand,1), 'k');
title('real ''r'' - rand ''k''')

subplot(2,2,4)
hold on
plot(room_ds.t, nansum(h_ipos,1) - nansum(h_ipos_rand,1));

% figure(8);
ai = nanmean(acc_ipos,1);
hi = nanmean(h_ipos,1);
ai_c = nanmean(acc_ipos - rand_acc_ipos,1);
hi_c = nanmean(h_ipos - h_ipos_rand,1);
naninds = isnan(ai.*hi);
% (corr(ai(~naninds)', hi(~naninds)'))


combo = [h_ipos(:,~naninds); NaN(5, sum(~naninds)); acc_ipos(:,~naninds)];

%%

figure(10); clf; hold on; 
plot(h.ms.room.pfield_decode.time_average, h.ms.room.pfield_decode.decode_dist, 'r'); 
plot(h.ms.room.pfield_decode.time_average, h.ms.arena.pfield_decode.decode_dist, 'b');
% yyaxis('right')
% plot([0 room_ds.t(end)], [0 0], 'k')
% plot(room_ds.t, hi, 'm-', 'LineWidth', 2)
% plot(room_ds.t, hpc_roomipos, 'm-', 'LineWidth', 2)
% plot(room_ds.t, hpc_arenaipos*-1, 'm-', 'LineWidth', 2)

%%
figure(7); clf;
% subplot(3,1,1)
% % plot(room_ds.t, ang_dist, 'k');
% % ylim([-.5*pi max(ang_dist)*1.2])
% % yyaxis('right')
subplot(2,1,2)
hold on
plot(room_ds.t, ai-ai_c, 'm-');
plot(room_ds.t, hi-hi_c, 'b-');
ylim([min([ai,hi])*1.2 max([ai,hi])*5])
subplot(2,1,1)
hold on
plot(room_ds.t, ai, 'm-');
plot(room_ds.t, hi, 'b-');
ylim([min([ai,hi])*1.2 max([ai,hi])*5])
% plot(atemp.t, ai-hi);

% subplot(3,1,2)
% hold on
% plot(room_ds.t, ai_c);
% plot(room_ds.t, hi_c);
% subplot(3,1,3)
% hold on
% plot(room_ds.t, hi_c - ai_c);
% plot(room_ds.t, hi-ai);
%%
figure(12); clf;
% subplot(3,1,1)
plot(room_ds.t, ang_dist, 'k');
ylim([-.5*pi max(ang_dist)*1.2])
yyaxis('right')
hold on
plot(room_ds.t, a_ensemb, 'm-');
plot(room_ds.t, hi, 'b-');
ylim([min([ai,hi])*1.2 max([ai,hi])*5])
