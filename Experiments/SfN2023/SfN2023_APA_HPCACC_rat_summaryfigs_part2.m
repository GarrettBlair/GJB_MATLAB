%% CAUSALITY
addpath(genpath('C:\Users\gjb326\Documents\MATLAB\Barnett_and_Seth_MV_Granger_Causality'))
temp = NaN(40, 2);
acc_GC_F = temp;
acc_GC_pval = temp;
acc_GC_shuffprob = temp;

hpc_GC_F = temp;
hpc_GC_pval = temp;
hpc_GC_shuffprob = temp;

animal_idx = 0;
for animalLoop = [1,3] % HPCACC24500 & HPCACC24503
    animal_name = animals{animalLoop};
    processed_dir = sprintf('%s%s/processed_files/', experiment_folder{1}, animal_name);
    fprintf('\n\nREPROCESS files for %s in folder  \n%s...\n', animal_name, processed_dir)
    temp = dir([processed_dir '*@placecells*']);
    nfiles = length(temp);
    this_idx = 0;
    animal_idx = animal_idx+1;
    for sessionLoop = 1:2:nfiles-1
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
                break
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
                [outstr, randstr] = granger_causality_Seth2015([a_ipos; h_ipos], .25, {'ACC' 'HPC'}, [], 10, false, false);
                this_idx = this_idx+1;
                acc_GC_F(this_idx, animal_idx) = outstr.F(2,1);%
                hpc_GC_F(this_idx, animal_idx) = outstr.F(1,2);%
                acc_GC_pval(this_idx, animal_idx) = outstr.pval(2,1);%
                hpc_GC_pval(this_idx, animal_idx) = outstr.pval(1,2);%
                acc_GC_shuffprob(this_idx, animal_idx) = outstr.prob(2,1);%
                hpc_GC_shuffprob(this_idx, animal_idx) = outstr.prob(1,2);%
            end
        else
            disp(hpc_acc_exist)
            
        end
    end
end
rmpath(genpath('C:\Users\gjb326\Documents\MATLAB\Barnett_and_Seth_MV_Granger_Causality'))

% d = acc_GC_F-hpc_GC_F;
%%
a = acc_GC_F;
h = hpc_GC_F;
d = (a-h)./(a+h);
entr_min = numEntr(:,[1 3])./(sessTime(:,[1 3])./60000);
rec_day = exp_day(:,[1 3]);
rectype = exp_type(:,[1 3]);

figure(99); clf

for i = 1:2
    type = rectype(:,i);
    
    ne = entr_min(:,i);
%     valid = eval(valid_string); % ~isnan(ne) & ( strcmp(type, 'TR') );%| strcmp(type, 'CON') );%.*ya.*yh
    valid = ~isnan(ne) & ( strcmp(type, 'TR') );%| strcmp(type, 'CON') );%.*ya.*yh
    xs = 1:length(ne(valid));
    
    subplot(2,2,1); hold on
    scatter(entr_min(valid,i), a(valid,i), acc_scatter_style{:});
    scatter(entr_min(valid,i), h(valid,i), hpc_scatter_style{:});
    axis([-1 4 -.01 .05])
    
    subplot(2,2,3); hold on
    scatter(entr_min(valid,i), d(valid,i), 'ko');
    plot([-1 5], [0 0], 'k:')
    % axis([0 4 -.04 .04])
    axis([-1 4 -1.1 1.1])
    
    subplot(2,2,2); hold on
    plot(xs, a(valid,i), acc_plot_style{:});
    plot(xs, h(valid,i), hpc_plot_style{:});
    axis([-1 20 -.01 .05])
    
    subplot(2,2,4); hold on
    plot([-1 40], [0 0], 'k:')
    plot(xs, d(valid,i), 'k-');
    axis([-1 20 -1.1 1.1])
end

%%
figure(100); clf

for i = 1:2
    type = rectype(:,i);
    
    ne = entr_min(:,i);
    valid = ~isnan(ne) & ( strcmp(type, 'CON') );%| strcmp(type, 'CON') );%.*ya.*yh
%     valid = eval(valid_string); % ~isnan(ne) & ( strcmp(type, 'TR') );%| strcmp(type, 'CON') );%.*ya.*yh
    xs = 1:length(ne(valid));
    
    subplot(2,2,1); hold on
    scatter(entr_min(valid,i), a(valid,i), acc_scatter_style{:});
    scatter(entr_min(valid,i), h(valid,i), hpc_scatter_style{:});
    axis([-1 4 -.01 .05])
    
    subplot(2,2,3); hold on
    scatter(entr_min(valid,i), d(valid,i), 'ko');
    plot([-1 5], [0 0], 'k:')
    % axis([0 4 -.04 .04])
    axis([-1 4 -1.1 1.1])
    
    subplot(2,2,2); hold on
    plot(xs, a(valid,i), acc_plot_style{:});
    plot(xs, h(valid,i), hpc_plot_style{:});
    axis([-1 10 -.01 .05])
    
    subplot(2,2,4); hold on
    plot([-1 40], [0 0], 'k:')
    plot(xs, d(valid,i), 'k-');
    axis([-1 10 -1.1 1.1])
end