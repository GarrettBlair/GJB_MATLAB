% animals = {'HPCACC24500', 'HPCACC24502', 'HPCACC24503'};% animals = {'Acc20832', 'Acc19947'};
% experiment_folder = {'D:\GarrettBlair\APA\'};

animals = {'HPCACC24500', 'HPCACC24502', 'HPCACC24503', 'Acc20832', 'Acc19947', 'Hipp18240'};% animals = {'Acc20832', 'Acc19947'};
% experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/';
experiment_folder = [];
experiment_folder{1} = 'D:/GarrettBlair/APA/';
experiment_folder{2} = 'D:/APA recordings/';


vars = {'pfield_decode.decode_dist.*bin_w', 'pfield_decode.decode_dist_shuffle_median.*bin_w', 'svm_decoding.pred_err' , 'svm_decoding.rand_err_median',...
    'split_corr', 'pcell_stats.infoPerSpike', 'pcell_stats.infoProb<=.05', 'pcell_stats.splitcorrProb<=.05', 'split_corr', 'pcell_stats.infoPerSpike'};
var_outname = {'BayDec_mean', 'BayDec_mean_rand', 'SVMDec_mean', 'SVMDec_mean_rand', ...
    'pfield_stability', 'pfield_info_av', 'percent_place', 'percent_stable', 'pfield_stability', 'pfield_info'};

%                         HPC.room.BayDec_mean(thisind, animalLoop) = nanmedian(r.pfield_decode.decode_dist .*bin_w);
%                         HPC.room.BayDec_mean_rand(thisind, animalLoop) = nanmedian(r.pfield_decode.decode_dist_shuffle_median .*bin_w);
%                         HPC.room.SVMDec_mean(thisind, animalLoop) = nanmedian(r.svm_decoding.pred_err);
%                         HPC.room.SVMDec_mean_rand(thisind, animalLoop) = nanmedian(r.svm_decoding.rand_err_median);
%                         HPC.room.pfield_stability_av(thisind, animalLoop) = nanmedian(r.split_corr);
%                         HPC.room.pfield_info_av(thisind, animalLoop) = nanmedian(r.pcell_stats.infoPerSpike);
%                         HPC.room.percent_place(thisind, animalLoop) = nanmean(r.pcell_stats.infoProb<=.05);
%                         HPC.room.percent_stable(thisind, animalLoop) = nanmean(r.pcell_stats.splitcorrProb<=.05);
%                         HPC.room.pfield_stability{thisind, animalLoop} = r.split_corr;
%                         HPC.room.pfield_info{thisind, animalLoop} = r.pcell_stats.infoPerSpike;

top_struct = {'room', 'arena'};
analysis_version = 'v1.42';
this_ver = str2double(analysis_version(2:end));

time_origin = datetime(2020, 1, 1, 0, 0, 0); %

numAnimals = length(animals);

if true
    %%
    temp = NaN(40,numAnimals);
    HPC = [];
    ACC = [];
    HPC.room.BayDec_mean = temp;
    HPC.room.BayDec_mean_rand = temp;
    HPC.room.SVMDec_mean = temp;
    HPC.room.SVMDec_mean_rand = temp;
    HPC.room.pfield_stability_av = temp;
    HPC.room.pfield_info_av = temp;
    HPC.room.percent_place = temp;
    HPC.room.percent_stable = temp;
    
    temp = cell(40,numAnimals);
    
    HPC.room.pfield_stability = temp;
    HPC.room.pfield_info = temp;
    HPC.arena  = HPC.room;
    ACC.room   = HPC.room;
    ACC.arena  = HPC.room;
    
    
    % for i = 1:length(top_struct)
    %     for j = 1:length(vars)
    %                                         eval(sprintf('%s_%s = nanmedian(ms.%s.%s)', top_struct{i}, var_outname{j}, top_struct{i}, vars{j}))
    %         varstr = sprintf('%s_%s = NaN(100,4);', top_struct{i}, var_outname{j}, top_struct{i}, vars{j});
    %         eval(varstr);
    %     end
    % end
    exp_day = NaN(40,numAnimals);
    exp_type = cell(40,numAnimals);
    exp_num = NaN(40,numAnimals);
    rec_region = zeros(40,numAnimals); % 1==HPC, 2==ACC
    numEntr = NaN(40,numAnimals);
    sessTime = NaN(40,numAnimals);
    %
    for animalLoop = 1:numAnimals
        hpc_sesscount = 0;
        acc_sesscount = 0;
        for expLoop = 1:length(experiment_folder)
            animal_name = animals{animalLoop};
            processed_dir = sprintf('%s%s/processed_files/', experiment_folder{expLoop}, animal_name);
            fprintf('\n\nREPROCESS files for %s in folder  \n%s...\n', animal_name, processed_dir)
            temp = dir([processed_dir '*@placecells*']);
            nfiles = length(temp);
            for sessionLoop = 1:nfiles
                %%
                fname = temp(sessionLoop).name;
                processedFile   = [temp(sessionLoop).folder '/' fname];
                prev_version    = load(processedFile, 'analysis_version');
                %
                prev_ver = str2double(prev_version.analysis_version(2:end));
                if prev_ver >= this_ver
                    
                    fprintf('~~~file: %s...', fname)
                    fprintf(' analysis_version : %1.2f\n', prev_ver)
                    clearvars ms
                    load(processedFile, 'ms');
                    [sessDate, trialname, trial_type, trial_num] = recdata_from_parentdir(ms.parentDir);
                    
                    if ~isfield(ms, 'cameraName')
                        ms.cameraName = 'MiniLFOV';
                    end
                    if contains(ms.cameraName, 'ACC_miniscope')
                        rec_region(sessionLoop, animalLoop) = 2;
                        thisregion='ACC';
                        acc_sesscount = acc_sesscount+1;
                        thisind = acc_sesscount;
                    elseif contains(ms.cameraName, 'HPC_miniscope')
                        rec_region(sessionLoop, animalLoop) = 1;
                        thisregion='HPC';
                        hpc_sesscount = hpc_sesscount+1;
                        thisind = hpc_sesscount;
                    elseif contains(ms.cameraName, 'MiniLFOV') && contains(ms.parentDir, 'Acc')
                        rec_region(sessionLoop, animalLoop) = 2;
                        thisregion='ACC';
                        acc_sesscount = acc_sesscount+1;
                        thisind = acc_sesscount;
                    elseif contains(ms.cameraName, 'MiniLFOV') && contains(ms.parentDir, 'Hipp')
                        rec_region(sessionLoop, animalLoop) = 1;
                        thisregion='HPC';
                        hpc_sesscount = hpc_sesscount+1;
                        thisind = hpc_sesscount;
                    elseif contains(ms.cameraName, 'MiniLFOV') && contains(ms.parentDir, 'HPC')
                        rec_region(sessionLoop, animalLoop) = 1;
                        thisregion='HPC';
                        hpc_sesscount = hpc_sesscount+1;
                        thisind = hpc_sesscount;
                    end
                    
                    %                     a = ms.room.svm_decoding;
                    bin_w = median(abs(diff(ms.params.pos_bins)));
                    
                    exp_day(thisind, animalLoop) = (days(sessDate - time_origin));
                    exp_type{thisind, animalLoop} = trial_type;
                    exp_num(thisind, animalLoop) = trial_num;%thisind;%
                    numEntr(thisind, animalLoop) = length(ms.room.entranceTimes);
                    sessTime(thisind, animalLoop) = ms.timestamps(end) - ms.timestamps(1);
                    r = ms.room;
                    a = ms.arena;
                    
                    for vloop = 1:length(vars)
                        if iscell(eval(sprintf('%s.room.%s', thisregion, var_outname{vloop})))
                            str1 = sprintf('%s.room.%s{thisind, animalLoop}  = r.%s;', thisregion, var_outname{vloop}, vars{vloop});
                            str2 = sprintf('%s.arena.%s{thisind, animalLoop} = a.%s;', thisregion, var_outname{vloop}, vars{vloop});
                        else
                            str1 = sprintf('%s.room.%s(thisind, animalLoop)  = nanmean(r.%s);', thisregion, var_outname{vloop}, vars{vloop});
                            str2 = sprintf('%s.arena.%s(thisind, animalLoop) = nanmean(a.%s);', thisregion, var_outname{vloop}, vars{vloop});
                        end
%                         disp(str1)
                        eval(str1)
                        eval(str2)
                    end
                elseif any(strfind(temp(sessionLoop).date, '19-Oct-2023')) && prev_ver==1.2
                    fprintf('~~~PROCESSED file RERUN skipped: %s\n', fname)
                    fprintf('\tchanging analysis_version : %s\n', prev_version.analysis_version)
                    tempdata = load(processedFile);
                    analysis_version = 'v1.42'
                    save(processedFile,'analysis_version', '-append')
                else
                    fprintf('~~~PROCESSED file RERUN skipped: %s\n', fname)
                    fprintf('\tsame analysis_version : %s\n', prev_version.analysis_version)
                    %                     skipfiles = skipfiles + 1;
                end
                %
                %                 fprintf(' Done! \n\n')
            end % % RERUN PROCESSED FILES
        end
        exp_day(:,animalLoop) = exp_day(:,animalLoop) - min(exp_day(:,animalLoop));
    end
    
end
%%
cyan = [.2 .9  .8];
magenta = [.8 .5  .9];

hpc_scatter_style = {25, 'marker', 'o', 'MarkerFaceColor', cyan, 'MarkerEdgeColor', cyan./1.2, 'MarkerFaceAlpha', .6};
acc_scatter_style = {25, 'marker', 'd', 'MarkerFaceColor', magenta, 'MarkerEdgeColor', magenta./1.2, 'MarkerFaceAlpha', .6};
hpc_plot_style    = {'-', 'Color', cyan, 'LineWidth', 1};
acc_plot_style    = {'-', 'Color', magenta, 'LineWidth', 1};


ACC.room.allentr = [];
HPC.room.allentr = [];
ACC.room.allcorr = [];
HPC.room.allcorr = [];
ACC.room.allbits = [];
HPC.room.allbits = [];
ACC.room.allplace = [];
HPC.room.allplace = [];
ACC.arena.allentr = [];
HPC.arena.allentr = [];
ACC.arena.allcorr = [];
HPC.arena.allcorr = [];
ACC.arena.allbits = [];
HPC.arena.allbits = [];
ACC.arena.allplace = [];
HPC.arena.allplace = [];

all_a_entr  = [];
all_h_entr  = [];
all_a_time  = [];
all_h_time  = [];
all_a_bits  = [];
all_h_bits  = [];
all_a_corr  = [];
all_h_corr  = [];
all_a_place = [];
all_h_place = [];
valid_string = "~isnan(ne) & ( strcmp(type, 'TR')  | strcmp(type, 'CON') )";
% a_symbols = {'o' 'square' 'diamond' 'pentagram' 'hexagram' 'x'};
for strloop = 1:2
    if strloop ==1
        tempA = ACC.room;
        tempH = HPC.room;
        figure(91); clf
        set(gcf, 'Color', [1 .9 .9], 'Position', [120 550 800 400])
    elseif strloop == 2
        tempA = ACC.arena;
        tempH = HPC.arena;
        figure(92); clf
        set(gcf, 'Color', [.9 .9 1], 'Position', [950 550 800 400])
    end
    for animalLoop = 1:numAnimals
        %
        ne = numEntr(:,animalLoop)./(sessTime(:,animalLoop)./60000);
%         ne(ne<=1) = NaN;
        type = exp_type(:,animalLoop);
        valid = eval(valid_string) ; % ne<1 & ( strcmp(type, 'TR')  | strcmp(type, 'CON') );%.*ya.*yh
%         valid = strcmp(type, 'OF');%.*ya.*yh
        xv = ne(valid);
        [ya, xa] = gb_cell2vec(tempA.pfield_info(:,animalLoop), valid, xv);
        [yh, xh] = gb_cell2vec(tempH.pfield_info(:,animalLoop), valid, xv);
        all_a_entr = cat(1, all_a_entr, xa);
        all_h_entr = cat(1, all_h_entr, xh);
        %
        %     valid = (~isnan(ne));%.*ya.*yh
        subplot(2,4,1); hold on
        axis([-1 5 -.1 1.6])
        all_a_bits  = cat(1, all_a_bits, ya);
        all_h_bits  = cat(1, all_h_bits, yh);
        sub_err_bars(xa, ya, [ ], .05, 'quantile', acc_plot_style, acc_scatter_style)
        %         plot(xv, ya_av, acc_plot_style{:})
        set(gca, 'YTick', 0:.4:1.6, 'XTick', 0:4)
        xlabel('Entrances/Minute')
        ylabel('Bits/spike')
        
        %         subplot(2,4,2); hold on
        axis([-1 5 -.1 1.6])
        sub_err_bars(xh, yh, [ ], .05, 'quantile', hpc_plot_style, hpc_scatter_style)
        set(gca, 'YTick', 0:.4:1.6, 'XTick', 0:4)
        xlabel('Entrances/Minute')
        ylabel('Bits/spike')
        
        ya = tempA.percent_place(valid,animalLoop); xa = xv;
        yh = tempH.percent_place(valid,animalLoop); xh = xv;
        subplot(2,4,3); hold on
        axis([-1 5 -.1 1.1])
        all_a_place  = cat(1, all_a_place, ya);
        all_h_place  = cat(1, all_h_place, yh);
        scatter(xa, ya, acc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', 0:4)
        xlabel('Entrances/Minute')
        ylabel('% Place cells')
        
        %         subplot(2,4,4); hold on
        axis([-1 5 -.1 1.1])
        scatter(xh, yh, hpc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', 0:4)
        xlabel('Entrances/Minute')
        ylabel('% Place cells')
        
        [ya, xa] = gb_cell2vec(tempA.pfield_stability(:,animalLoop), valid, xv);
        [yh, xh] = gb_cell2vec(tempH.pfield_stability(:,animalLoop), valid, xv);
        subplot(2,4,5); hold on
        axis([-1 5 -.1 1.1])
        all_a_corr  = cat(1, all_a_corr, ya);
        all_h_corr  = cat(1, all_h_corr, yh);
        sub_err_bars(xa, ya, [ ], .05, 'quantile', acc_plot_style, acc_scatter_style)
        set(gca, 'YTick', 0:.2:1, 'XTick', 0:4)
        xlabel('Entrances/Minute')
        ylabel('Within Sess Corr')
        
        %         subplot(2,4,6); hold on
        axis([-1 5 -.1 1.1])
        sub_err_bars(xh, yh, [ ], .05, 'quantile', hpc_plot_style, hpc_scatter_style)
        set(gca, 'YTick', 0:.2:1, 'XTick', 0:4)
        xlabel('Entrances/Minute')
        ylabel('Within Sess Corr')
        
        ya = tempA.percent_stable(valid,animalLoop); xa = xv;
        yh = tempH.percent_stable(valid,animalLoop); xh = xv;
        subplot(2,4,7); hold on
        axis([-1 5 -.1 1.1])
        all_a_place  = cat(1, all_a_place, ya);
        all_h_place  = cat(1, all_h_place, yh);
        scatter(xa, ya, acc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', 0:4)
        xlabel('Entrances/Minute')
        ylabel('% Stable cells')
        
        %         subplot(2,4,8); hold on
        axis([-1 5 -.1 1.1])
        scatter(xh, yh, hpc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', 0:4)
        xlabel('Entrances/Minute')
        ylabel('% Stable cells')
    end
    if strloop ==1
        ACC.room.allentr = all_a_entr;
        HPC.room.allentr = all_h_entr;
        ACC.room.allcorr = all_a_corr;
        HPC.room.allcorr = all_h_corr;
        ACC.room.allbits = all_a_bits;
        HPC.room.allbits = all_h_bits;
        ACC.room.allplace = all_a_place;
        HPC.room.allplace = all_h_place;
    elseif strloop == 2
        ACC.arena.allentr = all_a_entr;
        HPC.arena.allentr = all_h_entr;
        ACC.arena.allcorr = all_a_corr;
        HPC.arena.allcorr = all_h_corr;
        ACC.arena.allbits = all_a_bits;
        HPC.arena.allbits = all_h_bits;
        ACC.arena.allplace = all_a_place;
        HPC.arena.allplace = all_h_place;
    end
    all_a_entr  = [];
    all_h_entr  = [];
    all_a_time  = [];
    all_h_time  = [];
    all_a_bits  = [];
    all_h_bits  = [];
    all_a_corr  = [];
    all_h_corr  = [];
    all_a_place = [];
    all_h_place = [];
end
%%
figure(1231); clf;
bins = [0:.025:1.0];
xs = [1:20:length(bins)];
h = histcounts2(HPC.room.allcorr, HPC.arena.allcorr, bins, bins, 'Normalization', 'count');
a = histcounts2(ACC.room.allcorr, ACC.arena.allcorr, bins, bins, 'Normalization', 'count');
subplot(2,2,1); hold on
p = pcolor(h); hold on; plot([1,length(bins)-1], [1,length(bins)-1], 'Color', cyan)
p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis([-5 length(bins)+5 -5 length(bins)+5])
axis square
colorbar
title('CA1')
ylabel('Pfield corr (Room)')
xlabel('Pfield corr (Arena)')
colormap viridis

subplot(2,2,2); hold on
p = pcolor(a); hold on; plot([1,length(bins)-1], [1,length(bins)-1], 'Color', magenta)
p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
title('ACC')
ylabel('Pfield corr (Room)')
xlabel('Pfield corr (Arena)')
axis square
colorbar
axis([-5 length(bins)+5 -5 length(bins)+5])

bins = [-1.0:.025:1.0];
subplot(2,1,2); hold on
hh = histogram(HPC.room.allcorr-HPC.arena.allcorr, bins, 'Normalization', 'count', 'FaceColor', cyan);
ha = histogram(ACC.room.allcorr-ACC.arena.allcorr, bins, 'Normalization', 'count', 'FaceColor', magenta);
plot([nanmedian(HPC.room.allcorr-HPC.arena.allcorr), nanmedian(HPC.room.allcorr-HPC.arena.allcorr)], [-100 800], '-', 'Color', cyan/1.5, 'LineWidth', 2)
plot([nanmedian(ACC.room.allcorr-ACC.arena.allcorr), nanmedian(ACC.room.allcorr-ACC.arena.allcorr)], [-100 800], '--', 'Color', magenta/1.5, 'LineWidth', 2)
axis([-1.1 1.1 -100 800])

% hh = histogram(HPC.room.allcorr-HPC.arena.allcorr, bins, 'Normalization', 'count', 'FaceColor', cyan);
% ha = histogram(ACC.room.allcorr-ACC.arena.allcorr, bins, 'Normalization', 'count', 'FaceColor', magenta);
% plot([nanmedian(HPC.room.allcorr-HPC.arena.allcorr), nanmedian(HPC.room.allcorr-HPC.arena.allcorr)], [0 700], 'b', 'LineWidth', 2)
% plot([nanmedian(ACC.room.allcorr-ACC.arena.allcorr), nanmedian(ACC.room.allcorr-ACC.arena.allcorr)], [0 700], 'm', 'LineWidth', 2)
legend([hh, ha], {'CA1' 'ACC'})
ylabel('Number of cells')
xlabel('Pfield corr (Room-Arena)')
% [c,p] = nancorr(all_h_bits, all_h_entr)
% [c,p] = nancorr(all_a_bits, all_a_entr)
% [c,p] = nancorr(all_h_bits, all_h_entr)
% [c,p] = nancorr(all_h_bits, all_h_entr)
% [h,p] = kstest(ACC.room.allcorr - ACC.arena.allcorr)
% [h,p] = kstest(HPC.room.allcorr - HPC.arena.allcorr)
%%
figure(1241); clf;
bins = [0:.025:1.0];
xs = [1:20:length(bins)];
h = histcounts2(HPC.room.allbits, HPC.arena.allbits, bins, bins, 'Normalization', 'count');
a = histcounts2(ACC.room.allbits, ACC.arena.allbits, bins, bins, 'Normalization', 'count');
subplot(2,2,1); hold on
p = pcolor(h); hold on; plot([1,length(bins)-1], [1,length(bins)-1], 'Color', cyan)
p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis([-5 length(bins)+5 -5 length(bins)+5])
axis square
colorbar
title('CA1')
ylabel('Pfield bits (Room)')
xlabel('Pfield bits (Arena)')
colormap plasma

subplot(2,2,2); hold on
p = pcolor(a); hold on; plot([1,length(bins)-1], [1,length(bins)-1], 'Color', magenta)
p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis square
colorbar
title('ACC')
ylabel('Pfield bits (Room)')
xlabel('Pfield bits (Arena)')
axis([-5 length(bins)+5 -5 length(bins)+5])

subplot(2,1,2); hold on
bins = [-1.0:.025:1.0];
hh = histogram(HPC.room.allbits-HPC.arena.allbits, bins, 'Normalization', 'count', 'FaceColor', cyan);
ha = histogram(ACC.room.allbits-ACC.arena.allbits, bins, 'Normalization', 'count', 'FaceColor', magenta);
plot([nanmedian(HPC.room.allbits-HPC.arena.allbits), nanmedian(HPC.room.allbits-HPC.arena.allbits)], [-100 1600], '-', 'Color', cyan/1.5, 'LineWidth', 2)
plot([nanmedian(ACC.room.allbits-ACC.arena.allbits), nanmedian(ACC.room.allbits-ACC.arena.allbits)], [-100 1600], '--', 'Color', magenta/1.5, 'LineWidth', 2)
axis([-1.1 1.1 -100 1500])
legend([hh, ha], {'CA1' 'ACC'})
ylabel('Number of cells')
xlabel('Pfield bits (Room-Arena)')
    
% [h,p] = kstest(ACC.room.allbits-ACC.arena.allbits)
% [h,p] = kstest(HPC.room.allbits-HPC.arena.allbits)
%%
for strloop = 1:2
    if strloop ==1
        tempA = ACC.room;
        tempH = HPC.room;
        figure(93); clf
        set(gcf, 'Color', [1 .9 .9], 'Position', [120 80 800 400])
    elseif strloop == 2
        tempA = ACC.arena;
        tempH = HPC.arena;
        figure(94); clf
        set(gcf, 'Color', [.9 .9 1], 'Position', [950 80 800 400])
    end
    for animalLoop = 1:numAnimals
        %%
        ne = round(exp_day(:,animalLoop));
        type = exp_type(:,animalLoop);
%         valid = ~isnan(ne) & ( strcmp(type, 'TR') );% | strcmp(type, 'CON') );%.*ya.*yh
        valid = eval(valid_string);
        xv = 1:length(ne(valid));
%         xv = round(ne(valid));
        ya_av = tempA.pfield_info_av(valid,animalLoop);
        yh_av = tempH.pfield_info_av(valid,animalLoop);
        [ya, xa] = gb_cell2vec(tempA.pfield_info(:,animalLoop), valid, xv);
        [yh, xh] = gb_cell2vec(tempH.pfield_info(:,animalLoop), valid, xv);
        all_a_time = cat(1, all_a_time, xa);
        all_h_time = cat(1, all_h_time, xh);
        
        subplot(2,4,1); hold on
        axis([-1 12 -.1 1.6])
        %         sub_err_bars(xa, ya, [], .2, 'quantile', acc_plot_style, {acc_scatter_style{:}, 'MarkerFaceColor', magenta./1.4, 'MarkerFaceAlpha', .6})
        sub_err_bars(xa, ya, [], .2, 'quantile', acc_plot_style, acc_scatter_style)
        plot(xv, ya_av, acc_plot_style{:})
        set(gca, 'YTick', 0:.4:1.6, 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('Bits/spike')
        
        %         subplot(2,4,2); hold on
        axis([-1 12 -.1 1.6])
        sub_err_bars(xh, yh, [], .2, 'quantile', hpc_plot_style, hpc_scatter_style)
        plot(xv, yh_av, hpc_plot_style{:})
        set(gca, 'YTick', 0:.4:1.6, 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('Bits/spike')
        
        ya_av = tempA.pfield_stability_av(valid,animalLoop);
        yh_av = tempH.pfield_stability_av(valid,animalLoop);
        [ya, xa] = gb_cell2vec(tempA.pfield_stability(:,animalLoop), valid, xv);
        [yh, xh] = gb_cell2vec(tempH.pfield_stability(:,animalLoop), valid, xv);
        subplot(2,4,5); hold on
        axis([-1 12 -.1 1.1])
        sub_err_bars(xa, ya, [], .2, 'quantile', acc_plot_style, acc_scatter_style)
        plot(xv, ya_av, acc_plot_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('Within Sess Corr')
        
        %         subplot(2,4,6); hold on
        axis([-1 12 -.1 1.1])
        sub_err_bars(xh, yh, [], .2, 'quantile', hpc_plot_style, hpc_scatter_style)
        plot(xv, yh_av, hpc_plot_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('Within Sess Corr')
        
        ya = tempA.percent_place(valid,animalLoop); xa = xv;
        yh = tempH.percent_place(valid,animalLoop); xh = xv;
        subplot(2,4,3); hold on
        axis([-1 12 -.1 1.1])
        plot(xa, ya, acc_plot_style{:})
        scatter(xa, ya, acc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('% Place cells')
        
        %         subplot(2,4,4); hold on
        axis([-1 12 -.1 1.1])
        plot(xa, yh, hpc_plot_style{:})
        scatter(xh, yh, hpc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('% Place cells')
        
        ya = tempA.percent_stable(valid,animalLoop); xa = xv;
        yh = tempH.percent_stable(valid,animalLoop); xh = xv;
        subplot(2,4,7); hold on
        axis([-1 12 -.1 1.1])
        plot(xa, ya, acc_plot_style{:})
        scatter(xa, ya, acc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        xlabel('Entrances/Minute')
        ylabel('% Stable cells')
        
        %         subplot(2,4,8); hold on
        axis([-1 12 -.1 1.1])
        plot(xa, yh, hpc_plot_style{:})
        scatter(xh, yh, hpc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        xlabel('Entrances/Minute')
        ylabel('% Stable cells')
    end
end
%%
figure(95); clf
set(gcf, 'Color', [1 .9 .9])
for strloop = 1:2
    if strloop ==1
        tempA = ACC.room;
        tempH = HPC.room;
        figure(95); clf
        set(gcf, 'Color', [1 .9 .9], 'Position', [350 380 550 400])
    elseif strloop == 2
        tempA = ACC.arena;
        tempH = HPC.arena;
        figure(96); clf
        set(gcf, 'Color', [.9 .9 1], 'Position', [950 380 550 400])
    end
    for animalLoop = 1:numAnimals
        %%
        %%%%%%%%%%%%%%% BY PERFORMANCE %%%%%%%%%%%%%%%%%%%%%%
        type = exp_type(:,animalLoop);
        ne = numEntr(:,animalLoop)./(sessTime(:,animalLoop)./60000);
%         valid = ~isnan(ne) & ( strcmp(type, 'TR') );%| strcmp(type, 'CON') );%.*ya.*yh
        valid = eval(valid_string);
        xs = ne(valid);
        
        ya = tempA.BayDec_mean(:,animalLoop);
        yh = tempH.BayDec_mean(:,animalLoop);
        yar = tempA.BayDec_mean_rand(:,animalLoop);
        yhr = tempH.BayDec_mean_rand(:,animalLoop);
        ya = 100*ya./yar;
        yh = 100*yh./yhr;
%         ya = 100*(yar-ya)./yar;
%         yh = 100*(yhr-yh)./yhr;
        subplot(2,2,1); hold on
        plot([-2 6], [100 100], 'k:')
        axis([-1 5 0 120])
        scatter(xs, yh(valid), hpc_scatter_style{:})
        scatter(xs, ya(valid), acc_scatter_style{:})
        set(gca, 'YTick', [0:25:100], 'XTick', 0:4)
        xlabel('Entrances/Minute')
        ylabel(sprintf('Bayesian performance\n^{real}/_{random}'))
        
        ya = tempA.SVMDec_mean(:,animalLoop);
        yh = tempH.SVMDec_mean(:,animalLoop);
        yar = tempA.SVMDec_mean_rand(:,animalLoop);
        yhr = tempH.SVMDec_mean_rand(:,animalLoop);
        ya = 100*ya./yar;
        yh = 100*yh./yhr;
        subplot(2,2,2); hold on
        plot([-2 6], [100 100], 'k:')
        axis([-1 5 0 120])
        scatter(xs, yh(valid), hpc_scatter_style{:})
        scatter(xs, ya(valid), acc_scatter_style{:})
        set(gca, 'YTick', [0:25:100], 'XTick', 0:4)
        xlabel('Entrances/Minute')
        ylabel('SVM performance')
        
        %%%%%%%%%%%%%%% BY TIME %%%%%%%%%%%%%%%%%%%%%%
        ne = round(exp_day(:,animalLoop));
        valid = ~isnan(ne) & ( strcmp(type, 'TR') );% | strcmp(type, 'CON') );%.*ya.*yh
        xs = 1:length(ne(valid));
        %     xs = ne(valid);%linspace(1, nanmax(ne), length(ne(valid)));
        
        ya = tempA.BayDec_mean(:,animalLoop);
        yh = tempH.BayDec_mean(:,animalLoop);
        yar = tempA.BayDec_mean_rand(:,animalLoop);
        yhr = tempH.BayDec_mean_rand(:,animalLoop);
        ya = 100*ya./yar;
        yh = 100*yh./yhr;
        subplot(2,2,3); hold on
        plot([-2 40], [100 100], 'k:')
        axis([-1 12 0 120])
        plot(xs, ya(valid), acc_plot_style{:})
        plot(xs, yh(valid), hpc_plot_style{:})
        set(gca, 'YTick', [0:25:100], 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel(sprintf('Bayesian performance\n^{real}/_{random}'))
        
        ya = tempA.SVMDec_mean(:,animalLoop);
        yh = tempH.SVMDec_mean(:,animalLoop);
        yar = tempA.SVMDec_mean_rand(:,animalLoop);
        yhr = tempH.SVMDec_mean_rand(:,animalLoop);
        ya = 100*ya./yar;
        yh = 100*yh./yhr;
        subplot(2,2,4); hold on
        plot([-2 40], [100 100], 'k:')
        axis([-1 12 0 120])
        plot(xs, ya(valid), acc_plot_style{:})
        plot(xs, yh(valid), hpc_plot_style{:})
        set(gca, 'YTick', [0:25:100], 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('SVM performance')
    end
end



%% CAUSALITY
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
    axis([-1 10 -.01 .05])
    
    subplot(2,2,4); hold on
    plot([-1 40], [0 0], 'k:')
    plot(xs, d(valid,i), 'k-');
    axis([-1 10 -1.1 1.1])
end

%%
figure(100); clf

for i = 1:2
    type = rectype(:,i);
    
    ne = entr_min(:,i);
    valid = ~isnan(ne) & ( strcmp(type, 'CON') );%| strcmp(type, 'CON') );%.*ya.*yh
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

%%









function [yys, xxs] = gb_cell2vec(data_str, valid, xs)
yys = [];
xxs = [];
xi = find(valid);
for i = 1:length(xi)%size(data_str,1)
    yys = cat(1, yys, data_str{xi(i)});
    xxs = cat(1, xxs, xs(i) + 0*data_str{xi(i)});
end
% xxs = xxs(~isnan(xxs));
% yys = yys(~isnan(yys));
end
function sub_err_bars(x, y, q_vals, barwidth, barmethod, plotstyle_arg, scatterstyle_arg)
s = barwidth;
xa = x;
ya = y;
xv = unique(xa);
for i = 1:length(xv)
    if strcmp(barmethod, 'quantile')
        q  = quantile(ya(xa==xv(i)), .5);
        if ~isempty(q_vals)
            q1 = quantile(ya(xa==xv(i)), q_vals(1));
            q2 = quantile(ya(xa==xv(i)), q_vals(2));
        end
    elseif strcmp(barmethod, 'std')
        q  = nanmean(ya(xa==xv(i)));
        if ~isempty(q_vals)
            q1 = q - q_vals(1)*nanstd(ya(xa==xv(i)));
            q2 = q + q_vals(2)*nanstd(ya(xa==xv(i)));
        end
    end
    hold on
    if ~isempty(q_vals)
        plot([xv(i), xv(i)], [q1, q2], plotstyle_arg{:});
        plot([xv(i)-s, xv(i)+s], [q1, q1], plotstyle_arg{:});
        plot([xv(i)-s, xv(i)+s], [q2, q2], plotstyle_arg{:});
    end
    scatter(xv(i), q, scatterstyle_arg{:});
end
end



