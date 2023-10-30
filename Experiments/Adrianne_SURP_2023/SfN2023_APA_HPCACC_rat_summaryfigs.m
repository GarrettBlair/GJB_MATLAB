% animals = {'HPCACC24500', 'HPCACC24502', 'HPCACC24503'};% animals = {'Acc20832', 'Acc19947'};
% experiment_folder = {'D:\GarrettBlair\APA\'};

animals = {'HPCACC24500', 'HPCACC24502', 'HPCACC24503', 'Acc20832', 'Acc19947', 'Hipp18240'};% animals = {'Acc20832', 'Acc19947'};
% ~! HPCACC24502 is very arena-preferential (pfield bits and stability)
% animals = {'HPCACC24500', 'HPCACC24503', 'Acc20832', 'Acc19947', 'Hipp18240'};% animals = {'Acc20832', 'Acc19947'};
% experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/';
experiment_folder = [];
experiment_folder{1} = 'D:/GarrettBlair/APA/';
% experiment_folder{2} = 'D:/APA recordings/';


vars = {'pfield_decode.decode_dist.*bin_w', 'pfield_decode.decode_dist_shuffle_median.*bin_w', 'svm_decoding.pred_err' , 'svm_decoding.rand_err_median',...
    'split_corr', 'pcell_stats.infoPerSpike', 'pcell_stats.coherence', 'pcell_stats.infoProb<=.05', 'pcell_stats.splitcorrProb<=.05', 'split_corr',...
    'pcell_stats.infoPerSpike', 'pcell_stats.coherence', 'ipos_roompref'};
var_outname = {'BayDec_mean', 'BayDec_mean_rand', 'SVMDec_mean', 'SVMDec_mean_rand', ...
    'pfield_stability_av', 'pfield_info_av', 'pfield_coh_av', 'percent_place', 'percent_stable', 'pfield_stability',...
    'pfield_info', 'pfield_coherence', 'ipos_roompref'};

top_struct = {'room', 'arena'};
analysis_version = 'v1.42';
this_ver = str2double(analysis_version(2:end));

time_origin = datetime(2020, 1, 1, 0, 0, 0); %

numAnimals = length(animals);

if true
    %%
    temp = NaN(40,numAnimals);
    HPC = [];
    HPC.room.BayDec_mean = temp;
    HPC.room.BayDec_mean_rand = temp;
    HPC.room.SVMDec_mean = temp;
    HPC.room.SVMDec_mean_rand = temp;
    HPC.room.pfield_stability_av = temp;
    HPC.room.pfield_info_av = temp;
    HPC.room.pfield_coh_av = temp;
    HPC.room.percent_place = temp;
    HPC.room.percent_stable = temp;
    HPC.room.ipos_roompref = temp;
    HPC.ncells = temp;
    
    temp = cell(40,numAnimals);
    
    HPC.room.pfield_stability = temp;
    HPC.room.pfield_info = temp;
    HPC.room.pfield_coherence = temp;
    HPC.arena  = HPC.room;
    ACC  = HPC;
    
    
    exp_day = NaN(40,numAnimals);
    exp_type = cell(40,numAnimals);
    exp_num = NaN(40,numAnimals);
    rec_region = zeros(40,numAnimals); % 1==HPC, 2==ACC
    numEntr = NaN(40,numAnimals);
    sessTime = NaN(40,numAnimals);
    prop_inside = NaN(40,numAnimals);
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
                    prevparams = load(processedFile, 'params');
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
                    [theta, rho] = cart2pol(r.x, r.y);
%                     prop_inside = sum(theta>pi/3 & theta<2*pi/3)/length(theta);
                    prop_inside(thisind, animalLoop) = sum(theta>=pi/3 & theta<=2*pi/3)/length(theta);
                    
                    r.pcell_stats.coherence = pfield_coherence_calc(r.pfields, r.vmap);
                    a.pcell_stats.coherence = pfield_coherence_calc(a.pfields, a.vmap);
%                     concordant_pfield_change(ms, prevparams.params);
                    ipos = abs(r.momentary_pos_info) - abs(a.momentary_pos_info);
                    ipos_mean= nanmean(ipos,1);
                    r.ipos_roompref = nansum(ipos_mean>0)./sum(~isnan(ipos_mean));
                    a.ipos_roompref = nansum(ipos_mean<0)./sum(~isnan(ipos_mean));
%                     goodinds = any(~isnan(ipos),1);
%                     ipos = zscore( ipos(:, goodinds) );
%                     c = NaN(size(ipos,1), 1); p=c;
%                     for i =1:size(ipos,1); [c(i), p(i)] = circ_corrcl(theta(goodinds), ipos(i,:)); end
                    if strcmp(thisregion, 'HPC')
                        HPC.ncells(thisind, animalLoop) = size(ms.spks,1);
                    elseif strcmp(thisregion, 'ACC')
                        ACC.ncells(thisind, animalLoop) = size(ms.spks,1);
                    end
                    for vloop = 1:length(vars)
                        if iscell(eval(sprintf('%s.room.%s', thisregion, var_outname{vloop})))
                            str1 = sprintf('%s.room.%s{thisind, animalLoop}  = r.%s;', thisregion, var_outname{vloop}, vars{vloop});
                            str2 = sprintf('%s.arena.%s{thisind, animalLoop} = a.%s;', thisregion, var_outname{vloop}, vars{vloop});
                        else
                            str1 = sprintf('%s.room.%s(thisind, animalLoop)  = nanmean(r.%s);', thisregion, var_outname{vloop}, vars{vloop});
                            str2 = sprintf('%s.arena.%s(thisind, animalLoop) = nanmean(a.%s);', thisregion, var_outname{vloop}, vars{vloop});
                        end
                        eval(str1)
                        eval(str2)
                    end
                else
                    fprintf('~~~PROCESSED file RERUN skipped: %s\n', fname)
                    fprintf('\tsame analysis_version : %s\n', prev_version.analysis_version)
                end
            end % % RERUN PROCESSED FILES
        end
        exp_day(:,animalLoop) = exp_day(:,animalLoop) - min(exp_day(:,animalLoop));
    end
    
end

%%
cyan = [.2 .9  .8];
magenta = [.8 .5  .9];

hpc_scatter_style = {40, 'marker', 'o', 'MarkerFaceColor', cyan, 'MarkerEdgeColor', cyan./1.2, 'MarkerFaceAlpha', .6};
acc_scatter_style = {40, 'marker', 'o', 'MarkerFaceColor', magenta, 'MarkerEdgeColor', magenta./1.2, 'MarkerFaceAlpha', .6};
hpc_plot_style    = {'-', 'Color', cyan, 'LineWidth', 1};
acc_plot_style    = {'-', 'Color', magenta, 'LineWidth', 1};

entrperMin = numEntr./(sessTime./60000);
isTR = strcmp(exp_type, 'TR');
isCON = strcmp(exp_type, 'CON');
valid = ~isnan(entrperMin) & (isTR|isCON);

HPC.room.BaysPerf   = 1 - HPC.room.BayDec_mean./HPC.room.BayDec_mean_rand;
HPC.arena.BaysPerf  = 1 - HPC.arena.BayDec_mean./HPC.arena.BayDec_mean_rand;
HPC.room.SVMPerf    = 1 - HPC.room.SVMDec_mean./HPC.room.SVMDec_mean_rand;
HPC.arena.SVMPerf   = 1 - HPC.arena.SVMDec_mean./HPC.arena.SVMDec_mean_rand;
ACC.room.BaysPerf   = 1 - ACC.room.BayDec_mean./ACC.room.BayDec_mean_rand;
ACC.arena.BaysPerf  = 1 - ACC.arena.BayDec_mean./ACC.arena.BayDec_mean_rand;
ACC.room.SVMPerf    = 1 - ACC.room.SVMDec_mean./ACC.room.SVMDec_mean_rand;
ACC.arena.SVMPerf   = 1 - ACC.arena.SVMDec_mean./ACC.arena.SVMDec_mean_rand;

var_outname = {'BaysPerf', 'SVMPerf', 'pfield_stability_av', 'pfield_info_av',...
    'percent_place', 'percent_stable', 'pfield_coh_av', 'ipos_roompref'};
for rrr = 1:nvars
    eval(sprintf('anova_%s = [];', var_outname{rrr}))
    eval(sprintf('anova_%s_label = [];', var_outname{rrr}))
    eval(sprintf('subj_anova_%s = [];', var_outname{rrr}))
    eval(sprintf('subj_anova_%s_label = [];', var_outname{rrr}))
end
nvars = length(var_outname);
reg = {'HPC', 'ACC'};
framecolor = {'r', 'b'};
maxys = ones(nvars,1) + .1;
maxys(1:2) = [.8 .8];
maxys(4) = [1.3];
maxys(7) = [.65];
minys = zeros(nvars,1) - .1;
% top_struct = {'room', 'arena'};
figure(9854); clf
figure(9855); clf

anova_label_num = 0;
% 1 = HPC_room_var
% 2 = HPC_arena_var
% 3 = ACC_room_var
% 4 = ACC_arena_var
for r = 1:length(reg)
    figure(9845+r); clf
    if r==1; regcolor = cyan; style = hpc_scatter_style;  elseif r==2; regcolor = magenta;  style = acc_scatter_style; end
    reg_offset = ( r./(1+length(reg)) - .5 ) ;
    for rr = 1:length(top_struct)
        frame_offset = ( rr./(1+length(top_struct)) - .5 )/2 ;
    anova_label_num = anova_label_num+1;
        for rrr = 1:nvars
                eval(sprintf('temp = %s.%s.%s;', reg{r}, top_struct{rr}, var_outname{rrr}))
                eval(sprintf('temp_anova = anova_%s;', var_outname{rrr}));
                eval(sprintf('temp_anova_label = anova_%s_label;', var_outname{rrr}));
                eval(sprintf('subj_temp_anova = subj_anova_%s;', var_outname{rrr}));
                eval(sprintf('subj_temp_anova_label = subj_anova_%s_label;', var_outname{rrr}));
                
                tempv = temp;
                tempv(~valid) = NaN;
                allsess_vals = tempv(~isnan(tempv));
                temp_anova = cat(1, temp_anova, allsess_vals);
                temp_anova_label = cat(1, temp_anova_label, allsess_vals*0 + anova_label_num);
                eval(sprintf('anova_%s = temp_anova;', var_outname{rrr}));
                eval(sprintf('anova_%s_label = temp_anova_label;', var_outname{rrr}));

                all_av = nanmean(tempv(:));
                all_var = nanstd(tempv(:));
                ncomps = sum(~isnan(tempv(:)));
                if rrr==1
                    eval(sprintf('%s.%s.ncomps = ncomps;', reg{r}, top_struct{rr}));
                end
                
                av = nanmean(tempv,1);
                vars = nanstd(tempv);
                xs = find(~isnan(av)) + frame_offset;
                av = av(~isnan(av));
                vars = vars(~isnan(vars));
                subj_temp_anova = cat(1, subj_temp_anova, av');
                subj_temp_anova_label = cat(1, subj_temp_anova_label, av'*0 + anova_label_num);
                eval(sprintf('subj_anova_%s = subj_temp_anova;', var_outname{rrr}));
                eval(sprintf('subj_anova_%s_label = subj_temp_anova_label;', var_outname{rrr}));
                
                ind = sub2ind([ceil(nvars/2),2],rrr);
                figure(9845+r);
                subplot(ceil(nvars/2), 2, ind); hold on
                for i = 1:length(xs)
                    hold on
                    errorbar(xs(i), av(i), vars(i), 'Color', 'k')
                end
                scatter(xs, av, 20, 'o', 'MarkerEdgeColor', framecolor{rr}, 'MarkerFaceColor', framecolor{rr})
                axis([.5 6.5 minys(rrr), maxys(rrr)])
                varname = sprintf('%s.%s', reg{r}, var_outname{rrr});
                title(varname, 'Interpreter', 'none', 'Color', regcolor/1.5)
                set(gca, 'XTick', [1:6], 'Color', 'none')
                
                figure(9855);
                subplot(ceil(nvars/2), 2, ind); hold on
                xs = zeros(length(xs),1)+reg_offset + rr;
                bar(xs(1), mean(av), 'FaceColor', regcolor/1.5, 'BarWidth', .25)
                scatterxs = gb_rand_jitter(av, 15);
                scatter(scatterxs+xs+.05, av, style{:}, 'MarkerEdgeColor', 'k')
                errorbar(xs(1), mean(av), std(av), 'k', 'LineWidth', 2)
                axis([.5 2.5 minys(rrr), maxys(rrr)])
                set(gca, 'XTick', [1 2], 'XTickLabel', {'room' 'arena'}, 'Color', 'none')
                varname = sprintf('%s', var_outname{rrr});
                title(varname, 'Interpreter', 'none', 'Color', 'k')
                
                figure(9854);
                subplot(ceil(nvars/2), 2, ind); hold on
                xs = zeros(ncomps,1)+reg_offset + rr;
                bar(xs(1), all_av, 'FaceColor', regcolor/1.5, 'BarWidth', .25)
                scatterxs = gb_rand_jitter(tempv(~isnan(tempv)), 15);
                scatter(scatterxs+xs + .05, tempv(~isnan(tempv)), style{:});%, 'MarkerEdgeColor', 'k')
                errorbar(xs(1), all_av, all_var, 'k', 'LineWidth', 2)
                axis([.5 2.5 minys(rrr), maxys(rrr)])
                set(gca, 'XTick', [1 2], 'XTickLabel', {'room' 'arena'}, 'Color', 'none')
                varname = sprintf('%s', var_outname{rrr});
                title(varname, 'Interpreter', 'none', 'Color', 'k')
        end
    end
end
%%
ncellsubjorder = [6 2 1 3 4 5];
c = HPC.ncells;
c(~valid) = NaN;
nc = nanmean(c,1);
vc = nanstd(c,1);
figure(114); clf; hold on
bar([1:numAnimals]-.15, nc(ncellsubjorder), 'FaceColor', cyan, 'BarWidth', .2)
for i = 1:length(nc)
    errorbar(i-.15, nc(ncellsubjorder(i)), vc(ncellsubjorder(i)), 'Color', 'k')
end
c = ACC.ncells;
c(~valid) = NaN;
nc = nanmean(c,1);
vc = nanstd(c,1);
bar([1:numAnimals]+.15, nc(ncellsubjorder), 'FaceColor', magenta, 'BarWidth', .2)
for i = 1:length(nc)
    errorbar(i+.15, nc(ncellsubjorder(i)), vc(ncellsubjorder(i)), 'Color', 'k')
end
axis([.15 numAnimals+.85 0 900])
an = animals(ncellsubjorder);
set(gca, 'YTick', [0:100:900], 'XTick', 1:numAnimals, 'XTickLabel', an(:), 'XTickLabelRotation', -90);

nsessvald = sum(valid,1);
for i = 1:length(nsessvald)
    text(i-.25, -20, sprintf('%d sess', nsessvald(ncellsubjorder(i))), 'Rotation', -90, 'FontSize', 8)
end
%%
label = anova_pfield_info_av_label;
ishpc = label==1 | label==2;
isroom = label==1 | label==3;
effects = [];

% [p,tbl,stats] = kruskalwallis(anova_pfield_info_av, anova_pfield_info_av_label)
% c = multcompare(stats)
effects.region.spatial_info = ranksum(anova_pfield_info_av(ishpc), anova_pfield_info_av(~ishpc));
effects.frame.spatial_info = ranksum(anova_pfield_info_av(isroom), anova_pfield_info_av(~isroom));

% [p,tbl,stats] = kruskalwallis(anova_percent_place, anova_percent_place_label)
% c = multcompare(stats)
effects.region.perc_place = ranksum(anova_percent_place(ishpc), anova_percent_place(~ishpc));
effects.frame.perc_place = ranksum(anova_percent_place(isroom), anova_percent_place(~isroom));



% [p,tbl,stats] = kruskalwallis(anova_pfield_stability_av, anova_pfield_stability_av_label)
% c = multcompare(stats)
effects.region.field_corr = ranksum(anova_pfield_stability_av(ishpc), anova_pfield_stability_av(~ishpc));
effects.frame.field_corr = ranksum(anova_pfield_stability_av(isroom), anova_pfield_stability_av(~isroom));


% [p,tbl,stats] = kruskalwallis(anova_percent_stable, anova_percent_stable_label)
% c = multcompare(stats)
effects.region.perc_stable = ranksum(anova_percent_stable(ishpc), anova_percent_stable(~ishpc));
effects.frame.perc_stable = ranksum(anova_percent_stable(isroom), anova_percent_stable(~isroom));


[p,tbl,stats] = kruskalwallis(anova_BaysPerf, anova_BaysPerf_label);%, 'off');
c = multcompare(stats)
effects.region.Bayes = ranksum(anova_BaysPerf(ishpc), anova_BaysPerf(~ishpc));
effects.frame.Bayes = ranksum(anova_BaysPerf(isroom), anova_BaysPerf(~isroom));


% [p,tbl,stats] = kruskalwallis(anova_SVMPerf, anova_SVMPerf_label)
% c = multcompare(stats)
effects.region.SVM = ranksum(anova_SVMPerf(ishpc), anova_SVMPerf(~ishpc));
effects.frame.SVM = ranksum(anova_SVMPerf(isroom), anova_SVMPerf(~isroom));

% [p,tbl,stats] = kruskalwallis(anova_pfield_coh_av, anova_pfield_coh_av_label)
% c = multcompare(stats)
effects.region.field_coh = ranksum(anova_pfield_coh_av(ishpc), anova_pfield_coh_av(~ishpc));
effects.frame.field_coh = ranksum(anova_pfield_coh_av(isroom), anova_pfield_coh_av(~isroom));
effects.region_frame.field_coh = ranksum(anova_pfield_coh_av(ishpc&~isroom), anova_pfield_coh_av(~ishpc&~isroom));


%%
% cyan = [.2 .9  .8];
% magenta = [.8 .5  .9];
% 
hpc_scatter_style = {25, 'marker', 'o', 'MarkerFaceColor', cyan, 'MarkerEdgeColor', cyan./1.2, 'MarkerFaceAlpha', .6};
acc_scatter_style = {25, 'marker', 'd', 'MarkerFaceColor', magenta, 'MarkerEdgeColor', magenta./1.2, 'MarkerFaceAlpha', .6};
hpc_plot_style    = {'-', 'Color', cyan, 'LineWidth', 1};
acc_plot_style    = {'-', 'Color', magenta, 'LineWidth', 1};


ACC.room.allentr = [];
HPC.room.allentr = [];
ACC.room.allcorr = [];
HPC.room.allcorr = [];
ACC.room.allcoh  = [];
HPC.room.allcoh  = [];
ACC.room.allbits = [];
HPC.room.allbits = [];
ACC.room.allplace = [];
HPC.room.allplace = [];
ACC.arena.allentr = [];
HPC.arena.allentr = [];
ACC.arena.allcorr = [];
HPC.arena.allcorr = [];
ACC.arena.allcoh  = [];
HPC.arena.allcoh  = [];
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
all_a_coh   = [];
all_h_coh   = [];
all_a_place = [];
all_h_place = [];
valid_string = "~isnan(ne) & ( strcmp(type, 'TR')  | strcmp(type, 'CON') )";
% valid_string = "ne<1 & ( strcmp(type, 'TR')  | strcmp(type, 'CON') )";
% a_symbols = {'o' 'square' 'diamond' 'pentagram' 'hexagram' 'x'};
xticks = 0:1:4;
maxx = max(xticks) + median(abs(diff(xticks)));
minx = min(xticks) - median(abs(diff(xticks)));
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
%         ne = prop_inside(:,animalLoop);
%         valid = strcmp(type, 'OF');%.*ya.*yh
        xv = ne(valid);
        [ya, xa] = gb_cell2vec(tempA.pfield_info(:,animalLoop), valid, xv);
        [yh, xh] = gb_cell2vec(tempH.pfield_info(:,animalLoop), valid, xv);
        all_a_entr = cat(1, all_a_entr, xa);
        all_h_entr = cat(1, all_h_entr, xh);
        %
        %     valid = (~isnan(ne));%.*ya.*yh
        subplot(2,4,1); hold on
        axis([minx maxx -.1 1.6])
        all_a_bits  = cat(1, all_a_bits, ya);
        all_h_bits  = cat(1, all_h_bits, yh);
        sub_err_bars(xa, ya, [ ], .05, 'quantile', acc_plot_style, acc_scatter_style)
        %         plot(xv, ya_av, acc_plot_style{:})
        set(gca, 'YTick', 0:.4:1.6, 'XTick', xticks)
        xlabel('Entrances/Minute')
        ylabel('Bits/spike')
        
        %         subplot(2,4,2); hold on
        axis([minx maxx -.1 1.6])
        sub_err_bars(xh, yh, [ ], .05, 'quantile', hpc_plot_style, hpc_scatter_style)
        set(gca, 'YTick', 0:.4:1.6, 'XTick', xticks)
        xlabel('Entrances/Minute')
        ylabel('Bits/spike')
        
        ya = tempA.percent_place(valid,animalLoop); xa = xv;
        yh = tempH.percent_place(valid,animalLoop); xh = xv;
        subplot(2,4,3); hold on
        axis([minx maxx -.1 1.1])
        all_a_place  = cat(1, all_a_place, ya);
        all_h_place  = cat(1, all_h_place, yh);
        scatter(xa, ya, acc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', xticks)
        xlabel('Entrances/Minute')
        ylabel('% Place cells')
        
        %         subplot(2,4,4); hold on
        axis([minx maxx -.1 1.1])
        scatter(xh, yh, hpc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', xticks)
        xlabel('Entrances/Minute')
        ylabel('% Place cells')
        
        [ya, xa] = gb_cell2vec(tempA.pfield_stability(:,animalLoop), valid, xv);
        [yh, xh] = gb_cell2vec(tempH.pfield_stability(:,animalLoop), valid, xv);
        subplot(2,4,5); hold on
        axis([minx maxx -.1 1.1])
        all_a_corr  = cat(1, all_a_corr, ya);
        all_h_corr  = cat(1, all_h_corr, yh);
        sub_err_bars(xa, ya, [ ], .05, 'quantile', acc_plot_style, acc_scatter_style)
        set(gca, 'YTick', 0:.2:1, 'XTick', xticks)
        xlabel('Entrances/Minute')
        ylabel('Within Sess Corr')
        
        %         subplot(2,4,6); hold on
        axis([minx maxx -.1 1.1])
        sub_err_bars(xh, yh, [ ], .05, 'quantile', hpc_plot_style, hpc_scatter_style)
        set(gca, 'YTick', 0:.2:1, 'XTick', xticks)
        xlabel('Entrances/Minute')
        ylabel('Within Sess Corr')
        
        ya = tempA.percent_stable(valid,animalLoop); xa = xv;
        yh = tempH.percent_stable(valid,animalLoop); xh = xv;
        subplot(2,4,7); hold on
        axis([minx maxx -.1 1.1])
        all_a_place  = cat(1, all_a_place, ya);
        all_h_place  = cat(1, all_h_place, yh);
        scatter(xa, ya, acc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', xticks)
        xlabel('Entrances/Minute')
        ylabel('% Stable cells')
        
        %         subplot(2,4,8); hold on
        axis([minx maxx -.1 1.1])
        scatter(xh, yh, hpc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', xticks)
        xlabel('Entrances/Minute')
        ylabel('% Stable cells')
        
        [ya, ~] = gb_cell2vec(tempA.pfield_coherence(:,animalLoop), valid, xv);
        [yh, ~] = gb_cell2vec(tempH.pfield_coherence(:,animalLoop), valid, xv);
%         ya = tempA.pfield_coherence(valid,animalLoop);
%         yh = tempH.pfield_coherence(valid,animalLoop);
        all_a_coh  = cat(1, all_a_coh, ya);
        all_h_coh  = cat(1, all_h_coh, yh);

    end
    if strloop ==1
        ACC.room.allentr = all_a_entr;
        HPC.room.allentr = all_h_entr;
        ACC.room.allcorr = all_a_corr;
        HPC.room.allcorr = all_h_corr;
        ACC.room.allcoh  = all_a_coh;
        HPC.room.allcoh  = all_h_coh;
        ACC.room.allbits = all_a_bits;
        HPC.room.allbits = all_h_bits;
        ACC.room.allplace = all_a_place;
        HPC.room.allplace = all_h_place;
    elseif strloop == 2
        ACC.arena.allentr = all_a_entr;
        HPC.arena.allentr = all_h_entr;
        ACC.arena.allcorr = all_a_corr;
        HPC.arena.allcorr = all_h_corr;
        ACC.arena.allcoh = all_a_coh;
        HPC.arena.allcoh = all_h_coh;
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
    all_a_coh  = [];
    all_h_coh  = [];
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
imagesc(h); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', cyan)
% p.EdgeColor = 'none';
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
imagesc(a); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', magenta)
% p.EdgeColor = 'none';
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
plot([nanmedian(HPC.room.allcorr-HPC.arena.allcorr), nanmedian(HPC.room.allcorr-HPC.arena.allcorr)], [-100 1600], '-', 'Color', cyan/1.5, 'LineWidth', 2)
plot([nanmedian(ACC.room.allcorr-ACC.arena.allcorr), nanmedian(ACC.room.allcorr-ACC.arena.allcorr)], [-100 1600], '--', 'Color', magenta/1.5, 'LineWidth', 2)
axis([-1.1 1.1 -100 1600])

% hh = histogram(HPC.room.allcorr-HPC.arena.allcorr, bins, 'Normalization', 'count', 'FaceColor', cyan);
% ha = histogram(ACC.room.allcorr-ACC.arena.allcorr, bins, 'Normalization', 'count', 'FaceColor', magenta);
% plot([nanmedian(HPC.room.allcorr-HPC.arena.allcorr), nanmedian(HPC.room.allcorr-HPC.arena.allcorr)], [0 700], 'b', 'LineWidth', 2)
% plot([nanmedian(ACC.room.allcorr-ACC.arena.allcorr), nanmedian(ACC.room.allcorr-ACC.arena.allcorr)], [0 700], 'm', 'LineWidth', 2)
legend([hh, ha], {'CA1' 'ACC'})
ylabel('Number of cells')
xlabel('Frame preference, corr (Room-Arena)')
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
imagesc(h); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', cyan)
% p.EdgeColor = 'none';
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
imagesc(a); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', magenta)
% p.EdgeColor = 'none';
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
xlabel('Frame preference, bits (Room-Arena)')
%%

figure(1261); clf;
bins = [0:.025:1.0];
xs = [1:20:length(bins)];
h = histcounts2(HPC.arena.allcoh, HPC.room.allcoh, bins, bins, 'Normalization', 'count');
a = histcounts2(ACC.arena.allcoh, ACC.room.allcoh, bins, bins, 'Normalization', 'count');
subplot(2,2,1); hold on
imagesc(h); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', cyan)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis([-5 length(bins)+5 -5 length(bins)+5])
axis square
colorbar
title('CA1')
xlabel('Pfield coherence (Room)')
ylabel('Pfield coherence (Arena)')
colormap magma

subplot(2,2,2); hold on
imagesc(a); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', magenta)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis square
colorbar
title('ACC')
xlabel('Pfield coherence (Room)')
ylabel('Pfield coherence (Arena)')
axis([-5 length(bins)+5 -5 length(bins)+5])

subplot(2,1,2); hold on
bins = [-1.0:.025:1.0];
hh = histogram(HPC.room.allcoh-HPC.arena.allcoh, bins, 'Normalization', 'count', 'FaceColor', cyan);
ha = histogram(ACC.room.allcoh-ACC.arena.allcoh, bins, 'Normalization', 'count', 'FaceColor', magenta);
plot([nanmedian(HPC.room.allcoh-HPC.arena.allcoh), nanmedian(HPC.room.allcoh-HPC.arena.allcoh)], [-100 1600], '-', 'Color', cyan/1.5, 'LineWidth', 2)
plot([nanmedian(ACC.room.allcoh-ACC.arena.allcoh), nanmedian(ACC.room.allcoh-ACC.arena.allcoh)], [-100 1600], '--', 'Color', magenta/1.5, 'LineWidth', 2)
axis([-1.1 1.1 -100 1500])
legend([hh, ha], {'CA1' 'ACC'})
ylabel('Number of cells')
xlabel('Frame preference, coherence (Room-Arena)')
[p, h] = signtest(HPC.room.allcoh - HPC.arena.allcoh)
[p, h] = signtest(ACC.room.allcoh - ACC.arena.allcoh)
% [h,p] = kstest(ACC.room.allcoh - ACC.arena.allcoh)
[p] = ranksum(ACC.room.allcoh - ACC.arena.allcoh, HPC.room.allcoh - HPC.arena.allcoh)
%% CORR AND SPATIAL
figure(1251); clf;
bins = [0:.025:1.0];
xs = [1:20:length(bins)];
h = histcounts2(HPC.room.allcorr, HPC.room.allbits, bins, bins, 'Normalization', 'count');
a = histcounts2(ACC.room.allcorr, ACC.room.allbits, bins, bins, 'Normalization', 'count');
subplot(2,2,1); hold on
imagesc(h); %plot([1,length(bins)-1], [1,length(bins)-1], 'Color', cyan)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis([-5 length(bins)+5 -5 length(bins)+5])
axis square
colorbar
title('CA1')
ylabel('Pfield corr (Room)')
xlabel('Pfield bits (Room)')
colormap viridis

subplot(2,2,2); hold on
imagesc(a); %plot([1,length(bins)-1], [1,length(bins)-1], 'Color', magenta)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
title('ACC')
ylabel('Pfield corr (Room)')
xlabel('Pfield bits (Room)')
axis square
colorbar
axis([-5 length(bins)+5 -5 length(bins)+5])

h = histcounts2(HPC.arena.allcorr, HPC.arena.allbits, bins, bins, 'Normalization', 'count');
a = histcounts2(ACC.arena.allcorr, ACC.arena.allbits, bins, bins, 'Normalization', 'count');
subplot(2,2,3); hold on
imagesc(h); %plot([1,length(bins)-1], [1,length(bins)-1], 'Color', cyan)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis([-5 length(bins)+5 -5 length(bins)+5])
axis square
colorbar
title('CA1')
ylabel('Pfield corr (Arena)')
xlabel('Pfield bits (Arena)')
colormap viridis

subplot(2,2,4); hold on
imagesc(a); %plot([1,length(bins)-1], [1,length(bins)-1], 'Color', magenta)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
title('ACC')
ylabel('Pfield corr (Arena)')
xlabel('Pfield bits (Arena)')
axis square
colorbar
axis([-5 length(bins)+5 -5 length(bins)+5])

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



%%
% SfN2023_APA_HPCACC_rat_summaryfigs_part2;
%%
accfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files/2023_06_15_H14_25_45_TR13_@placecells_ACC_miniscope2.mat';
hpcfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files/2023_06_15_H14_25_45_TR13_@placecells_HPC_miniscope1.mat';
% accfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files/2023_07_10_H15_35_43_OF24_@placecells_ACC_miniscope2.mat';
% hpcfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files/2023_07_10_H15_35_43_OF24_@placecells_HPC_miniscope1.mat';
hpc_ms = load(hpcfile);
acc_ms = load(accfile);

r = hpc_ms.ms.room;
a = hpc_ms.ms.arena;
r.pcell_stats.coherence = pfield_coherence_calc(r.pfields, r.vmap);
a.pcell_stats.coherence = pfield_coherence_calc(a.pfields, a.vmap);
% pull out cell examples
figure(1001); clf; hold on
subplot(1,2,1); hold on
shksx = interp1(hpc_ms.ms.timestamps, r.x, r.entranceTimes, 'linear');
shksy = interp1(hpc_ms.ms.timestamps, r.y, r.entranceTimes, 'linear');
p = patch([0, -60*cos(pi/3), 60*cos(pi/3)], [0 60*sin(pi/3) 60*sin(pi/3)], 'r'); p.FaceAlpha = .3; p.EdgeColor='none';
plot(r.x,r.y, 'Color', [.4 .4 .4], 'LineWidth', .5); axis([-44 44 -44 44])
scatter(shksx,shksy, 50, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 .5 .5])
axis square
set(gca, 'Color', 'none')

subplot(1,2,2); hold on
shksx = interp1(hpc_ms.ms.timestamps, a.x, r.entranceTimes, 'linear');
shksy = interp1(hpc_ms.ms.timestamps, a.y, r.entranceTimes, 'linear');
plot(a.x,a.y, 'Color', [.4 .4 .4], 'LineWidth', .5); axis([-44 44 -44 44])
scatter(shksx,shksy, 50, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 .5 .5])
axis square
set(gca, 'Color', 'none')

ncells = 4;

goodc = (r.pcell_stats.coherence >= quantile(r.pcell_stats.coherence, .75));
goodi = (r.pcell_stats.infoPerSpike >= quantile(r.pcell_stats.infoPerSpike, .75));
goodcr = (r.split_corr >= quantile(r.split_corr, .75));
good_r = find(goodc &  goodi & goodcr);

goodc = (a.pcell_stats.coherence >= quantile(a.pcell_stats.coherence, .75));
goodi = (a.pcell_stats.infoPerSpike >= quantile(a.pcell_stats.infoPerSpike, .75));
goodcr = (a.split_corr >= quantile(a.split_corr, .75));
good_a = find(goodc &  goodi & goodcr);
%
r_pf = [];
a_pf = [];
% rord = randperm(length(good_r));
rord = [11    12     6     8     9    10     5     3]; % for repro
% aord = randperm(length(good_a));
aord = [4     8     2     7     6     5     9    10]; % for repro
for i = 1:ncells
% r_pf = [squeeze(r.pfields_smooth_split1(good_r(i), :,:)); squeeze(r.pfields_smooth_split2(good_r(i), :,:))]; % squeeze(r.pfields_smooth(good, :,:));
% a_pf = [squeeze(a.pfields_smooth_split1(good_a(i), :,:)); squeeze(a.pfields_smooth_split2(good_a(i), :,:))]; % squeeze(a.pfields_smooth(good, :,:));
temp1 = squeeze(r.pfields_smooth(good_r(rord(i)), :,:));
temp2 = squeeze(a.pfields_smooth(good_r(rord(i)), :,:));
temp = cat(2, temp1, NaN(size(temp1,1),2), temp2);
temp = normalize_matrix(temp);
r_pf = [r_pf; NaN(4, size(temp,2)); temp];

temp1 = squeeze(r.pfields_smooth(good_a(aord(i)), :,:));
temp2 = squeeze(a.pfields_smooth(good_a(aord(i)), :,:));
temp = cat(2, temp1, NaN(size(temp1,1),2), temp2);
temp = normalize_matrix(temp)*1;
a_pf = [a_pf; NaN(4, size(temp,2)); temp];

end
a_pf = cat(2, a_pf, NaN(size(a_pf,1),1));
r_pf = cat(2, r_pf, NaN(size(r_pf,1),1));

% subplot(9,2,[3:3+(ncells-1)*2]);
figure(1002); clf
p = pcolor([r_pf a_pf]);
axis image off
p.EdgeColor = 'none';
set(gca, 'Color', 'none')
colormap viridis

hpc_ms.params.skip_ensemble = true;
% [h_ipos] = Fenton_ipos(hpc_ms.ms, .25, 'arena', hpc_ms.params);
h_ipos = abs(hpc_ms.ms.room.momentary_pos_info) - abs(hpc_ms.ms.arena.momentary_pos_info);
a_ipos = abs(acc_ms.ms.room.momentary_pos_info) - abs(acc_ms.ms.arena.momentary_pos_info);
h_coni = hpc_ms.ms.room.conjoint_ipos_av - hpc_ms.ms.arena.conjoint_ipos_av;
a_coni = acc_ms.ms.room.conjoint_ipos_av - acc_ms.ms.arena.conjoint_ipos_av;

goodinds = ~any(isnan([h_ipos; a_ipos]),1);
h_ipos = zscore( h_ipos(:, goodinds) );
a_ipos = zscore( a_ipos(:, goodinds) );
% [c, p] = xcorr(mean(h_ipos,1)', mean(a_ipos,1)');
% [c, p] = corr(mean(h_ipos,1)', mean(a_ipos,1)', 'type', 'Kendall');
% [c, p] = corr([h_ipos; a_ipos(1:10,:)*0; a_ipos]');
allipos = [h_ipos; a_ipos];
% allipos(abs(allipos)<1) = 0;
% c = NaN(size(allipos,1)); p = c;
% for i = 1:size(allipos,1)
%     for j = 1:size(allipos,1)
%     a = allipos(i,:);
%     b = allipos(j,:);
%     if sum(~isnan(a.*b), 2) >= 50;
%     [c(i,j), p(i,j)] = nancorr(a,b);
%     c(j,i) = c(i,j);
%     p(j,i) = p(i,j);
%     end
%     end
% end
[c, p] = corr(allipos');
allipos = [h_ipos; a_ipos];
% [c, p] = corr([hpc_ms.ms.room.svm_decoding.spks_bin; acc_ms.ms.room.svm_decoding.spks_bin]');
c((eye(size(c,1)))==1)= NaN;
area_label = [zeros(size(h_ipos,1),1); ones(size(a_ipos,1),1)];
cm = nansum(abs(c),1);
% cm = nansum(p<.0001,1);
[~, cmord] = sort(c(cm==max(cm),:));
% c(end:-1:size(c,1)-size(a_ipos,1), :) = 2;
% c(:, end:-1:size(c,1)-size(a_ipos,1)) = 2;
c = c(cmord,:);
c = c(:,cmord);
allipos = allipos(cmord,:);
area_label = area_label(cmord);
% c = cat(2, cat(1, c, area_label'), [area_label;0]);
figure; imagesc(c, [-.4 .4])
set(gca, 'XTick', find(area_label==0), 'YTick', find(area_label==1), 'TickDir', 'out', 'XTickLabel', 'h', 'YTickLabel', 'a');
figure; imagesc(allipos, [-2 2]);
set(gca, 'YTick', find(area_label==1), 'TickDir', 'out', 'YTickLabel', 'a');
% SHOW WHICH REGIONS ARE INTERSPERSED
r = acc_ms.ms.room;
a = acc_ms.ms.arena;
r.pcell_stats.coherence = pfield_coherence_calc(r.pfields, r.vmap);
a.pcell_stats.coherence = pfield_coherence_calc(a.pfields, a.vmap);
% pull out cell examples


goodc = (r.pcell_stats.coherence >= quantile(r.pcell_stats.coherence, .75));
goodi = (r.pcell_stats.infoPerSpike >= quantile(r.pcell_stats.infoPerSpike, .75));
goodcr = (r.split_corr >= quantile(r.split_corr, .75));
good_r = find(goodc &  goodi & goodcr);

goodc = (a.pcell_stats.coherence >= quantile(a.pcell_stats.coherence, .75));
goodi = (a.pcell_stats.infoPerSpike >= quantile(a.pcell_stats.infoPerSpike, .75));
goodcr = (a.split_corr >= quantile(a.split_corr, .75));
good_a = find(goodc &  goodi & goodcr);
%
r_pf = [];
a_pf = [];
% rord = randperm(length(good_r))
rord = [3     4     1     2     5]; % for repro
% aord = randperm(length(good_a))
aord = [8     6     2     4     5     3     7     1]; % for repro
for i = 1:ncells
% r_pf = [squeeze(r.pfields_smooth_split1(good_r(i), :,:)); squeeze(r.pfields_smooth_split2(good_r(i), :,:))]; % squeeze(r.pfields_smooth(good, :,:));
% a_pf = [squeeze(a.pfields_smooth_split1(good_a(i), :,:)); squeeze(a.pfields_smooth_split2(good_a(i), :,:))]; % squeeze(a.pfields_smooth(good, :,:));
temp1 = squeeze(r.pfields_smooth(good_r(rord(i)), :,:));
temp2 = squeeze(a.pfields_smooth(good_r(rord(i)), :,:));
temp = cat(2, temp1, NaN(size(temp1,1),2), temp2);
temp = normalize_matrix(temp);
r_pf = [r_pf; NaN(4, size(temp,2)); temp];

temp1 = squeeze(r.pfields_smooth(good_a(aord(i)), :,:));
temp2 = squeeze(a.pfields_smooth(good_a(aord(i)), :,:));
temp = cat(2, temp1, NaN(size(temp1,1),2), temp2);
temp = normalize_matrix(temp)*1;
a_pf = [a_pf; NaN(4, size(temp,2)); temp];

end
a_pf = cat(2, a_pf, NaN(size(a_pf,1),1));
r_pf = cat(2, r_pf, NaN(size(r_pf,1),1));
figure(1003); clf
% subplot(9,2,[12:12+(ncells-1)*2]);
p = pcolor([r_pf a_pf]);
axis image off
p.EdgeColor = 'none';
colormap plasma
set(gca, 'Color', 'none')
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


