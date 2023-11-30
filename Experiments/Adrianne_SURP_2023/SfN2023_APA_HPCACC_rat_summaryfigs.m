% analysis ideas pre-sfn:
% - look at activiy peri-avoidance and shocks
% - extract whole frame fluorescence to get the big flashing activiy?
% - make sure that interpolated frames from bad frames are excluded from
% analysis
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
analysis_version = 'v2.01';
this_ver = str2double(analysis_version(2:end));

time_origin = datetime(2020, 1, 1, 0, 0, 0); %

numAnimals = length(animals);
sfep_occthresh = 3; 
sfep_smoothingsize = 3; % gaussian smoothing radius for fields
plotSFEP = false;
reprocessonly = true;
%
% file: Hipp18240 - 2022_09_30_H18_50_38_TR27_@placecells.mat
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
    HPC.room.sfep.nullabsdiff = temp; % summed difference in the SFEP map compared to a spatial only firing model
    HPC.room.sfep.nulldiff = temp; % summed abs difference in the SFEP map compared to a spatial only firing model
    HPC.room.sfep.mapcorr = temp; % corr of the SFEP map compared to a spatial only input
    HPC.room.sfep.coherence = temp;
    HPC.room.sfep.bits = temp;
    HPC.room.sfep.coherence_infer = temp;
    HPC.room.sfep.bits_infer = temp;
    HPC.ncells = temp;
    HPC.ipos_runsz = temp;
    HPC.ipos_av = temp;
    HPC.avoidpeth = [NaN(121,1)];
    HPC.shkpeth = [NaN(121,1)];
    
    temp = cell(40,numAnimals);
    
    HPC.fr = temp;
    HPC.room.sfep.map = temp;
    HPC.room.sfep.infer_map = temp;
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
    %%
    for animalLoop = 1:numAnimals
%         if any(animalLoop==[1 6])
%             plotSFEP = true;
%         else
%             plotSFEP = false;
%         end
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
                file_ver = str2double(prev_version.analysis_version(2:end));
                
                if file_ver < this_ver
                    rerun_ipos = true;
                else
                    rerun_ipos = false;
                end
                    fprintf('~~~file: %s...', fname)
                    fprintf(' analysis_version : %1.2f\n', file_ver)
                    clearvars ms
                    load(processedFile, 'ms');
                    prevparams = load(processedFile, 'params', 'analysis_version');
                    [sessDate, trialname, trial_type, trial_num] = recdata_from_parentdir(ms.parentDir);
                    if ~isempty(strfind(ms.parentDir, 'Hipp18240/2022_09_08/16_09_56_HAB/'))
                        trial_num = 0;
                    end
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
                    if rerun_ipos == true
                        disp('old ver, rerun ipos')
                        disp(prev_version)
                        params = ms.params;
                        params.skip_ensemble = false;
                        
                        
                        params.yaw_bins = -pi:pi/8:pi;
                        params.rho_bins = [0 10 20 30 45];
                        params.num_random_shuffle_pcell  = 100;
                        params.num_random_shuffle_decode = 100;
                        params.min_spd_thresh   = -1;
                        params.min_samples      = 5; % samples
                        params.occupancy_thresh = 2; % seconds
                        params.parfor_progbar = false;
                        params.ipos_int_time = .25;
                        ms.room.svm_decoding = []; ms.arena.svm_decoding = [];
                        fprintf('\t\tstart time: \t%s  \n', datetime())
                        [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
                        fprintf('\t\tdecode end time: \t\t%s  \n', datetime())
                        [ms.room.momentary_pos_info, ~, ms.room.conjoint_ipos_min, ms.room.conjoint_ipos_av] = ...
                            Fenton_ipos_testing(ms, params.ipos_int_time, 'room_polar',  params);
                        [ms.arena.momentary_pos_info, ~, ms.arena.conjoint_ipos_min, ms.arena.conjoint_ipos_av] = ...
                            Fenton_ipos_testing(ms, params.ipos_int_time, 'room_polar', params);
%                         if isfield(ms.room, 'momentary_pos_info_px')
%                             disp('rm momentary_pos_info_px')
%                             ms.room  = rmfield(ms.room, 'momentary_pos_info_px');
%                             ms.arena = rmfield(ms.arena, 'momentary_pos_info_px');
%                         end
                        disp(analysis_version)
                        fprintf('\t\tipos end time: \t\t%s  \n', datetime())
                        fprintf('SAVING - %s\n\n', processedFile)
                        save(processedFile, 'ms', 'analysis_version', '-append');
%                         save(processedFile, 'ms', 'analysis_version', '-append', '-v7.3');
                    end
                    if reprocessonly==false
                    r = ms.room;
                    a = ms.arena;
                    if strcmp(trial_type, 'TR')
                        shock_zone_center = pi/2; % typical room shock configuration
                        shock_zone_size = pi/6; % size in rad from center to edge
                        distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure
                    elseif strcmp(trial_type, 'CON')
                        shock_zone_center = 3*pi/2; % typical room shock configuration
                        shock_zone_size = pi/6; % size in rad from center to edge
                        distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure
                    end
                    ipos = abs(r.momentary_pos_info) - abs(a.momentary_pos_info);
                    
                    [avoid_ind, shks_ind] = shock_zone_proximity(ms,...
                        shock_zone_center, shock_zone_size, distance_entrance_size, false);
                    isavoid = false(size(r.svm_decoding.x)); isavoid(r.svm_decoding.spks_bin_group(avoid_ind)) = true;
                    isshk = false(size(r.svm_decoding.x)); isshk(r.svm_decoding.spks_bin_group(shks_ind)) = true;
                    ipos_mean= nanmean(ipos,1);
                    [avoid_peth] = gb_PETH(ipos_mean, isavoid, 90, 30);
                    [shk_peth] = gb_PETH(ipos_mean, isshk, 90, 30);
                    
                    
                    [theta, rho] = cart2pol(r.x, r.y);
%                     prop_inside = sum(theta>pi/3 & theta<2*pi/3)/length(theta);
                    prop_inside(thisind, animalLoop) = sum(theta>=pi/3 & theta<=2*pi/3)/length(theta);
                    
                    [~, r.pcell_stats.coherence] = pfield_coherence_calc(r.pfields, r.vmap); % [coh, coh_z]
                    [~, a.pcell_stats.coherence] = pfield_coherence_calc(a.pfields, a.vmap); % [coh, coh_z]
%                     concordant_pfield_change(ms, prevparams.params);
                    goodinds = any(~isnan(ipos),1);
                    iposz = zscore( ipos(:, goodinds) );
                    ipos_mean= nanmean(ipos(:, goodinds),1);
                    [~,~,stats] = runstest(ipos_mean);
                    ipos_runsz = stats.z;
                    ipos_av = nanmean(abs(ipos_mean));
                    

                    
                    [sfep, room_inferred_ipos, arena_inferred_ipos] = SFEP_map_quantification(ms, sfep_smoothingsize, true, true);
%                     [sfep, room_inferred_ipos, arena_inferred_ipos] = CONIPOS_map_quantification(ms, smoothing_size, false, true);
                    if ~isfield(ms.room, 'inferred_ipos') && ~isfield(ms.arena, 'inferred_ipos')
                        ms.room.inferred_ipos  = room_inferred_ipos;
                        ms.arena.inferred_ipos = arena_inferred_ipos;
                        ms.room.pcell_stats.coherence = r.pcell_stats.coherence;
                        ms.arena.pcell_stats.coherence = a.pcell_stats.coherence;
                        save(processedFile, 'ms', '-append')
                    end
                    if strcmp(trial_type, 'CON') % rotate the map so the shock is on top as in TR
                        sfep.room_map           = imrotate(sfep.room_map, 180);
                        sfep.room_map_smooth    = imrotate(sfep.room_map_smooth, 180);
                        sfep.arena_map          = imrotate(sfep.arena_map, 180);
                        sfep.arena_map_smooth   = imrotate(sfep.arena_map_smooth, 180);
                    end
                    room_sfep_map           = sfep.room_map;
                    room_sfep_map_smooth    = sfep.room_map_smooth;
                    room_sfep_coherence     = sfep.room_coh;
                    room_sfep_bits          = sfep.room_bits;
                    room_sfep_coherence_infer= sfep.infer.room_coh;
                    room_sfep_bits_infer    = sfep.infer.room_bits;
                    arena_sfep_map          = sfep.arena_map;
                    arena_sfep_map_smooth   = sfep.arena_map_smooth;
                    arena_sfep_coherence    = sfep.arena_coh;
                    arena_sfep_bits         = sfep.arena_bits;
                    arena_sfep_coherence_infer= sfep.infer.arena_coh;
                    arena_sfep_bits_infer   = sfep.infer.arena_bits;
                    
                    
                    r.ipos_roompref = nansum(ipos_mean>0)./sum((ipos_mean<0) | (ipos_mean>0));
                    a.ipos_roompref = nansum(ipos_mean<0)./sum((ipos_mean<0) | (ipos_mean>0));
                    %                     c = NaN(size(ipos,1), 1); p=c;
                    %                     for i =1:size(ipos,1); [c(i), p(i)] = circ_corrcl(theta(goodinds), ipos(i,:)); end
                    if plotSFEP==true
                        figure(animalLoop*1000+trial_num);
                        subplot(4,1,1)
                        plot(theta(goodinds)); axis tight
                        subplot(4,1,2:4);
                        imagesc(iposz, [-8 8]);
                        colormap redblue
                        figure(animalLoop*100+trial_num);
                        colormap magma
                        subplot(3,4,6); plot(r.x,r.y,'k'); axis([-44 44 -44 44]);  axis off
                        subplot(3,4,8); plot(a.x,a.y,'k'); axis([-44 44 -44 44]);  axis off
                        subplot(3,4,10); imagesc(r.vmap); set(gca,'YDir', 'normal'); axis image off;
                        subplot(3,4,12); imagesc(a.vmap); set(gca,'YDir', 'normal'); axis image off;
                    end
                    if strcmp(thisregion, 'HPC')
                        if plotSFEP==true
                            subplot(3,1,1); hold on; plot(ms.room.svm_decoding.t(goodinds), ipos_mean, 'b');
                            subplot(3,4,5); imagesc(room_sfep_map_smooth); set(gca,'YDir', 'normal'); axis image off;
                            ylabel('HPC')
                            title(sprintf('bits=%2.2f  coh=%2.2f', room_sfep_bits, room_sfep_coherence))
                            subplot(3,4,7); imagesc(arena_sfep_map_smooth); set(gca,'YDir', 'normal'); axis image off;
                            title(sprintf('bits=%2.2f  coh=%2.2f', arena_sfep_bits, arena_sfep_coherence))
                        end
                        HPC.avoidpeth = [HPC.avoidpeth avoid_peth];
                        HPC.shkpeth = [HPC.shkpeth shk_peth];
                        HPC.fr{thisind, animalLoop} = nansum(ms.spks>0,2)./sessTime(thisind, animalLoop);
                        HPC.ncells(thisind, animalLoop) = size(ms.spks,1);
                        HPC.ipos_runsz(thisind, animalLoop) = ipos_runsz;
                        HPC.ipos_av(thisind, animalLoop) = ipos_av;
                        HPC.room.sfep.map{thisind, animalLoop}  = room_sfep_map_smooth;
                        HPC.arena.sfep.map{thisind, animalLoop} = arena_sfep_map_smooth;
                        HPC.room.sfep.infer_map{thisind, animalLoop}  = sfep.infer.room_map_smooth;
                        HPC.arena.sfep.infer_map{thisind, animalLoop} = sfep.infer.arena_map_smooth;
                        HPC.room.sfep.mapcorr(thisind, animalLoop)  = sfep.room_corr;
                        HPC.arena.sfep.mapcorr(thisind, animalLoop) = sfep.arena_corr;
                        HPC.room.sfep.nulldiff(thisind, animalLoop)  = sfep.room_diff;
                        HPC.arena.sfep.nulldiff(thisind, animalLoop) = sfep.arena_diff;
                        HPC.room.sfep.nullabsdiff(thisind, animalLoop)  = sfep.room_absdiff;
                        HPC.arena.sfep.nullabsdiff(thisind, animalLoop) = sfep.arena_absdiff;
                        HPC.room.sfep.coherence(thisind, animalLoop)  = room_sfep_coherence;
                        HPC.arena.sfep.coherence(thisind, animalLoop) = arena_sfep_coherence;
                        HPC.room.sfep.bits(thisind, animalLoop)  = room_sfep_bits;
                        HPC.arena.sfep.bits(thisind, animalLoop) = arena_sfep_bits;
                        HPC.room.sfep.coherence_infer(thisind, animalLoop)  = room_sfep_coherence_infer;
                        HPC.arena.sfep.coherence_infer(thisind, animalLoop) = arena_sfep_coherence_infer;
                        HPC.room.sfep.bits_infer(thisind, animalLoop)  = room_sfep_bits_infer;
                        HPC.arena.sfep.bits_infer(thisind, animalLoop) = arena_sfep_bits_infer;
                    elseif strcmp(thisregion, 'ACC')
                        if plotSFEP==true
                            subplot(3,1,1); hold on; plot(ms.room.svm_decoding.t(goodinds), ipos_mean, 'm');
                            subplot(3,4,9); imagesc(room_sfep_map_smooth); set(gca,'YDir', 'normal'); axis image off;
                            ylabel('ACC')
                            title(sprintf('bits=%2.2f  coh=%2.2f', room_sfep_bits, room_sfep_coherence))
                            subplot(3,4,11); imagesc(arena_sfep_map_smooth); set(gca,'YDir', 'normal'); axis image off;
                            title(sprintf('bits=%2.2f  coh=%2.2f', arena_sfep_bits, arena_sfep_coherence))
                        end
                        ACC.avoidpeth = [ACC.avoidpeth avoid_peth];
                        ACC.shkpeth = [ACC.shkpeth shk_peth];
                        ACC.fr{thisind, animalLoop} = nansum(ms.spks>0,2)./sessTime(thisind, animalLoop);
                        ACC.ncells(thisind, animalLoop) = size(ms.spks,1);
                        ACC.ipos_runsz(thisind, animalLoop) = ipos_runsz;
                        ACC.ipos_av(thisind, animalLoop) = ipos_av;
                        ACC.room.sfep.map{thisind, animalLoop}  = room_sfep_map_smooth;
                        ACC.arena.sfep.map{thisind, animalLoop} = arena_sfep_map_smooth;
                        ACC.room.sfep.infer_map{thisind, animalLoop}  = sfep.infer.room_map_smooth;
                        ACC.arena.sfep.infer_map{thisind, animalLoop} = sfep.infer.arena_map_smooth;
                        ACC.room.sfep.mapcorr(thisind, animalLoop)  = sfep.room_corr;
                        ACC.arena.sfep.mapcorr(thisind, animalLoop) = sfep.arena_corr;
                        ACC.room.sfep.nulldiff(thisind, animalLoop)  = sfep.room_diff;
                        ACC.arena.sfep.nulldiff(thisind, animalLoop) = sfep.arena_diff;
                        ACC.room.sfep.nullabsdiff(thisind, animalLoop)  = sfep.room_absdiff;
                        ACC.arena.sfep.nullabsdiff(thisind, animalLoop) = sfep.arena_absdiff;
                        ACC.room.sfep.coherence(thisind, animalLoop)  = room_sfep_coherence;
                        ACC.arena.sfep.coherence(thisind, animalLoop) = arena_sfep_coherence;
                        ACC.room.sfep.bits(thisind, animalLoop)  = room_sfep_bits;
                        ACC.arena.sfep.bits(thisind, animalLoop) = arena_sfep_bits;
                        ACC.room.sfep.coherence_infer(thisind, animalLoop)  = room_sfep_coherence_infer;
                        ACC.arena.sfep.coherence_infer(thisind, animalLoop) = arena_sfep_coherence_infer;
                        ACC.room.sfep.bits_infer(thisind, animalLoop)  = room_sfep_bits_infer;
                        ACC.arena.sfep.bits_infer(thisind, animalLoop) = arena_sfep_bits_infer;
                    end
                        drawnow
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
                    end
%                     for vloop = 1:length(sfep_vars)
%                         if iscell(eval(sprintf('%s.room.%s', thisregion, sfep_var_outname{vloop})))
%                             str1 = sprintf('%s.room.%s{thisind, animalLoop}  = r.%s;', thisregion, sfep_var_outname{vloop}, sfep_vars{vloop});
%                             str2 = sprintf('%s.arena.%s{thisind, animalLoop} = a.%s;', thisregion, sfep_var_outname{vloop}, sfep_vars{vloop});
%                         else
%                             str1 = sprintf('%s.room.%s(thisind, animalLoop)  = nanmean(r.%s);', thisregion, sfep_var_outname{vloop}, sfep_vars{vloop});
%                             str2 = sprintf('%s.arena.%s(thisind, animalLoop) = nanmean(a.%s);', thisregion, sfep_var_outname{vloop}, sfep_vars{vloop});
%                         end
%                         eval(str1)
%                         eval(str2)
%                     end
%                 else
%                     fprintf('~~~PROCESSED file RERUN skipped: %s\n', fname)
%                     fprintf('\tsame analysis_version : %s\n', prev_version.analysis_version)
%                 end
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
hpc_plot_style    = {'-', 'Color', cyan, 'LineWidth', 1}; %#ok<*NASGU>
acc_plot_style    = {'-', 'Color', magenta, 'LineWidth', 1};

entrperMin = numEntr./(sessTime./60000);
isTR = strcmp(exp_type, 'TR');
isCON = strcmp(exp_type, 'CON');
all_valid = ~isnan(entrperMin) & (isTR|isCON);

HPC.room.BaysPerf   = 1 - HPC.room.BayDec_mean./HPC.room.BayDec_mean_rand;
HPC.arena.BaysPerf  = 1 - HPC.arena.BayDec_mean./HPC.arena.BayDec_mean_rand;
HPC.room.SVMPerf    = 1 - HPC.room.SVMDec_mean./HPC.room.SVMDec_mean_rand;
HPC.arena.SVMPerf   = 1 - HPC.arena.SVMDec_mean./HPC.arena.SVMDec_mean_rand;
ACC.room.BaysPerf   = 1 - ACC.room.BayDec_mean./ACC.room.BayDec_mean_rand;
ACC.arena.BaysPerf  = 1 - ACC.arena.BayDec_mean./ACC.arena.BayDec_mean_rand;
ACC.room.SVMPerf    = 1 - ACC.room.SVMDec_mean./ACC.room.SVMDec_mean_rand;
ACC.arena.SVMPerf   = 1 - ACC.arena.SVMDec_mean./ACC.arena.SVMDec_mean_rand;

var_outname = {'BaysPerf', 'SVMPerf', 'pfield_stability_av', 'pfield_info_av',...
    'percent_place', 'percent_stable', 'pfield_coh_av'}; %, 'ipos_roompref', 'sfep.coherence', 'sfep.bits'};
nvars = length(var_outname);
for rrr = 1:nvars
    eval(sprintf('anova_%s = [];', var_outname{rrr}))
%     eval(sprintf('anova_%s_label = [];', var_outname{rrr}))
    eval(sprintf('subj_anova_%s = [];', var_outname{rrr}))
%     eval(sprintf('subj_anova_%s_label = [];', var_outname{rrr}))
end
reg = {'HPC', 'ACC'};
framecolor = {'r', 'b'};
maxys = ones(nvars,1) + .1;
maxys(1:2) = [.8 .8];
maxys(4) = [1.3]; %#ok<*NBRAK>
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
            %%
                eval(sprintf('temp = %s.%s.%s;', reg{r}, top_struct{rr}, var_outname{rrr}))
                eval(sprintf('temp_anova = anova_%s;', var_outname{rrr}));
%                 eval(sprintf('temp_anova_label = anova_%s_label;', var_outname{rrr}));
                eval(sprintf('subj_temp_anova = subj_anova_%s;', var_outname{rrr}));
%                 eval(sprintf('subj_temp_anova_label = subj_anova_%s_label;', var_outname{rrr}));
                
                tempv = temp;
                tempv(~all_valid) = NaN;
                tempv(~isfinite(tempv)) = NaN;
                allsess_vals = tempv(:);
                temp_anova = cat(1, temp_anova, allsess_vals);
%                 temp_anova_label = cat(1, temp_anova_label, allsess_vals*0 + anova_label_num);
                eval(sprintf('anova_%s = temp_anova;', var_outname{rrr}));
%                 eval(sprintf('anova_%s_label = temp_anova_label;', var_outname{rrr}));

                all_av = nanmean(tempv(:));
                all_var = nanstd(tempv(:));
                ncomps = sum(isfinite(tempv(:)));
                if rrr==1
                    eval(sprintf('%s.%s.ncomps = ncomps;', reg{r}, top_struct{rr}));
                end
                
                av = nanmean(tempv,1);
                vars = nanstd(tempv);
                xs = find(isfinite(av)) + frame_offset;
                av = av(isfinite(av));
                vars = vars(isfinite(vars));
                subj_temp_anova = cat(1, subj_temp_anova, av');
%                 subj_temp_anova_label = cat(1, subj_temp_anova_label, av'*0 + anova_label_num);
                eval(sprintf('subj_anova_%s = subj_temp_anova;', var_outname{rrr}));
%                 eval(sprintf('subj_anova_%s_label = subj_temp_anova_label;', var_outname{rrr}));
                
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
                scatterxs = gb_rand_jitter(tempv(isfinite(tempv)), 15);
                scatter(scatterxs+xs + .05, tempv(isfinite(tempv)), style{:});%, 'MarkerEdgeColor', 'k')
                errorbar(xs(1), all_av, all_var, 'k', 'LineWidth', 2)
                axis([.5 2.5 minys(rrr), maxys(rrr)])
                set(gca, 'XTick', [1 2], 'XTickLabel', {'room' 'arena'}, 'Color', 'none')
                varname = sprintf('%s', var_outname{rrr});
                title(varname, 'Interpreter', 'none', 'Color', 'k')
        end
    end
end
drawnow;
%%
HPC.room.sfep.cohdiff = HPC.room.sfep.coherence - HPC.room.sfep.coherence_infer;
HPC.arena.sfep.cohdiff = HPC.arena.sfep.coherence - HPC.arena.sfep.coherence_infer;
ACC.room.sfep.cohdiff = ACC.room.sfep.coherence - ACC.room.sfep.coherence_infer;
ACC.arena.sfep.cohdiff = ACC.arena.sfep.coherence - ACC.arena.sfep.coherence_infer;
HPC.room.sfep.bitdiff = HPC.room.sfep.bits - HPC.room.sfep.bits_infer;
HPC.arena.sfep.bitdiff = HPC.arena.sfep.bits - HPC.arena.sfep.bits_infer;
ACC.room.sfep.bitdiff = ACC.room.sfep.bits - ACC.room.sfep.bits_infer;
ACC.arena.sfep.bitdiff = ACC.arena.sfep.bits - ACC.arena.sfep.bits_infer;
% HPC.room.sfep.nullcorr = 1 - HPC.room.sfep.mapcorr;
% HPC.arena.sfep.nullcorr = 1 - HPC.arena.sfep.mapcorr;
% ACC.room.sfep.nullcorr = 1 - ACC.room.sfep.mapcorr;
% ACC.arena.sfep.nullcorr = 1 - ACC.arena.sfep.mapcorr;
% var_outname = {'sfep.coherence', 'sfep.bits', 'sfep.coherence_infer', 'sfep.bits_infer',...
%     'sfep.nulldiff', 'sfep.nullabsdiff', 'sfep.mapcorr', 'ipos_roompref'};
var_outname = {'sfep.coherence', 'sfep.bits', 'sfep.cohdiff', 'sfep.bitdiff',...
    'sfep.nulldiff', 'sfep.mapcorr', 'ipos_roompref'};
nvars = length(var_outname);
anova_label = [];
subj_anova_label = [];
for rrr = 1:nvars
    eval(sprintf('anova_%s = [];', var_outname{rrr}))
%     eval(sprintf('anova_%s_label = [];', var_outname{rrr}))
    eval(sprintf('subj_anova_%s = [];', var_outname{rrr}))
%     eval(sprintf('subj_anova_%s_label = [];', var_outname{rrr}))
end
reg = {'HPC', 'ACC'};
framecolor = {'r', 'b'};
maxys = ones(nvars,1) + .1;
maxys([2,4]) = [1.3]; % bits
maxys([5]) = [.7]; % sfep map diff
minys = zeros(nvars,1) - .1;
minys([5]) = [-.7]; % sfep map diff
minys(3) = [-.5]; % sfep map diff
minys([7]) = [-.5]; % sfep map diff
% top_struct = {'room', 'arena'};
figure(9954); clf
figure(9955); clf

anova_label_num = 0;
% 1 = HPC_room_var
% 2 = HPC_arena_var
% 3 = ACC_room_var
% 4 = ACC_arena_var
for r = 1:length(reg)
    figure(9945+r); clf
    if r==1; regcolor = cyan; style = hpc_scatter_style;  elseif r==2; regcolor = magenta;  style = acc_scatter_style; end
    reg_offset = ( r./(1+length(reg)) - .5 ) ;
    for rr = 1:length(top_struct)
        frame_offset = ( rr./(1+length(top_struct)) - .5 )/2 ;
        anova_label_num = anova_label_num+1;
        for rrr = 1:nvars
            %%
                eval(sprintf('temp = %s.%s.%s;', reg{r}, top_struct{rr}, var_outname{rrr}))
                eval(sprintf('temp_anova = anova_%s;', var_outname{rrr}));
%                 eval(sprintf('temp_anova_label = anova_%s_label;', var_outname{rrr}));
                eval(sprintf('subj_temp_anova = subj_anova_%s;', var_outname{rrr}));
%                 eval(sprintf('subj_temp_anova_label = subj_anova_%s_label;', var_outname{rrr}));
                
                tempv = temp;
                tempv(~all_valid) = NaN;
                tempv(~isfinite(tempv)) = NaN;
                allsess_vals = tempv(:);
                temp_anova = cat(1, temp_anova, allsess_vals);
                eval(sprintf('anova_%s = temp_anova;', var_outname{rrr}));
                
                if rrr==1
                eval(sprintf('temp_anova_label = anova_label;'));
                temp_anova_label = cat(1, temp_anova_label, allsess_vals*0 + anova_label_num);
                eval(sprintf('anova_label = temp_anova_label;'));
                
                eval(sprintf('subj_temp_anova_label = subj_anova_label;'));
                subj_temp_anova_label = cat(1, subj_temp_anova_label, av'*0 + anova_label_num);
                eval(sprintf('subj_anova_label = subj_temp_anova_label;'));
                end

                all_av = nanmean(tempv(:));
                all_var = nanstd(tempv(:));
                ncomps = sum(isfinite(tempv(:)));
                if rrr==1
                    eval(sprintf('%s.%s.ncomps = ncomps;', reg{r}, top_struct{rr}));
                end
                
                av = nanmean(tempv,1);
                vars = nanstd(tempv);
                xs = find(isfinite(av)) + frame_offset;
                av = av(isfinite(av));
                vars = vars(isfinite(vars));
                subj_temp_anova = cat(1, subj_temp_anova, av');
%                 subj_temp_anova_label = cat(1, subj_temp_anova_label, av'*0 + anova_label_num);
                eval(sprintf('subj_anova_%s = subj_temp_anova;', var_outname{rrr}));
%                 eval(sprintf('subj_anova_%s_label = subj_temp_anova_label;', var_outname{rrr}));
                
                ind = sub2ind([ceil(nvars/2),2],rrr);
                figure(9945+r);
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
                
                figure(9955);
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
                
                figure(9954);
                subplot(ceil(nvars/2), 2, ind); hold on
                xs = zeros(ncomps,1)+reg_offset + rr;
                bar(xs(1), all_av, 'FaceColor', regcolor/1.5, 'BarWidth', .25)
                scatterxs = gb_rand_jitter(tempv(isfinite(tempv)), 15);
                scatter(scatterxs+xs + .05, tempv(isfinite(tempv)), style{:});%, 'MarkerEdgeColor', 'k')
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
c(~all_valid) = NaN;
nc = nanmean(c,1);
vc = nanstd(c,1);
figure(1114); clf; hold on
bar([1:numAnimals]-.15, nc(ncellsubjorder), 'FaceColor', cyan, 'BarWidth', .2)
for i = 1:length(nc)
    errorbar(i-.15, nc(ncellsubjorder(i)), vc(ncellsubjorder(i)), 'Color', 'k')
end

c = ACC.ncells;
c(~all_valid) = NaN;
nc = nanmean(c,1);
vc = nanstd(c,1);
bar([1:numAnimals]+.15, nc(ncellsubjorder), 'FaceColor', magenta, 'BarWidth', .2)
for i = 1:length(nc)
    errorbar(i+.15, nc(ncellsubjorder(i)), vc(ncellsubjorder(i)), 'Color', 'k')
end
axis([.15 numAnimals+.85 0 900])
an = animals(ncellsubjorder);
set(gca, 'YTick', [0:200:800], 'XTick', 1:numAnimals, 'XTickLabel', an(:), 'XTickLabelRotation', -90);

nsessvald = sum(all_valid,1);
for i = 1:length(nsessvald)
    text(i-.25, -20, sprintf('%d sess', nsessvald(ncellsubjorder(i))), 'Rotation', -90, 'FontSize', 8)
end
ylabel('Average cells recoded per session')
%%%%%%%%%%%%%%%%%%%%%%%%%%%% IPOS mean abs av 

figure(1226); clf; 
subplot(1,2,1); hold on;
tempv = HPC.ipos_av;
tempv(~all_valid) = NaN; hv = tempv;
all_av = nanmean(tempv(:),1);
all_var = nanstd(tempv(:),1);
regcolor = cyan; style = hpc_scatter_style;  
reg_offset = ( 1./(1+length(reg)) - .5 ) ;
scatterxs = gb_rand_jitter(tempv(isfinite(tempv)), 15);
xs = zeros(length(scatterxs),1)+ reg_offset;
bar(xs(1), all_av, 'FaceColor', regcolor/1.5, 'BarWidth', .25)
scatter(scatterxs+xs + .05, tempv(isfinite(tempv)), style{:});%, 'MarkerEdgeColor', 'k')
errorbar(xs(1), all_av, all_var, 'k', 'LineWidth', 2)

tempv = ACC.ipos_av;
tempv(~all_valid) = NaN; av = tempv;
all_av = nanmean(tempv(:),1);
all_var = nanstd(tempv(:),1);
regcolor = magenta;  style = acc_scatter_style;
reg_offset = ( 2./(1+length(reg)) - .5 ) ;
scatterxs = gb_rand_jitter(tempv(isfinite(tempv)), 15);
xs = zeros(length(scatterxs),1)+ reg_offset;
bar(xs(1), all_av, 'FaceColor', regcolor/1.5, 'BarWidth', .25)
scatter(scatterxs+xs + .05, tempv(isfinite(tempv)), style{:});%, 'MarkerEdgeColor', 'k')
errorbar(xs(1), all_av, all_var, 'k', 'LineWidth', 2)
axis([-.5 .5 -.01, .05])

subplot(1,2,2); hold on;
tempv = HPC.ipos_av;
tempv(~all_valid) = NaN;
tempv = nanmean(tempv,1); hvs = tempv;
all_av = nanmean(tempv(:));
all_var = nanstd(tempv(:));
regcolor = cyan; style = hpc_scatter_style;  
reg_offset = ( 1./(1+length(reg)) - .5 ) ;
scatterxs = gb_rand_jitter(tempv(isfinite(tempv)), 10);
xs = zeros(length(scatterxs),1)+ reg_offset;
bar(xs(1), all_av, 'FaceColor', regcolor/1.5, 'BarWidth', .25)
scatter(scatterxs+xs + .05, tempv(isfinite(tempv)), style{:}, 'MarkerEdgeColor', 'k')
errorbar(xs(1), all_av, all_var, 'k', 'LineWidth', 2)

tempv = ACC.ipos_av;
tempv(~all_valid) = NaN;
tempv = nanmean(tempv,1); avs = tempv;
all_av = nanmean(tempv(:));
all_var = nanstd(tempv(:));
regcolor = magenta;  style = acc_scatter_style;
reg_offset = ( 2./(1+length(reg)) - .5 ) ;
scatterxs = gb_rand_jitter(tempv(isfinite(tempv)), 10);
xs = zeros(length(scatterxs),1)+ reg_offset;
bar(xs(1), all_av, 'FaceColor', regcolor/1.5, 'BarWidth', .25)
scatter(scatterxs+xs + .05, tempv(isfinite(tempv)), style{:}, 'MarkerEdgeColor', 'k')
errorbar(xs(1), all_av, all_var, 'k', 'LineWidth', 2)
axis([-.5 .5 -.01, .05])

ranksum(hv(~isnan(hv)),av(~isnan(av)))
ranksum(hvs(~isnan(hvs)),avs(~isnan(avs)))
% set(gca, 'XTick', [1 2], 'XTickLabel', {'room' 'arena'}, 'Color', 'none')


figure(1226); clf; hold on
bar([1:numAnimals]-.15, nc(ncellsubjorder), 'FaceColor', cyan, 'BarWidth', .2)
for i = 1:length(nc)
    errorbar(i-.15, nc(ncellsubjorder(i)), vc(ncellsubjorder(i)), 'Color', 'k')
end
c = ACC.ipos_av;
c(~all_valid) = NaN;
nc = nanmean(c,1);
vc = nanstd(c,1);
bar([1:numAnimals]+.15, nc(ncellsubjorder), 'FaceColor', magenta, 'BarWidth', .2)
for i = 1:length(nc)
    errorbar(i+.15, nc(ncellsubjorder(i)), vc(ncellsubjorder(i)), 'Color', 'k')
end
axis([.15 numAnimals+.85 0 .05])
an = animals(ncellsubjorder);
set(gca, 'YTick', [0:.01:.05], 'XTick', 1:numAnimals, 'XTickLabel', an(:), 'XTickLabelRotation', -90);

nsessvald = sum(all_valid,1);
for i = 1:length(nsessvald)
    text(i-.25, -20, sprintf('%d sess', nsessvald(ncellsubjorder(i))), 'Rotation', -90, 'FontSize', 8)
end
ylabel('Average ipos power recoded per session')
%%
hpc_fr = [];
acc_fr = [];
for i = 1:size(HPC.fr,1)
    for j = 1:size(HPC.fr,2)
        if all_valid(i,j)==true
            hpc_fr = cat(1, hpc_fr, HPC.fr{i,j});
            acc_fr = cat(1, acc_fr, ACC.fr{i,j});
        end
    end
end
hpc_fr = hpc_fr*60;
acc_fr = acc_fr*60;
% figure; histogram(acc_fr); hold on; histogram(hpc_fr)
figure(1915); clf; hold on
set(gcf, 'Position', [411   524   576   287])
bins = [0:.001:.07];
hh = histogram(hpc_fr, bins, 'Normalization', 'count', 'FaceColor', cyan);
ah = histogram(acc_fr, bins, 'Normalization', 'count', 'FaceColor', magenta);
% hh = histogram(hpc_fr, bins, 'Normalization', 'cdf', 'FaceColor', cyan);
% ah = histogram(acc_fr, bins, 'Normalization', 'cdf', 'FaceColor', magenta);
plot([nanmedian(hpc_fr), nanmedian(hpc_fr)], [-100 1600], '-', 'Color', cyan/1.5, 'LineWidth', 2)
plot([nanmedian(acc_fr), nanmedian(acc_fr)], [-100 1600], '--', 'Color', magenta/1.5, 'LineWidth', 2)
[p, h, stats] = ranksum(hpc_fr, acc_fr);
axis([-.005 .071 -100 1600])
% axis([-.005 .071 -.1 1.1])

hpc_sfeproom = [];
acc_sfeproom = [];
hpc_sfeparena = [];
acc_sfeparena = [];
gv = all_valid;%entrperMin<1 & (isTR);
for i = 1:size(HPC.fr,2)
    for j = 1:size(HPC.fr,1)
        if gv(j,i)==true
            hpc_sfeproom = cat(3, hpc_sfeproom, HPC.room.sfep.map{j,i});
            acc_sfeproom = cat(3, acc_sfeproom, ACC.room.sfep.map{j,i});
            hpc_sfeparena = cat(3, hpc_sfeparena, HPC.arena.sfep.map{j,i});
            acc_sfeparena = cat(3, acc_sfeparena, ACC.arena.sfep.map{j,i});
%             figure; imagesc(squeeze(HPC.room.sfep.map{j,i}))
        end
    end
end
figure(1916); clf; lims = [0 .8];
subplot(2,2,1); i = squeeze(nanmean(hpc_sfeproom,3)); ii = imagesc(i, lims);  ii.AlphaData = double(~isnan(i)); 
set(gca, 'YDir', 'normal'); axis off image; colorbar; title('HPC room')
subplot(2,2,2); i = squeeze(nanmean(acc_sfeproom,3)); ii = imagesc(i, lims);  ii.AlphaData = double(~isnan(i)); 
set(gca, 'YDir', 'normal'); axis off image; colorbar; title('ACC room')
subplot(2,2,3); i = squeeze(nanmean(hpc_sfeparena,3)); ii = imagesc(i, lims);  ii.AlphaData = double(~isnan(i)); 
set(gca, 'YDir', 'normal'); axis off image; colorbar; title('HPC arena')
subplot(2,2,4); i = squeeze(nanmean(acc_sfeparena,3)); ii = imagesc(i, lims);  ii.AlphaData = double(~isnan(i)); 
set(gca, 'YDir', 'normal'); axis off image; colorbar; title('ACC room')
colormap magma

figure(1917); clf; hold on 
histogram(hpc_sfeproom(:), [0:.025:1], 'Normalization', 'probability', 'FaceColor', cyan, 'EdgeColor', 'r', 'FaceAlpha', .3)
histogram(hpc_sfeparena(:), [0:.025:1], 'Normalization', 'probability', 'FaceColor', cyan, 'EdgeColor', 'b', 'FaceAlpha', .3)
histogram(acc_sfeproom(:), [0:.025:1], 'Normalization', 'probability', 'FaceColor', magenta, 'EdgeColor', 'r', 'FaceAlpha', .3)
histogram(acc_sfeparena(:), [0:.025:1], 'Normalization', 'probability', 'FaceColor', magenta, 'EdgeColor', 'b', 'FaceAlpha', .3)
% HPCACC24500 sessns 9 16 20 25
% Hipp18240 sessns 21 27
%%
label = anova_label;
ishpc = label==1 | label==2;
isroom = label==1 | label==3;
effects = [];

% [p,tbl,stats] = kruskalwallis(anova_pfield_info_av, anova_pfield_info_av_label)
[p,tbl,stats] = kruskalwallis(anova_pfield_info_av, anova_label); %#ok<ASGLU>
c = multcompare(stats);
effects.region.spatial_info = ranksum(anova_pfield_info_av(ishpc), anova_pfield_info_av(~ishpc));
effects.frame.spatial_info = ranksum(anova_pfield_info_av(isroom), anova_pfield_info_av(~isroom));

[p,tbl,stats] = kruskalwallis(anova_percent_place, anova_label); %#ok<ASGLU>
% c = multcompare(stats)
effects.region.perc_place = ranksum(anova_percent_place(ishpc), anova_percent_place(~ishpc));
effects.frame.perc_place = ranksum(anova_percent_place(isroom), anova_percent_place(~isroom));



% [p,tbl,stats] = kruskalwallis(anova_pfield_stability_av, anova_label); %#ok<ASGLU>
% c = multcompare(stats)
effects.region.field_corr = ranksum(anova_pfield_stability_av(ishpc), anova_pfield_stability_av(~ishpc));
effects.frame.field_corr = ranksum(anova_pfield_stability_av(isroom), anova_pfield_stability_av(~isroom));


% [p,tbl,stats] = kruskalwallis(anova_percent_stable, anova_label); %#ok<ASGLU>
% c = multcompare(stats)
effects.region.perc_stable = ranksum(anova_percent_stable(ishpc), anova_percent_stable(~ishpc));
effects.frame.perc_stable = ranksum(anova_percent_stable(isroom), anova_percent_stable(~isroom));


% [p,tbl,stats] = kruskalwallis(anova_BaysPerf, anova_label);%, 'off'); %#ok<ASGLU>
% c = multcompare(stats)
effects.region.Bayes = ranksum(anova_BaysPerf(ishpc), anova_BaysPerf(~ishpc));
effects.frame.Bayes = ranksum(anova_BaysPerf(isroom), anova_BaysPerf(~isroom));


[p,tbl,stats] = kruskalwallis(anova_SVMPerf, anova_label); %#ok<ASGLU>
c = multcompare(stats);
effects.region.SVM = ranksum(anova_SVMPerf(ishpc), anova_SVMPerf(~ishpc));
effects.frame.SVM = ranksum(anova_SVMPerf(isroom), anova_SVMPerf(~isroom));


[p,tbl,stats] = kruskalwallis(anova_pfield_coh_av, anova_label); %#ok<ASGLU>
c = multcompare(stats);
effects.region.field_coh = ranksum(anova_pfield_coh_av(ishpc), anova_pfield_coh_av(~ishpc));
effects.frame.field_coh = ranksum(anova_pfield_coh_av(isroom), anova_pfield_coh_av(~isroom));
effects.region_frame.field_coh = ranksum(anova_pfield_coh_av(ishpc&~isroom), anova_pfield_coh_av(~ishpc&~isroom));



[p,tbl,stats] = kruskalwallis(anova_ipos_roompref, anova_label); %#ok<ASGLU>
c = multcompare(stats);
effects.frame.sfep_pref  = ranksum(anova_ipos_roompref(isroom), anova_ipos_roompref(~isroom));

% [~,tbl,stats] = kruskalwallis(anova_sfep.coherence, anova_label); %#ok<ASGLU>
% c = multcompare(stats);
% effects.frame.sfep_coh  = ranksum(anova_sfep.coherence(isroom), anova_sfep.coherence(~isroom));
% effects.region.sfep_coh = ranksum(anova_sfep.coherence(ishpc), anova_sfep.coherence(~ishpc));
% 
% 
% [p,tbl,stats] = kruskalwallis(anova_sfep.bits, anova_label); %#ok<ASGLU>
% c = multcompare(stats);
% effects.frame.sfep_bits  = ranksum(anova_sfep.bits(isroom), anova_sfep.bits(~isroom));
% effects.region.sfep_bits = ranksum(anova_sfep.bits(ishpc), anova_sfep.bits(~ishpc));

% % [p,tbl,stats] = kruskalwallis(anova_sfep.bits-anova_sfep.bits_infer, anova_label); %#ok<ASGLU>
% % c = multcompare(stats);

[p,tbl,stats] = kruskalwallis(anova_sfep.mapcorr, anova_label); %#ok<ASGLU>
c = multcompare(stats);
effects.frame.mapcorr = ranksum(anova_sfep.mapcorr(isroom), anova_sfep.mapcorr(~isroom));
effects.regionxframe.mapcorr = ranksum(anova_sfep.mapcorr(ishpc&isroom), anova_sfep.mapcorr(~ishpc&isroom));
% signrank(anova_sfep.nulldiff(ishpc&isroom))
% signrank(anova_sfep.nulldiff(ishpc&~isroom))
% signrank(anova_sfep.nulldiff(~ishpc&isroom))
% signrank(anova_sfep.nulldiff(~ishpc&~isroom))

[p,tbl,stats] = kruskalwallis(anova_ipos_roompref, anova_label); %#ok<ASGLU>
c = multcompare(stats);

[p,tbl,stats] = kruskalwallis([HPC.ipos_runsz(all_valid); ACC.ipos_runsz(all_valid)], [HPC.ipos_runsz(all_valid)*0; ACC.ipos_runsz(all_valid)*0+1]); %#ok<ASGLU>
c = multcompare(stats);


% effects.region.field_coh = ranksum(anova_pfield_coh_av(ishpc), anova_pfield_coh_av(~ishpc));
% effects.frame.field_coh = ranksum(anova_pfield_coh_av(isroom), anova_pfield_coh_av(~isroom));
% effects.region_frame.field_coh = ranksum(anova_pfield_coh_av(ishpc&~isroom), anova_pfield_coh_av(~ishpc&~isroom));
acc_sfepdiff = (ACC.room.ipos_roompref - ACC.arena.ipos_roompref);
acc_sfepdiff(~all_valid) = NaN;
hpc_sfepdiff = (HPC.room.ipos_roompref - HPC.arena.ipos_roompref);
hpc_sfepdiff(~all_valid) = NaN;
figure(1120); clf; hold on
marks = {'d' 'o' 'x'};
clr = light_colormap(jet(numAnimals), 2);
plot([-1 1], [-1 1], 'k')
ind = 0;
for i = 1:numAnimals
    if any(acc_sfepdiff(:,i).*hpc_sfepdiff(:,i))
        ind = ind+1;
        scatter(acc_sfepdiff(:,i), hpc_sfepdiff(:,i), 50, 'Marker', marks{ind},...
            'MarkerFaceColor', clr(i,:)/1.25, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .65)
    end
    
end
ind = 0;
for i = 1:numAnimals
    if any(acc_sfepdiff(:,i).*hpc_sfepdiff(:,i))
        ind = ind+1;
        scatter(nanmedian(acc_sfepdiff(:,i)), nanmedian(hpc_sfepdiff(:,i)), 100, 'Marker', marks{ind},...
            'MarkerFaceColor', clr(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .65)
    end
end
xlabel('ACC (room-arena)')
ylabel('HPC (room-arena)')
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
h = histcounts2(HPC.arena.allcorr, HPC.room.allcorr, bins, bins, 'Normalization', 'count');
a = histcounts2(ACC.arena.allcorr, ACC.room.allcorr, bins, bins, 'Normalization', 'count');
maxval = 10*ceil(max([a(:);h(:)])/10);
subplot(2,2,1); hold on
imagesc(h, [0 maxval]); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', cyan)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis([-5 length(bins)+5 -5 length(bins)+5])
axis square
cc = colorbar; cc.Ticks = [0 maxval/2 maxval];
title('CA1')
xlabel('Pfield corr (Room)')
ylabel('Pfield corr (Arena)')
% colormap viridis
colormap magma

subplot(2,2,2); hold on
imagesc(a, [0 maxval]); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', magenta)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
title('ACC')
xlabel('Pfield corr (Room)')
ylabel('Pfield corr (Arena)')
axis square
cc = colorbar; cc.Ticks = [0 maxval/2 maxval];
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
h = histcounts2(HPC.arena.allbits, HPC.room.allbits, bins, bins, 'Normalization', 'count');
a = histcounts2(ACC.arena.allbits, ACC.room.allbits, bins, bins, 'Normalization', 'count');
maxval = 10*ceil(max([h(:)])/10);
subplot(2,2,1); hold on
imagesc(h, [0 maxval]); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', cyan)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis([-5 length(bins)+5 -5 length(bins)+5])
axis square
cc = colorbar; cc.Ticks = [0 maxval/2 maxval];
title('CA1')
xlabel('Pfield bits (Room)')
ylabel('Pfield bits (Arena)')
% colormap plasma
colormap magma

subplot(2,2,2); hold on
maxval = 10*ceil(max([a(:)])/10);
imagesc(a, [0 maxval]); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', magenta)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis square
cc = colorbar; cc.Ticks = [0 maxval/2 maxval];
title('ACC')
xlabel('Pfield bits (Room)')
ylabel('Pfield bits (Arena)')
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
figure(1260); clf;
imagesc(rand(10,10), [0 1]); colormap redblue; colorbar

figure(1261); clf;
bins = [0:.025:1.0];
xs = [1:20:length(bins)];
h = histcounts2(HPC.arena.allcoh, HPC.room.allcoh, bins, bins, 'Normalization', 'count');
a = histcounts2(ACC.arena.allcoh, ACC.room.allcoh, bins, bins, 'Normalization', 'count');
maxval = 10*ceil(max([a(:);h(:)])/10);
subplot(2,2,1); hold on
imagesc(h, [0 maxval]); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', cyan)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis([-5 length(bins)+5 -5 length(bins)+5])
axis square
cc = colorbar; cc.Ticks = [0 maxval/2 maxval];
title('CA1')
xlabel('Pfield coherence (Room)')
ylabel('Pfield coherence (Arena)')
colormap magma

subplot(2,2,2); hold on
imagesc(a, [0 maxval]); plot([1,length(bins)-1], [1,length(bins)-1], 'Color', magenta)
% p.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XtickLabel', bins(xs), 'Color', 'none') % '0.2670    0.0049    0.3294')
set(gca, 'YTick', xs, 'YtickLabel', bins(xs))
axis square
cc = colorbar; cc.Ticks = [0 maxval/2 maxval];
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
%% CORR AND SPATIAL, unused
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
%% PLOT DECODING BY TIME
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
        plot(xv, ya_av, acc_plot_style{:})
        plot(xv, yh_av, hpc_plot_style{:})
        set(gca, 'YTick', 0:.4:1.6, 'XTick', [0:5:25])
        axis([-1 12 -.1 1.6])
%         sub_err_bars(xa, ya, [], .2, 'std', acc_plot_style, acc_scatter_style)
%         sub_err_bars(xh, yh, [], .2, 'std', hpc_plot_style, hpc_scatter_style)
        set(gca, 'YTick', 0:.4:1.6, 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('Bits/spike')
        
        ya_av = tempA.pfield_stability_av(valid,animalLoop);
        yh_av = tempH.pfield_stability_av(valid,animalLoop);
        [ya, xa] = gb_cell2vec(tempA.pfield_stability(:,animalLoop), valid, xv);
        [yh, xh] = gb_cell2vec(tempH.pfield_stability(:,animalLoop), valid, xv);
        
        subplot(2,4,5); hold on
        axis([-1 12 -.1 1.1])
        plot(xv, ya_av, acc_plot_style{:})
        plot(xv, yh_av, hpc_plot_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        axis([-1 12 -.1 1.1])
%         sub_err_bars(xh, yh, [], .2, 'std', hpc_plot_style, hpc_scatter_style)
%         sub_err_bars(xa, ya, [], .2, 'std', acc_plot_style, acc_scatter_style)
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('Within Sess Corr')
        
        ya = tempA.percent_place(valid,animalLoop); xa = xv;
        yh = tempH.percent_place(valid,animalLoop); xh = xv;
        subplot(2,4,3); hold on
        axis([-1 12 -.1 1.1])
        plot(xa, ya, acc_plot_style{:})
%         scatter(xa, ya, acc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('% Place cells')
        
        %         subplot(2,4,4); hold on
        axis([-1 12 -.1 1.1])
        plot(xa, yh, hpc_plot_style{:})
%         scatter(xh, yh, hpc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        xlabel('Recording Session')
        ylabel('% Place cells')
        
        ya = tempA.percent_stable(valid,animalLoop); xa = xv;
        yh = tempH.percent_stable(valid,animalLoop); xh = xv;
        subplot(2,4,7); hold on
        axis([-1 12 -.1 1.1])
        plot(xa, ya, acc_plot_style{:})
%         scatter(xa, ya, acc_scatter_style{:})
        set(gca, 'YTick', 0:.2:1, 'XTick', [0:5:25])
        xlabel('Entrances/Minute')
        ylabel('% Stable cells')
        
        %         subplot(2,4,8); hold on
        axis([-1 12 -.1 1.1])
        plot(xa, yh, hpc_plot_style{:})
%         scatter(xh, yh, hpc_scatter_style{:})
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
% HPCACC24500 sessns 9 16 20 25
% Hipp18240 sessns 21 27
% accfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files/2023_06_15_H14_25_45_TR13_@placecells_ACC_miniscope2.mat';
% hpcfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files/2023_06_15_H14_25_45_TR13_@placecells_HPC_miniscope1.mat';
% accfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_15_H11_24_36_TR12_@placecells_ACC_miniscope2.mat';
% hpcfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_15_H11_24_36_TR12_@placecells_HPC_miniscope1.mat';
accfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_28_H17_31_09_CON20_@placecells_ACC_miniscope2.mat';
hpcfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_28_H17_31_09_CON20_@placecells_HPC_miniscope1.mat';
hpc_ms = load(hpcfile);
acc_ms = load(accfile);

r = hpc_ms.ms.room;
a = hpc_ms.ms.arena;

%%
figure(1005); clf; 
subplot_tight(3,1,1, [.05 .1])
% sp = hpc_ms.ms.room.svm_decoding.spks_bin(1:50,:)>0;
hsp = hpc_ms.ms.neuron.C+hpc_ms.ms.neuron.YrA;
hsp = normalize_rows(hsp);
hsp = bin_spks_average(hsp(:, 1:10000), 5, false);
t = bin_spks_average(hpc_ms.ms.timestamps(1:10000)'./1000, 5, false);
t1 = find(t>=60,1); t2 = find(t>=120,1);
% hsp = 1-cat(3, sp, sp, sp);
stacked_traces(hsp(1:15,:), 2.5, {'Color', cyan/1.5, 'LineWidth', 2}); axis tight
set(gca, 'XTick', [t1 t2], 'XTickLabel', [60 120], 'YTick', [2 16], 'YTickLabel', [1 15])
subplot_tight(3,1,2, [.05 .1])
asp = acc_ms.ms.neuron.C+acc_ms.ms.neuron.YrA;
asp = normalize_rows(asp);
asp = bin_spks_average(asp(:, 1:10000), 5, false);
t = bin_spks_average(acc_ms.ms.timestamps(1:10000)'./1000, 5, false);
t1 = find(t>=60,1); t2 = find(t>=120,1);
% hsp = 1-cat(3, sp, sp, sp);
stacked_traces(asp(1:15,:), 2.5, {'Color', magenta/1.5, 'LineWidth', 2}); axis tight
set(gca, 'XTick', [t1 t2], 'XTickLabel', [60 120], 'YTick', [2 16], 'YTickLabel', [1 15])

subplot_tight(3,1,3, [.05 .1])
plot(acc_ms.ms.room.svm_decoding.t, nanmean(acc_ms.ms.room.svm_decoding.spks_bin,1), 'Color', magenta/1.2); hold on
plot(hpc_ms.ms.room.svm_decoding.t, nanmean(hpc_ms.ms.room.svm_decoding.spks_bin,1), 'Color', cyan/1.2); 
ylabel('% Cells Active'); xlabel('Time (sec)')
%%
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
%%
hpc_ms.params.skip_ensemble = true;
h_time = hpc_ms.ms.room.svm_decoding.t;
a_time = acc_ms.ms.room.svm_decoding.t;

% [ang_dist, inside_ang_ref] = angular_distance(ang1, ang_ref, bleed_size);

h_ipos = abs(hpc_ms.ms.room.momentary_pos_info) - abs(hpc_ms.ms.arena.momentary_pos_info);
a_ipos = abs(acc_ms.ms.room.momentary_pos_info) - abs(acc_ms.ms.arena.momentary_pos_info);
h_ipos_infer = hpc_ms.ms.room.inferred_ipos - hpc_ms.ms.arena.inferred_ipos;
a_ipos_infer = acc_ms.ms.room.inferred_ipos - acc_ms.ms.arena.inferred_ipos;
h_spks = hpc_ms.ms.room.svm_decoding.spks_bin; 
a_spks = acc_ms.ms.room.svm_decoding.spks_bin;
goodh = ~any(isnan([h_ipos]),1);
gooda = ~any(isnan([a_ipos]),1);
goodhi = ~any(isnan([h_ipos_infer]),1);
goodai = ~any(isnan([a_ipos_infer]),1);
if length(h_time)<length(a_time)
    t = a_time(goodh);
    h_ipos_match = NaN(size(h_ipos, 1), sum(gooda));
    h_ipos_match_infer = NaN(size(h_ipos_infer, 1), sum(goodai));
    h_spk_match = NaN(size(h_ipos, 1), sum(gooda));
    for i = 1:size(h_ipos)
        h_ipos_match(i,:) = interp1(h_time(goodh), h_ipos(i, goodh), t, 'linear')';
        h_ipos_match_infer(i,:) = interp1(h_time(goodhi), h_ipos_infer(i, goodhi), a_time(goodai), 'linear')';
        h_spk_match(i,:)  = interp1(h_time(goodh), h_spks(i, goodh), t, 'linear')';
    end
    a_ipos_match = a_ipos(:,gooda);
    a_ipos_match_infer = a_ipos_infer(:,goodai);
    a_spk_match = a_spks(:,gooda);
    
    h_ipos_mean = interp1(h_time(goodh), nanmean(h_ipos(:, goodh),1), t, 'linear')';
    a_ipos_mean = nanmean(a_ipos(:, gooda),1);
else
    t = h_time(goodh);
    a_ipos_match = NaN(size(a_ipos, 1), sum(goodh));
    a_ipos_match_infer = NaN(size(a_ipos_infer, 1), sum(goodhi));
    a_spk_match = NaN(size(a_ipos, 1), sum(goodh));
    for i = 1:size(a_ipos)
        a_ipos_match(i,:) = interp1(a_time(gooda), a_ipos(i, gooda), t, 'linear')';
        a_ipos_match_infer(i,:) = interp1(a_time(goodai), a_ipos_infer(i, goodai), h_time(goodhi), 'linear')';
        a_spk_match(i,:)  = interp1(a_time(gooda), a_spks(i, gooda), t, 'linear')';
    end
    h_ipos_match = h_ipos(:,goodh);
    h_ipos_match_infer = h_ipos_infer(:,goodhi);
    h_spk_match = h_spks(:,goodh);
    
    a_ipos_mean = interp1(a_time(gooda), nanmean(a_ipos(:, gooda),1), h_time(goodh), 'linear')';
    h_ipos_mean = nanmean(h_ipos(:, goodh),1);
end
% [h_ipos] = Fenton_ipos(hpc_ms.ms, .25, 'arena', hpc_ms.params);

goodinds = ~any(isnan([h_ipos_match; a_ipos_match]),1);
goodindsi = ~any(isnan([h_ipos_match_infer; a_ipos_match_infer]),1);
% h_ipos = zscore( h_ipos_match(:, goodinds) );
% a_ipos = zscore( a_ipos_match(:, goodinds) );
h_ipos = ( h_ipos_match(:, goodinds) );
a_ipos = ( a_ipos_match(:, goodinds) );
h_ipos_infer = ( h_ipos_match_infer(:, goodindsi) );
a_ipos_infer = ( a_ipos_match_infer(:, goodindsi) );
[hpc_acc_iposcorr, hpc_acc_iposcorrpval] = nancorr( nanmean(h_ipos,1), nanmean(a_ipos,1) );
[hpc_acc_iposcorr_infer, hpc_acc_iposcorrpval_infer] = nancorr( nanmean(h_ipos_infer,1), nanmean(a_ipos_infer,1) );
h_ipos = (bin_spks_average(h_ipos, 5, false));
a_ipos = (bin_spks_average(a_ipos, 5, false));
h_ipos_infer = (bin_spks_average(h_ipos_infer, 5, false));
a_ipos_infer = (bin_spks_average(a_ipos_infer, 5, false));

h_spks = ( h_spk_match(:, goodinds) );
a_spks = ( a_spk_match(:, goodinds) );

hsp = normalize_rows(bin_spks_average(h_spks, 5, false));
asp = normalize_rows(bin_spks_average(a_spks, 5, false));
tsub = bin_spks_average(t', 5, false);
t1 = find(tsub>=60,1); t2 = find(tsub>=120,1);

% [c, p] = xcorr(mean(h_ipos,1)', mean(a_ipos,1)');
% [c, p] = corr(mean(h_ipos,1)', mean(a_ipos,1)', 'type', 'Kendall');
% [c, p] = corr([h_ipos; a_ipos(1:10,:)*0; a_ipos]');
allipos = [h_ipos; a_ipos];
alliposshow = [h_ipos; NaN(10, size(h_ipos,2)); a_ipos];
%%
figure(1117); clf; set(gcf, 'Color', 'w', 'Position', [327   287   600   500])
subplot_tight(2,1,1);
imagesc(h_ipos(1:50, :), [-.5 .5]); colormap redblue; colorbar;
set(gca, 'XTick', [t1 t2], 'XTickLabel', [60 120], 'YTick', [.5 50.5], 'YTickLabel', [1 50])
subplot_tight(2,1,2);
imagesc(a_ipos(1:50, :), [-.5 .5]); colormap redblue; colorbar;
set(gca, 'XTick', [t1 t2], 'XTickLabel', [60 120], 'YTick', [.5 50.5], 'YTickLabel', [1 50])

figure(1118); clf; set(gcf, 'Color', 'w', 'Position', [527   287   600   500])
subplot_tight(2,1,1);
imagesc(hsp(1:50, :), [0 1]); colormap viridis; colorbar;
set(gca, 'XTick', [t1 t2], 'XTickLabel', [60 120], 'YTick', [.5 50.5], 'YTickLabel', [1 50])
subplot_tight(2,1,2);
imagesc(asp(1:50, :), [0 1]); colormap magma; colorbar;
set(gca, 'XTick', [t1 t2], 'XTickLabel', [60 120], 'YTick', [.5 50.5], 'YTickLabel', [1 50])

figure(1119); clf; set(gcf, 'Color', 'w', 'Position', [727   287   600   500])
subplot(2,1,1); hold on;
plot(nanmean(hsp,1), 'Color', cyan/1.2);
plot(nanmean(asp,1), 'Color', magenta/1.2);
axis tight
ylim([0 .15])
subplot(2,1,2); hold on;
plot([0 size(h_ipos,2)], [0 0], 'Color', 'k');
plot(nanmean(h_ipos,1), 'Color', cyan/1.2);
plot(nanmean(a_ipos,1), 'Color', magenta/1.2);
axis tight
ylim([-.1 .1])

figure; 
subplot(1,2,1);
% pie([HPC.room.ipos_roompref(exp_num(:,1)==20,1), HPC.arena.ipos_roompref(exp_num(:,1)==20,1)], {'room' 'arena'})
pie(100*[HPC.room.ipos_roompref(exp_num(:,1)==20,1), HPC.arena.ipos_roompref(exp_num(:,1)==20,1)])
title('HPC')
subplot(1,2,2);
% pie([ACC.room.ipos_roompref(exp_num(:,1)==20,1), ACC.arena.ipos_roompref(exp_num(:,1)==20,1)], {'room' 'arena'})
pie(100*[ACC.room.ipos_roompref(exp_num(:,1)==20,1), ACC.arena.ipos_roompref(exp_num(:,1)==20,1)])
title('ACC')


figure(1120); clf; set(gcf, 'Color', 'w', 'Position', [927   287   600   500])
subplot(2,1,1); hold on;
ipos_patch(h_ipos);
plot([1 size(h_ipos,2)], [0 0], 'k');
plot(nanmean(h_ipos,1), 'Color', cyan);
axis tight
ylim([-.1 .1])
subplot(2,1,2); hold on;
ipos_patch(a_ipos);
plot([1 size(a_ipos,2)], [0 0], 'k');
plot(nanmean(a_ipos,1), 'Color', magenta);
axis tight
ylim([-.1 .1])

figure(1121); clf; set(gcf, 'Color', 'w', 'Position', [1127   287   600   500])
subplot(3,1,1); hold on;
plot([0 size(h_ipos_infer,2)], [0 0], 'Color', 'k');
plot(nanmean(h_ipos_infer,1), 'Color', cyan/1.2);
plot(nanmean(a_ipos_infer,1), 'Color', magenta/1.2);
axis tight
ylim([-.1 .1])

subplot(3,1,2); hold on;
ipos_patch(h_ipos_infer);
plot([1 size(h_ipos_infer,2)], [0 0], 'k');
plot(nanmean(h_ipos_infer,1), 'Color', cyan);
axis tight
ylim([-.1 .1])
subplot(3,1,3); hold on;
ipos_patch(a_ipos_infer);
plot([1 size(a_ipos_infer,2)], [0 0], 'k');
plot(nanmean(a_ipos_infer,1), 'Color', magenta);
axis tight
ylim([-.1 .1])


%% trying to look at switches from one frame to another
h_ipos = ( h_ipos_match(:, goodinds) );
a_ipos = ( a_ipos_match(:, goodinds) );
mh = nanmean(h_ipos,1);
mhs = conv(mh, ones(5,1)./5, 'same');
% ma = mh*2;
ma = nanmean(a_ipos,1);
mas = conv(ma, ones(5,1)./5, 'same');
ssig = [(mhs(1:end-1)<0 & mhs(2:end)>0), 0];
[posa, posma] = gb_PETH(ma, ssig, 10, 10);
ssig = [(mhs(1:end-1)>0 & mhs(2:end)<0), 0];
[nega, negma] = gb_PETH(ma, ssig, 10, 10);
ssig = [(mas(1:end-1)<0 & mas(2:end)>0), 0];
[posh, posmh] = gb_PETH(mh, ssig, 10, 10);
ssig = [(mas(1:end-1)>0 & mas(2:end)<0), 0];
[negh, negmh] = gb_PETH(mh, ssig, 10, 10);
figure(120); clf;
subplot(1,2,1); hold on;
lims = [min([posma; posmh; negma; negmh])*1.2 max([posma; posmh; negma; negmh])*1.2];
plot([0,0], lims, 'k:')
% plot([-10:10], posa, 'Color', magenta);
% plot([-10:10], posh, 'Color', cyan);
shadedErrorBar([-10:10], posmh, nanstd(posh,[],2)./sqrt(size(posh,2)-1), 'lineProps',{'Color', cyan});
shadedErrorBar([-10:10], posma, nanstd(posa,[],2)./sqrt(size(posa,2)-1), 'lineProps',{'Color', magenta});
% plot([-10:10], posma, 'Color', magenta/1.2);
% plot([-10:10], posmh, 'Color', cyan/1.2);
axis tight

subplot(1,2,2); hold on;
plot([0,0], lims, 'k:')
shadedErrorBar([-10:10], negmh, nanstd(negh,[],2)./sqrt(size(negh,2)-1), 'lineProps',{'Color', cyan});
shadedErrorBar([-10:10], negma, nanstd(nega,[],2)./sqrt(size(nega,2)-1), 'lineProps',{'Color', magenta});
% plot([-10:10], negma, 'Color', magenta);
% plot([-10:10], negmh, 'Color', cyan);
axis tight
%% SFEP MAP PLOT
% hr_sfep = imrotate(HPC.room.sfep.map{(exp_num(:,1)==20), 1}, -180);
% ha_sfep = imrotate(HPC.arena.sfep.map{(exp_num(:,1)==20), 1}, -180);
% ar_sfep = imrotate(ACC.room.sfep.map{(exp_num(:,1)==20), 1}, -180);
% aa_sfep = imrotate(ACC.arena.sfep.map{(exp_num(:,1)==20), 1}, -180);
[hsfep] = SFEP_map_quantification(hpc_ms.ms, sfep_smoothingsize, false, false);
[asfep] = SFEP_map_quantification(acc_ms.ms, sfep_smoothingsize, false, false);

% hr_sfep = normalize_matrix(hsfep.room_map_smooth);
% ha_sfep = normalize_matrix(hsfep.arena_map_smooth);
% ar_sfep = normalize_matrix(asfep.room_map_smooth);
% aa_sfep = normalize_matrix(asfep.arena_map_smooth);
% hr_sfep_i = normalize_matrix(hsfep.infer.room_map_smooth);%  imrotate(HPC.room.sfep.infer_map{(exp_num(:,1)==20), 1}, -180);
% ha_sfep_i = normalize_matrix(hsfep.infer.arena_map_smooth);%  imrotate(HPC.room.sfep.infer_map{(exp_num(:,1)==20), 1}, -180);
% ar_sfep_i = normalize_matrix(asfep.infer.room_map_smooth);%  imrotate(HPC.room.sfep.infer_map{(exp_num(:,1)==20), 1}, -180);
% aa_sfep_i = normalize_matrix(asfep.infer.arena_map_smooth);%  imrotate(HPC.room.sfep.infer_map{(exp_num(:,1)==20), 1}, -180);
% hr_sfep = imrotate(hsfep.room_map, -180);
% ha_sfep = imrotate(hsfep.arena_map, -180);
% ar_sfep = imrotate(asfep.room_map, -180);
% aa_sfep = imrotate(asfep.arena_map, -180);
% hr_sfep_i = imrotate(hsfep.infer.room_map, -180);
% ha_sfep_i = imrotate(hsfep.infer.arena_map, -180);
% ar_sfep_i = imrotate(asfep.infer.room_map, -180);
% aa_sfep_i = imrotate(asfep.infer.arena_map, -180);
hr_sfep = imrotate(hsfep.room_map_smooth, -180);
ha_sfep = imrotate(hsfep.arena_map_smooth, -180);
ar_sfep = imrotate(asfep.room_map_smooth, -180);
aa_sfep = imrotate(asfep.arena_map_smooth, -180);
hr_sfep_i = imrotate(hsfep.infer.room_map_smooth, -180);
ha_sfep_i = imrotate(hsfep.infer.arena_map_smooth, -180);
ar_sfep_i = imrotate(asfep.infer.room_map_smooth, -180);
aa_sfep_i = imrotate(asfep.infer.arena_map_smooth, -180);


figure(10025); clf;
colormap magma
subplot(2,2,1); plot(r.x,r.y,'k'); axis([-44 44 -44 44]);  axis equal off
subplot(2,2,2); plot(a.x,a.y,'k'); axis([-44 44 -44 44]);  axis equal off
subplot(2,2,3); i = imrotate(r.vmap, -180); ii = imagesc(i, [0 20]); ii.AlphaData = double(~isnan(i)); set(gca,'YDir', 'normal'); axis image off; colorbar
subplot(2,2,4); i = imrotate(a.vmap, -180); ii = imagesc(i, [0 10]); ii.AlphaData = double(~isnan(i)); set(gca,'YDir', 'normal'); axis image off; colorbar

figure(10026); clf;
colormap magma
subplot(2,3,1); ii = imagesc(hr_sfep, [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(hr_sfep)); axis image off; colorbar
title('HPC')
subplot(2,3,2); ii = imagesc(hr_sfep_i, [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(hr_sfep_i)); axis image off; colorbar
subplot(2,3,3); ii = imagesc(normalize_matrix(hr_sfep - hr_sfep_i), [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(hr_sfep - hr_sfep_i)); axis image off; colorbar
subplot(2,3,4); ii = imagesc(ha_sfep, [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(ha_sfep)); axis image off; colorbar
subplot(2,3,5); ii = imagesc(ha_sfep_i, [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(ha_sfep_i)); axis image off; colorbar
subplot(2,3,6); ii = imagesc(normalize_matrix(ha_sfep - ha_sfep_i), [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(ha_sfep - ha_sfep_i)); axis image off; colorbar

figure(10027); clf;
colormap magma
subplot(2,3,1); ii = imagesc(ar_sfep, [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(ar_sfep)); axis image off; colorbar
title('ACC')
subplot(2,3,2); ii = imagesc(ar_sfep_i, [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(ar_sfep_i)); axis image off; colorbar
subplot(2,3,3); ii = imagesc(normalize_matrix(ar_sfep - ar_sfep_i), [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(ar_sfep - ar_sfep_i)); axis image off; colorbar
subplot(2,3,4); ii = imagesc(aa_sfep, [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(aa_sfep)); axis image off; colorbar
subplot(2,3,5); ii = imagesc(aa_sfep_i, [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(aa_sfep_i)); axis image off; colorbar
subplot(2,3,6); ii = imagesc(normalize_matrix(aa_sfep - aa_sfep_i), [0 1]); set(gca,'YDir', 'normal'); ii.AlphaData = double(~isnan(aa_sfep - aa_sfep_i)); axis image off; colorbar



%%
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
% allipos = [h_ipos; a_ipos];
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
%%
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



