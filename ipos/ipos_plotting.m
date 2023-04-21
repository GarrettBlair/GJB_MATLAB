%% 
warning('off','stats:runstest:ValuesOmitted')

% integration_times = [.1 .25 .5 1 3 5 30];
integration_times = [2.5 5 10 20];
% integration_times = [.1 .2 1 3 5 10 30 60 120];
% integration_times = [1];
% animals = {'Acc20832', 'Acc19947'};
% animals = {'Acc20832', 'Acc19947'};
% sess = {'HAB2' 'TR0' 'TR3' 'TR6' 'RET0' 'NEW0' 'NEW4'};

animals = {'Acc20832', 'Acc19947', 'Hipp18240'};
% sess = {'TR0' 'TR3' 'TR6' 'WTR10' 'TR11' 'TR17' 'DRK15'};
sess = {'TR0' 'TR3' 'TR6'};
numAnimals = length(animals);

% sess = {'HAB0', 'HAB1', 'HAB2', 'HAB3', 'TR0', 'TR1', 'TR2', 'TR3', 'TR4', 'TR5', 'TR6', 'TR7'};
% sess = {'RET0', 'NEW0', 'NEW1', 'NEW2', 'NEW3', 'NEW4'};
% sess = {'HAB2' 'TR0' 'TR3' 'TR6'};
% nsess = length(sess);

shock_zone_center = pi/2; % typical room shock configuration
shock_zone_size = pi/6; % size in rad from center to edge
distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure

smoothing_interval = 3; % sec for averaging smoothing of ipos for correlation

nantemp         = NaN(numAnimals, 100, length(integration_times));
celltemp        = cell(numAnimals, 100, length(integration_times));

rval            = nantemp;
rsig            = nantemp;
bestlag         = nantemp;
rval_shift      = nantemp;
rsig_shift      = nantemp;
prop_cells_sig  = nantemp;
ipos_runs_sig   = nantemp;%celltemp;
ipos_runs_z     = nantemp;%celltemp;

ddir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\';
%%
saveVars = {'animal' 'thissess' 'integration_val' 'fname' 'params' 'room_momentary_pos'...
    'arena_momentary_pos' 'ms_room_temp' 'ms_arena_temp' 'roomrep_pf_roomframe'...
    'arenarep_pf_arenaframe' 'roomrep_pf_roomframe_smooth' 'arenarep_pf_arenaframe_smooth'...
    'ipos_ra' 'ipos_smooth' 'ipos_smooth_corr' 'sz_distance' 'room_ensemble_prob' 'arena_ensemble_prob'};
plotting_ipos = false;
for aLoop = 1:numAnimals
    animal = animals{aLoop};
    fdir = sprintf('%s%s\\processed_files\\', ddir, animal);
    savedir = sprintf('%sipos_files\\', fdir);
    cd(fdir)
    % use this for all
%     processedFiles = dir('*@placecells.mat');
%     nsess = length(processedFiles);
    % use this for specific sessions listed in {sess}
    nsess = length(sess);
    for sessLoop = 1:nsess
        %%
        thissess = sess{sessLoop};
        temp = dir(['*' thissess '*']);
        fname = temp.name;
        clearvars ms* params ipos* ms_room_temp ms_arena_temp
%         fname = processedFiles(sessLoop).name;
%         us = find(fname=='_', 2, 'last');
%         thissess = fname(us(1)+1:us(2)-1);

        fprintf('\n~~~ file - %s', fname)
        load(fname, 'ms')
        fprintf('\n~~~ file - %s', ms.fileName)
        for intLoop = 1:length(integration_times)
            integration_val = integration_times(intLoop);
            saveName = sprintf('%s_%s_int_%03.2fsec.mat', animal, thissess, integration_val);
            saveFullFile = sprintf('%s%s', savedir, saveName);
            %%
            for ijk = 5:length(saveVars)
                eval(sprintf('%s = [];', saveVars{ijk}));
                %                 fprintf('\n%s = [];', saveVars{ijk});
            end
            if true % isfile(saveFullFile)==false %%%%%%%%%%%%%%%%%%
                fprintf('\n~~~ deltaT = %3.2f', integration_val)
                fignum = aLoop*10000 + sessLoop*100 + intLoop*10;
                figname = sprintf('%s - %s ,  dT=%3.2f', animal, thissess, integration_val);
                fprintf('\n~~~ %s', figname)
                
                ii = aLoop; jj = sessLoop; kk = intLoop;
                
                %             params = [];
                params.pos_bins = [-45, -36:4:36, 45]; % in cm, x and y
                params.occupancy_thresh = .250; % sec
                params.ipos_int_time = integration_val; % .2; % seconds, binning time for computing momentary spatial information
                params.ipos_smoothing_interval = smoothing_interval;
                
                [room_momentary_pos, ms_room_temp, room_ensemble_prob]   = Fenton_ipos(ms, params.ipos_int_time, 'room', params);
                [arena_momentary_pos, ms_arena_temp, arena_ensemble_prob] = Fenton_ipos(ms, params.ipos_int_time, 'arena', params);
                ensemble_diff = abs(room_ensemble_prob) - abs(arena_ensemble_prob);
                ensemble_diff(isnan(ensemble_diff)) = 0;
                ipos_time = ms_room_temp.t; % ms.room.ipos_time
                r = abs(room_momentary_pos);
                r(isnan(r)) = 0;
                a = abs(arena_momentary_pos);
                a(isnan(a)) = 0;
                rm = mean(r);
                am = mean(a);
                smooth_5sec_len = round(5/integration_val);
                smooth_20sec_len = round(20/integration_val);
                ipos_ra = r-a;
                ipos_mean = nanmean(ipos_ra,1);
                
                %             ipos_smooth = bin_spks(ipos_ra, smooth_5sec_len, true);
                ipos_smooth = average_spks_time(ipos_ra, smoothing_interval, ms_room_temp.t, true, 'mean');
                ipos_smooth_corr = corr(ipos_smooth');
                ipos_smooth_corr(find(eye(size(ipos_smooth,1)))) = NaN;
                %%
                ipos_mean = nanmean(ipos_ra,1);
                ipos_smooth_mean = nanmean(ipos_smooth);
%                 roomrep_pf_roomframe = make_occupancymap_2D(ms_room_temp.x,  ms_room_temp.y,  ipos_mean, params.pos_bins, params.pos_bins);
%                 arenarep_pf_arenaframe = make_occupancymap_2D(ms_arena_temp.x, ms_arena_temp.y, ipos_mean, params.pos_bins, params.pos_bins);
                roomrep_pf_roomframe_smooth = make_occupancymap_2D(ms_room_temp.x, ms_room_temp.y, ipos_smooth_mean, params.pos_bins, params.pos_bins);
                arenarep_pf_arenaframe_smooth = make_occupancymap_2D(ms_arena_temp.x, ms_arena_temp.y, ipos_smooth_mean, params.pos_bins, params.pos_bins);
                
                ensemblerep_pf_roomframe = make_occupancymap_2D(ms_room_temp.x,  ms_room_temp.y,  ensemble_diff, params.pos_bins, params.pos_bins);
                ensemblerep_pf_arenaframe = make_occupancymap_2D(ms_arena_temp.x, ms_arena_temp.y, ensemble_diff, params.pos_bins, params.pos_bins);
                
                
                [nsegs, nframes] = size(ipos_ra);
                %             ipos_runs_sig{ii,jj,kk} = NaN(nsegs, 1);
                %             ipos_runs_z{ii,jj,kk} = NaN(nsegs, 1);
                pp = NaN(nsegs,1);
                zz = NaN(nsegs,1);
                for i = 1:nsegs
                    [~, p, h] = runstest(ipos_ra(i,:));
                    pp(i) = -log10(p);
                    zz(i) = h.z;
                    %                 ipos_runs_sig{ii,jj,kk}(i) = -log10(p);
                    %                 ipos_runs_z{ii,jj,kk}(i) = h.z;
                end
                %             figure(10); subplot_tight(2, 12, 12*(aLoop-1) + sessLoop)
                %             plot(pp)
                %             drawnow
                ipos_runs_sig(ii,jj,kk) = median(pp);
                ipos_runs_z(ii,jj,kk) = median(zz);
                prop_cells_sig(ii, jj, kk) = sum(pp > 4) / nsegs; % 2 = p=.0001 for -log10 scale
                
                x = ms.room.x; y = ms.room.y; t = ms.timestamps./1000;
                [th, rth] = cart2pol(x,y);
                % figure(1); clf; polarhistogram(th,24);
                d = shock_zone_center - th;
                d = abs(mod(d + pi, 2*pi) - pi);
                sz_distance = average_spks_time(d', params.ipos_int_time, ms.timestamps./1000, false, 'mean');
                
                % rplot = mean(ipos_ra,1); rplot(rplot<0) = NaN;
                % aplot = mean(ipos_ra,1); aplot(rplot>=0) = NaN;
                
                if any(ipos_mean)
                    [rval(ii, jj, kk), rsig(ii, jj, kk)] = corr(sz_distance', ipos_mean');
                    [xcv, lag] = xcorr(sz_distance', ipos_mean', smooth_20sec_len); % can look 10 sec before/after
                    bl = find(max(abs(xcv)) == abs(xcv));
                    bestlag(ii, jj, kk) = integration_val*lag(bl);
                    sz_distance_shift = circshift(sz_distance, -1*lag(bl));
                    [rval_shift(ii, jj, kk), rsig_shift(ii, jj, kk)] = corr(sz_distance_shift', ipos_mean');
                else
                    rval(ii, jj, kk)        = 0;
                    rsig(ii, jj, kk)        = NaN;
                    bestlag(ii, jj, kk)     = 0;
                    rval_shift(ii, jj, kk)  = 0;
                    rsig_shift(ii, jj, kk)  = NaN;
                end
                % dd = average_spks_time(d', params.ipos_int_time, ms.timestamps./1000, false, 'mean');
                %
                % disp(ms.fileName)
                if plotting_ipos==true
                    fprintf('\n\t%d : %d', fignum+1, fignum+3)
                    %%
%                     shk = ms.room.entranceTimes./1000;
%                     shk_ind = shk*NaN;
%                     for i = 1:length(shk)
%                         if shk(i) < ipos_time
%                             shk_ind(i) = 1;
%                         elseif shk(i) >= ipos_time(end)
%                             shk_ind(i) = length(ipos_time);
%                         else
%                             shk_ind(i) = find(ipos_time(1:end-1)<=shk(i) & ipos_time(2:end)>shk(i), 1);
%                         end
%                     end
                    figure(fignum+1); clf;
                    set(gcf, 'Name', figname);
                    % plot(rplot, 'r')
                    % plot(aplot, 'b')
                    subplot_tight(2,1,1); hold on
                    plot(normalize_rows(ipos_mean)-.5, 'k')
                    plot((ensemble_diff)-1.5, 'k')
%                     plot(ipos_smooth_mean, 'b')
                    plot(normalize_rows(sz_distance)+.5, 'r')
                    ylim([ -2 2])
                    subplot_tight(2,4, 5)
                    imagesc(ensemble_pf_roomframe); colorbar; axis image off
                    title('ensemblerep roomframe')
                    subplot_tight(2,4, 6)
                    imagesc(roomrep_pf_roomframe_smooth); colorbar; axis image off
                    title('roomrep roomframe smoothed')
                    subplot_tight(2,4, 7)
                    imagesc(ensemblerep_pf_arenaframe); colorbar; axis image off
                    title('ensemblerep arenaframe')
                    subplot_tight(2,4, 8)
                    imagesc(arenarep_pf_arenaframe_smooth); colorbar; axis image off
                    title('roomrep arenaframe smoothed')
                    % imagesc(normalize_rows(ipos_ra));
                    %%
                    figure(fignum+2); clf;
                    set(gcf, 'Name', figname);
                    imagesc(ipos_smooth_corr)
                    set(gca,'CLim',[-.2 .5])
                    
                    figure(fignum+3); clf;
                    set(gcf, 'Name', figname);
                    %                 imagesc(normalize_rows(ipos_smooth));
%                     imagesc(normalize_rows(ipos_ra));
%                     imagesc(ipos_ra, [-4 4]);
                    imagesc(ipos_smooth, [-4 4]);
                    %                 imagesc(ipos_smooth, [-4 4]);
                    % imagesc(ipos_ra, [-2 2]);
                    if false
                        sec2plot = 60;
                        t = ipos_time./1000 - ipos_time(1)./1000;
                        tind = mod(t, 60); tind = [0 find(tind(1:end-1)>=sec2plot/2 & tind(2:end)<sec2plot/2)];
                        tlabel = round([0 ceil(t(tind(2:end)))]./60);
                        set(gca, 'XTick', tind, 'XTickLabel', tlabel)
                        hold on;
                        yyaxis('right')
                        hold on
                        for i = 1:length(shk)
                            plot([shk(i), shk(i)] , [0 nsegs], 'r--', 'LineWidth', 1.5)
                            
                        end
                        plot(t, normalize_rows(sz_distance)*nsegs, 'k', 'LineWidth', 1)
                    else
                        hold on
                        
                        for i = 1:length(shk)
                            plot([shk_ind(i), shk_ind(i)] , [0 nsegs], 'r--', 'LineWidth', 1.5)
                            
                        end
                        plot(normalize_rows(sz_distance)*nsegs, 'k', 'LineWidth', 1)
                    end
                    drawnow
                end
                fprintf('\n\t\t ~~~ Saving:\n%s', saveName);
                save(saveFullFile, saveVars{:});
            else
                fprintf('\n\t\t ~~~ FILE SKIPPED:\n%s', saveName);
            end
        end
    end
end
fprintf('\n\t~~~ DONE ~~~')

%%
% rval            = nantemp;
% rsig            = nantemp;
% bestlag         = nantemp;
% rval_shift      = nantemp;
% rsig_shift      = nantemp;
% prop_cells_sig  = nantemp;
% ipos_runs_sig   = celltemp;
% ipos_runs_z     = celltemp;
if false
vars = {'\  rval' '\  prop_cells_sig' '\  ipos_runs_sig'};
nvar = length(vars);
figure(8); clf
ni = length(integration_times);
for intLoop = 1:ni
    a = squeeze(rval(:,:,intLoop));
    ind = nvar*(intLoop-1) + 1;
    subplot(ni, nvar, ind); ind = ind+1;
    plot(a')
    title(sprintf('%s  %.2f',vars{1}, integration_times(intLoop)))
    
    a = squeeze(prop_cells_sig(:,:,intLoop));
    subplot(ni, nvar, ind); ind = ind+1;
    plot(a')
    title(sprintf('%s  %.2f',vars{2}, integration_times(intLoop)))
   
    a = squeeze(ipos_runs_sig(:,:,intLoop));
    subplot(ni, nvar, ind); ind = ind+1;
    plot(a')
    title(sprintf('%s  %.2f',vars{3}, integration_times(intLoop)))
    
%     a = squeeze(rval_shift(:,:,intLoop));
%     subplot(ni, nvar, ind); ind = ind+1;
%     plot(a')
%     title(sprintf('rval shift  %.2f', integration_times(intLoop)))
%     
%     a = squeeze(rsig_shift(:,:,intLoop));
%     subplot(ni, nvar, ind); ind = ind+1;
%     plot(a')
%     title(sprintf('rsig shift  %.2f', integration_times(intLoop)))
    
end
end