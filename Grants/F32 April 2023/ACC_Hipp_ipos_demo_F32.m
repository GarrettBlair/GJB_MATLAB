%%
clear
animal = {'Hipp18240', 'Acc19947', 'Acc20832', 'Hipp18240', 'Acc19947', 'Acc20832'};
fname = {'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Hipp18240_TR3_decoding_0.25sec.mat',...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Acc19947_TR3_decoding_0.25sec.mat',...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Acc20832_TR3_decoding_0.25sec.mat',...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Hipp18240_TR6_decoding_0.25sec.mat',...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Acc19947_TR6_decoding_0.25sec.mat',...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Acc20832_TR6_decoding_0.25sec.mat'};
names = {'Rat 1 Hipp', 'Rat 2 Acc', 'Rat 2 Acc', 'Rat 1 Hipp', 'Rat 2 Acc', 'Rat 2 Acc'};
figure(1); clf;
figure(2); clf

nsess = 6;
%% CALCULATE THE CORELATION OF SIGNALS
if false;
mean_spkcorr = NaN(nsess,1);
mean_iposcorr = NaN(nsess,1);
figure(1000); clf
for jj = 6:nsess
    %%
    load(fname{jj}, 'room_angle_spks', 'ipos_file')
    ip = load(ipos_file);
    temp = load(['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\' animal{jj} '\processed_files\' ip.fname]);
    
    ipos = ip.ipos_ra;
    spks = ip.ms_room_temp.spks_bin;
    % % ipos_sm = NaN(size(ipos));
    % % for i = 1:size(ipos,1)
    % %     ipos_sm(i,:) = conv(ipos(i,:), ones(4,1)./4, 'same');
    % % end
    % ipos_mean = nanmean(ipos, 1);
    % cellcorr = corr([ipos_mean; spks]', 'type', 'Pearson');
    % cellcorr = cellcorr(2,:,1);
    % iposcorr = corr([ipos_mean; ipos]', 'type', 'Pearson');
    % iposcorr = iposcorr(2,:,1);
    
    [cellcorr,  cellcorr_prob] = corr(spks', 'type', 'Kendall');
    cellcorr(eye(size(cellcorr,1))==1) = NaN;
    cellcorr_prob(eye(size(cellcorr,1))==1) = NaN;
    [iposcorr,  iposcorr_prob] = corr(ipos', 'type', 'Kendall');
    iposcorr(eye(size(iposcorr,1))==1) = NaN;
    iposcorr_prob(eye(size(iposcorr,1))==1) = NaN;
    
    mean_spkcorr(jj) = nanmean(cellcorr(:));
    mean_iposcorr(jj) = nanmean(iposcorr(:));
    subplot(nsess, 1, jj)
    imagesc(ipos, [-.5 .5])
    % subplot(nsess, 2, 2*(jj-1) + 2)
    % imagesc(ipos, [-1 1])
end
% save('D:\F32 FIGS\cellcorr_ipos_corr_kendall.mat', 'mean_spkcorr', 'mean_iposcorr', 'animal', 'fname', 'names')
end
%%
figure(1); clf;
figure(2); clf
nsess = 6;
dist_corr = NaN(nsess,1);
numShks = NaN(nsess,1);
decodeAccuracy = NaN(nsess,1);
decodeRand = NaN(nsess,1);
decodeRandStd = NaN(nsess,1);
arena_perc = NaN(nsess,1);
room_perc = NaN(nsess,1);
segs_corr = NaN(nsess,1);

for jj = 1:nsess
    %%
    load(fname{jj}, 'room_angle_spks', 'ipos_file')
    ip = load(ipos_file);
    temp = load(['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\' animal{jj} '\processed_files\' ip.fname]);
    % temp = load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\2023_01_09_H14_32_47_TR3_@placecells.mat');
    %
    % aaa = room_pos_spks;
    aaa = room_angle_spks;
    % aaa = room_angle_;
    % % % [angcounts, ~, ~,  angbin, radbin] = histcounts2(aaa.theta, aaa.r, aaa.abins, aaa.rbins);
    % % % valid_bin = ~isnan(aaa.pred_same.*aaa.pred_diff);
    % % %
    % % % [angbin_same, radbin_same] = ind2sub(size(angcounts), aaa.pred_same);
    % % % [angbin_diff, radbin_diff] = ind2sub(size(angcounts), aaa.pred_diff);
    % % % [x, y] = pol2cart( aaa.theta, aaa.r);
    % % %
    % % % same_x = x*NaN; same_y = x*NaN;
    % % % diff_x = x*NaN; diff_y = x*NaN;
    % % %
    % % % [same_x(valid_bin), same_y(valid_bin)] = ...
    % % %     pol2cart( aaa.bin_center{1}(angbin_same(valid_bin)), aaa.bin_center{2}(radbin_same(valid_bin)));
    % % % [diff_x(valid_bin), diff_y(valid_bin)] = ...
    % % %     pol2cart( aaa.bin_center{1}(angbin_diff(valid_bin)), aaa.bin_center{2}(radbin_diff(valid_bin)));
    
    [theta, r] = cart2pol( ip.ms_room_temp.x, ip.ms_room_temp.y );
    
    ipos = ip.ipos_ra;
    spks = ip.ms_room_temp.spks_bin;
    ipos_sm = NaN(size(ipos));
    for i = 1:size(ipos,1)
        ipos_sm(i,:) = conv(ipos(i,:), ones(4,1)./4, 'same');
    end
    % ipos_mean = nanmean(ipos_sm, 1);
    ipos_mean = nanmean(ipos, 1);
    [~, ipos_runstest(jj)] = runstest(ipos_mean);
    
    % ipos_mean = nanmean(ipos, 1);
    % ipos_mean = nanmean(ip.ipos_smooth, 1);
    room_perc(jj) = sum(ipos_mean>0)./length(ipos_mean);
    arena_perc(jj) = sum(ipos_mean<0)./length(ipos_mean);
    f = ipos_mean; f(f>0) = 0;
    t_sec = ip.ms_room_temp.t;
    szd = aaa.theta - pi/2;
    szd = abs(mod(szd + pi, 2*pi) - pi);
    e = temp.ms.room.entranceTimes;
    shks = temp.ms.room.shockTimes;
    
    figure(1)
    subplot(nsess,1, jj); cla
    hold on;
    % area(1:length(ipos_mean), ipos_mean, 'FaceColor', 'r', 'EdgeColor', 'none')
    % area(1:length(ipos_mean), f, 'FaceColor', 'b', 'EdgeColor', 'none')
    % ttick = mod(t_sec, 400);
    % ttick = find(ttick(1:end-1) > ttick(2:end))+1;
    % set(gca, 'XTick', ttick, 'XtickLabel', round(t_sec(ttick)), 'YTick', [-.1 .0 .1])
    % axis tight
    % xlim([0 find(t_sec>=1200,1)])
    % ylim([-.15 .35])
    % yyaxis('right')
    % plot(szd, 'k')
    % ylim([-4*pi 2*pi])
    % title(sprintf('%s  -  Room: %0.3f%%     Arena: %0.3f%%      ncells: %d', names{jj}, room_perc(jj), arena_perc(jj), size(ipos,1)))
    % drawnow
    area(t_sec, ipos_mean, 'FaceColor', 'r', 'EdgeColor', 'none')
    area(t_sec, f, 'FaceColor', 'b', 'EdgeColor', 'none')
    ttick = mod(t_sec, 20);
    ttick = find(ttick(1:end-1) > ttick(2:end))+1;
    % ttick = [640:60:1000];
    xlim([640 1000])
    % set(gca, 'XTick', t_sec(ttick), 'XtickLabel', round(t_sec(ttick)-t_sec(ttick(1))))
    set(gca, 'YTick', [-.1 .0 .1])
    % axis tight
    % xlim([0 1200])
    ylim([-.3 .4])
    yyaxis('right')
    rectangle('Position', [0 0 t_sec(end) pi/6], 'FaceColor', [.8 .8 .8], 'EdgeColor', 'none')
    scatter(shks./1000,   shks*0,   'md', 'MarkerFaceColor', 'm', 'MarkerFaceAlpha', .4)
    % scatter(t_sec,   decoded_angle,  5, 'ro', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .4)
    % plot(t_sec,   abs(room_angle_spks.pred_err),  'r')
    plot(t_sec,   szd,  'k-')
    % plot(t_sec, szd, 'k')
    % plot([0 t_sec(end)], [0 0]+pi/6, 'k-')
    % plot([0 t_sec(end)], [0 0], 'k-')
    ylim([-4*pi 4])
    % set(gca, 'XTick', (ttick), 'XtickLabel', round(t_sec(ttick)-t_sec(ttick(1))), 'YTick', [0 pi/2 pi])
    set(gca, 'XTick', t_sec(ttick), 'XtickLabel', round(t_sec(ttick)-t_sec(find(t_sec>=640,1))))
    set(gca, 'YTick',  [0 pi/2 pi], 'YColor', 'k', 'YTickLabel', {'0' 'pi' 'pi/2'})
    title(sprintf('%s  -  Room: %0.3f%%     Arena: %0.3f%%      ncells: %d', names{jj}, room_perc(jj), arena_perc(jj), size(ipos,1)))
    % set(gca, )
    drawnow
    
    
    
    c = corr([ipos_mean; ipos]');
    segs_corr(jj) = mean(c(2:end,1));
    dist_corr(jj) = corr(szd', ipos_mean');
    numShks(jj) = length(temp.ms.room.shockTimes);
    decodeAccuracy(jj) = nanmean(abs(room_angle_spks.pred_err));
    decodeRand(jj) = nanmean(abs(room_angle_spks.rand_err));
    decodeRandStd(jj) = nanstd(abs(room_angle_spks.rand_err));
    
    if isfield(temp.ms.room, 'entranceTimes')
        e = temp.ms.room.entranceTimes;
        shks = temp.ms.room.shockTimes;
        ei = interp1(ip.ms_room_temp.t.*1000, 1:length(ip.ms_room_temp.t), e, 'nearest');
        % shock_zone_center = pi/2; % typical room shock configuration
        % shock_zone_size = pi/6; % size in rad from center to edge
        % distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure
        decoded_angle = room_angle_spks.bin_center(room_angle_spks.pred_real) - pi/2;
        decoded_dist = abs(mod(decoded_angle + pi, 2*pi) - pi);
        
        % sz_pos = (mod(szd + pi, 2*pi) - pi);
        approach = (szd(1:end-1)>pi/3 & szd(2:end)<=pi/3);
        ai = find(approach);
        exits = (szd(1:end-1)<=pi/3 & szd(2:end)>pi/3);
        exi = find(exits);
        if exi(1)<ai(1)
            exi = exi(2:end);
        end
        if ai(end)>exi(end)
            ai = ai(1:end-1);
        end
        %%% remove avoid indices that are too close in time
        good_ai = ai;
        good_exi = exi;
        % j=0;
        % while any(diff(t_sec(good_ai)) < 5) && j<100
        %     close_ind = find(diff(t_sec(good_ai)) < 5, 1) + 1;
        %     good_ai = setdiff(good_ai, good_ai(close_ind));
        %     j=j+1;
        %     if j ==100
        %         disp('looping out')
        %     end
        % end
        
        for i = 1:length(ai)-1
            if t_sec(ai(i+1))- t_sec(ai(i)) < 5
                good_ai = setdiff(ai, ai(i+1));
                good_exi = setdiff(exi, exi(i+1));
            end
        end
        ai = good_ai;
        exi = good_exi;
        
        is_avoid = true(length(ai),1);
        closest_pos = NaN(length(ai),1);
        for i = 1:length(ai)
            if any(ei >= ai(i) & (ei <= exi(i)))
                is_avoid(i) = false;
            end
            closest_pos(i) = find(szd(ai(i):exi(i)) == min(szd(ai(i):exi(i))), 1) + ai(i) - 1;
        end
        avoid_ind = closest_pos(is_avoid);
        shk_ind = closest_pos(~is_avoid);
        % scatter(avoid_ind, avoid_ind*0, 'go')
        % scatter(shk_ind,   shk_ind*0,   'mx')
        % % i = 20;
        % % avoid_peth = gb_PETH(ipos_mean, avoid_ind, i, 4);
        % % fail_peth = gb_PETH(ipos_mean, shk_ind, i, 4);
        % % figure(2);
        % % subplot(nsess,1, jj)
        % % hold on;
        % % plot(mean(avoid_peth'), 'g'); plot(mean(fail_peth'), 'r')
        % %
        % % ttick = mod(0:24, 4);
        % % ttick = find(ttick==0);
        % % set(gca, 'XTick', ttick, 'XtickLabel', [-20:4:4]*ip.integration_val, 'YTick', [-.025 .0 .025])
        % % plot([i+1 i+1], [-.05 .05], 'k')
        figure(2);
        subplot(nsess, 1, jj)
        [m, v] = make_occupancymap_2D(ip.ms_room_temp.x, ip.ms_room_temp.y, ipos_mean, ip.params.pos_bins, ip.params.pos_bins);
        % [m, v] = make_occupancymap_2D(ip.ms_room_temp.x, ip.ms_room_temp.y, ipos_mean, ip.params.pos_bins, ip.params.pos_bins);
        imagesc(m./v, [-.1 .1])
        axis image
    end
    
    % shk_ind = ei;
    %%%%%%%%%%%%% show a spatial rep of the ipos
    %%%%%%%%%%%%% show a corelation of ipos with shock zone dist
end
%%
figure(101); clf
set(gcf, 'Position', [225 518 1301 358])
hold on;
area(t_sec, ipos_mean, 'FaceColor', 'r', 'EdgeColor', 'none')
% area(t_sec, f, 'FaceColor', 'b', 'EdgeColor', 'none')
% plot(t_sec, ipos_mean, 'Color', 'r')
% plot(t_sec, f, 'Color', 'b')
ttick = mod(t_sec, 20);
ttick = find(ttick(1:end-1) > ttick(2:end))+1;
% ttick = [640:60:1000];
xlim([640 1000])
% set(gca, 'XTick', t_sec(ttick), 'XtickLabel', round(t_sec(ttick)-t_sec(ttick(1))))
set(gca, 'XTick', t_sec(ttick), 'XtickLabel', round(t_sec(ttick)-t_sec(find(t_sec>=640,1))))
set(gca, 'YTick', [-.1 .0 .1])
% axis tight
% xlim([0 1200])
ylim([-.3 .4])
set(gcf,'units','centimeters','position',[4,4,24,8]);
set(gca,'FontSize',10,'FontName','Arial');
% export_fig 'D:\F32 FIGS\ipos_shockdist_iposonly' -eps -dALLOWPSTRANSPARENCY

figure(103); clf
hold on
set(gcf, 'Position', [225         760        1301         116])
% yyaxis('right')
rectangle('Position', [0 0 t_sec(end) pi/6], 'FaceColor', [.8 .8 .8], 'EdgeColor', 'none')
scatter(shks./1000,   shks*0,   'md', 'MarkerFaceColor', 'm')
xlim([640 1000])
% scatter(t_sec,   decoded_angle,  5, 'ro', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .4)
% plot(t_sec,   abs(room_angle_spks.pred_err),  'r')
plot(t_sec,   szd,  'k-')
% plot(t_sec, szd, 'k')
% plot([0 t_sec(end)], [0 0]+pi/6, 'k-')
% plot([0 t_sec(end)], [0 0], 'k-')
ylim([-8 4])
% set(gca, 'XTick', (ttick), 'XtickLabel', round(t_sec(ttick)-t_sec(ttick(1))), 'YTick', [0 pi/2 pi])
set(gca, 'XTick', t_sec(ttick), 'XtickLabel', round(t_sec(ttick)-t_sec(find(t_sec>=640,1))))
set(gca, 'YTick',  [0 pi/2 pi], 'YColor', 'k', 'YTickLabel', {'0' 'pi' 'pi/2'})
title(sprintf('%s  -  Room: %0.3f%%     Arena: %0.3f%%      ncells: %d', names{jj}, room_perc(jj), arena_perc(jj), size(ipos,1)))
% set(gca, )
drawnow
% set(gcf,'units','centimeters','position',[4,4,24,8]);
set(gcf,'units','centimeters','position',[4,4,24,8]);
set(gca,'FontSize',10,'FontName','Arial');
% export_fig 'D:\F32 FIGS\ipos_shockdist_distonly' -eps
% export_fig 'D:\F32 FIGS\ipos_shockdist' -eps -transparent % -noSAFER % -dALLOWPSTRANSPARENCY
% export_fig 'D:\F32 FIGS\ipos_shockdist' -pdf% -transparent % -noSAFER % -dALLOWPSTRANSPARENCY

%%
figure;
pie([arena_perc(6), room_perc(6)])

figure(3); clf
hold on
xs = [-.1 0 .1];
markers = {'o' 'd' 's'};
% subplo
for jj = 1:3
    scatter(1+xs(jj), arena_perc(jj), 30, 'Marker', markers{jj}, 'MarkerFaceColor', [0 0 .8], 'MarkerEdgeColor', 'none')
    scatter(2+xs(jj), room_perc(jj), 30, 'Marker', markers{jj}, 'MarkerFaceColor', [.8 0 0], 'MarkerEdgeColor', 'none')
end
for jj = 4:6
    scatter(3+xs(jj-3), arena_perc(jj), 30, 'Marker', markers{jj-3}, 'MarkerFaceColor', [0 0 .8], 'MarkerEdgeColor', 'none')
    scatter(4+xs(jj-3), room_perc(jj), 30, 'Marker', markers{jj-3}, 'MarkerFaceColor', [.8 0 0], 'MarkerEdgeColor', 'none')
end
% boxplot([arena_perc' room_perc'])
% boxplot(room_perc)
axis([ .5 4.5 0 1])
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'Early', 'Early', 'Late', 'Late'}, 'YTick', [0:.4:.8])
ylabel('Proportion of time in the preferred frame')
axis square
set(gcf,'units','centimeters','position',[4,4,24,8]);
set(gca,'FontSize',10,'FontName','Arial');
export_fig 'D:\F32 FIGS\ipos_reference_occ' -eps -dALLOWPSTRANSPARENCY

%%

figure(56); clf
subplot(2,1,1)
hold on
rectangle('Position', [0 0 t_sec(end) pi/6], 'FaceColor', [.8 .8 .8], 'EdgeColor', 'none')
scatter(shks./1000,   shks*0,   'md', 'MarkerFaceColor', 'm', 'MarkerFaceAlpha', .4)
% scatter(t_sec,   decoded_angle+pi/2,  5, 'r.', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .4)
% plot(t_sec,   (room_angle_spks.theta),  'k')
scatter(t_sec,   decoded_dist,  15, 'r.', 'MarkerFaceColor', 'r')
% scatter(t_sec,   decoded_dist,  15, 'ro', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none')
% plot(t_sec,   decoded_dist, 'r-')
plot(t_sec,   szd,  'k')
xlim([0 110])
ylim([-.4 pi+.4])
set(gca,'XTick',0:30:120,'FontName','Arial');

subplot(2,1,2)
hold on
xs = [-.1 0 .1];
markers = {'o' 'd' 's'};

for jj = 1:3
%     scatter(1+xs(jj), decodeRand(jj) - decodeAccuracy(jj), 30, 'Marker', markers{jj}, 'MarkerFaceColor', [0 .8 .8], 'MarkerEdgeColor', 'none')
%     scatter(2+xs(jj), decodeRand(jj+3) - decodeAccuracy(jj+3), 30, 'Marker', markers{jj}, 'MarkerFaceColor', [.8 0 .8], 'MarkerEdgeColor', 'none')
    plot([1+xs(jj), 1+xs(jj)], [decodeRand(jj) decodeAccuracy(jj)], 'k-')
    plot([2+xs(jj), 2+xs(jj)], [decodeRand(jj+3) decodeAccuracy(jj+3)], 'k-')
    scatter(1+xs(jj), decodeAccuracy(jj), 10, 'Marker', markers{jj}, 'MarkerFaceColor', [0 .8 .8], 'MarkerEdgeColor', 'none')
    scatter(2+xs(jj), decodeAccuracy(jj+3), 10, 'Marker', markers{jj}, 'MarkerFaceColor', [.8 0 .8], 'MarkerEdgeColor', 'none')
    scatter(1+xs(jj), decodeRand(jj), 10, 'Marker', markers{jj}, 'MarkerFaceColor', [.2 .2 .2], 'MarkerEdgeColor', 'none')
    scatter(2+xs(jj), decodeRand(jj+3), 10, 'Marker', markers{jj}, 'MarkerFaceColor', [.2 .2 .2], 'MarkerEdgeColor', 'none')
end
% rand_sd = nanmean(decodeRandStd(:));
% plot([0 3], [0 0], 'k-')
% boxplot([arena_perc' room_perc'])
% boxplot(room_perc)
axis square
% axis([ .5 2.5 -pi/2 pi/2])
axis([ .5 2.5 0-.2 pi/2 + .2])
set(gca, 'XTick', [1 2], 'XTickLabel', {'Early', 'Late'}, 'YTick', [0:pi/4:pi/2])
% ylabel('Decoding performance over chance')
set(gcf,'units','centimeters','position',[4,4,10,4]);
set(gca,'FontSize',10,'FontName','Arial');
export_fig 'D:\F32 FIGS\decoding_perf' -eps -dALLOWPSTRANSPARENCY


%%
animal = {'Hipp18240', 'Acc19947', 'Acc20832'};
% fname = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Acc20832_TR6_decoding_0.25sec.mat';
% fname = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Acc19947_TR6_decoding_0.25sec.mat';
% fname = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Hipp18240_TR6_decoding_0.25sec.mat';
% fname = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Hipp18240_TR3_decoding_0.25sec.mat';
% fname = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Acc20832_TR3_decoding_0.25sec.mat';
fname = {'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Hipp18240_TR6_decoding_0.25sec.mat',...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Acc19947_TR6_decoding_0.25sec.mat',...
    'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\spk2_decoding\Acc20832_TR6_decoding_0.25sec.mat'};
names = {'Rat 1 Hipp', 'Rat 2 Acc', 'Rat 2 Acc'};
% load(fname, 'room_angle_spks', 'ipos_file')
% ip = load(ipos_file);
% % temp = load(['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\' ip.fname]);
% % temp = load(['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\' ip.fname]);
% temp = load(['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\' ip.fname]);
figure(4); clf;
figure(5); clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TR6

for jj = 1:3
    %%
    load(fname{jj}, 'room_angle_spks', 'ipos_file')
    ip = load(ipos_file);
    temp = load(['C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\' animal{jj} '\processed_files\' ip.fname]);
    % temp = load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc20832\processed_files\2023_01_09_H14_32_47_TR3_@placecells.mat');
    %
    % aaa = room_pos_spks;
    aaa = room_angle_spks;
    % aaa = room_angle_;
    % % % [angcounts, ~, ~,  angbin, radbin] = histcounts2(aaa.theta, aaa.r, aaa.abins, aaa.rbins);
    % % % valid_bin = ~isnan(aaa.pred_same.*aaa.pred_diff);
    % % %
    % % % [angbin_same, radbin_same] = ind2sub(size(angcounts), aaa.pred_same);
    % % % [angbin_diff, radbin_diff] = ind2sub(size(angcounts), aaa.pred_diff);
    % % % [x, y] = pol2cart( aaa.theta, aaa.r);
    % % %
    % % % same_x = x*NaN; same_y = x*NaN;
    % % % diff_x = x*NaN; diff_y = x*NaN;
    % % %
    % % % [same_x(valid_bin), same_y(valid_bin)] = ...
    % % %     pol2cart( aaa.bin_center{1}(angbin_same(valid_bin)), aaa.bin_center{2}(radbin_same(valid_bin)));
    % % % [diff_x(valid_bin), diff_y(valid_bin)] = ...
    % % %     pol2cart( aaa.bin_center{1}(angbin_diff(valid_bin)), aaa.bin_center{2}(radbin_diff(valid_bin)));
    
    [theta, r] = cart2pol( ip.ms_room_temp.x, ip.ms_room_temp.y );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TR6
    
    ipos = ip.ipos_ra;
    ipos_sm = NaN(size(ipos));
    for i = 1:size(ipos,1)
        ipos_sm(i,:) = conv(ipos(i,:), ones(4,1)./4, 'same');
    end
    ipos_mean = nanmean(ipos_sm, 1);
    % ipos_mean = nanmean(ipos, 1);
    % ipos_mean = nanmean(ip.ipos_smooth, 1);
    room_perc(jj) = sum(ipos_mean>0)./length(ipos_mean);
    arena_perc(jj) = sum(ipos_mean<0)./length(ipos_mean);
    f = ipos_mean; f(f>0) = 0;
    t_sec = ip.ms_room_temp.t;
    figure(4)
    subplot(nsess,1, jj)
    hold on;
    area(1:length(ipos_mean), ipos_mean, 'FaceColor', 'r', 'EdgeColor', 'none')
    area(1:length(ipos_mean), f, 'FaceColor', 'b', 'EdgeColor', 'none')
    ttick = mod(t_sec, 400);
    ttick = find(ttick(1:end-1) > ttick(2:end))+1;
    set(gca, 'XTick', ttick, 'XtickLabel', round(t_sec(ttick)), 'YTick', [-.1 .0 .1])
    axis tight
    xlim([0 find(t_sec>=1200,1)])
    ylim([-.15 .15])
    title(sprintf('%s  -  Room: %0.3f%%     Arena: %0.3f%%      ncells: %d', names{jj}, room_perc(jj), arena_perc(jj), size(ipos,1)))
    drawnow
    if isfield(temp.ms.room, 'entranceTimes')
        e = temp.ms.room.entranceTimes;
        ei = interp1(ip.ms_room_temp.t.*1000, 1:length(ip.ms_room_temp.t), e, 'nearest');
        % shock_zone_center = pi/2; % typical room shock configuration
        % shock_zone_size = pi/6; % size in rad from center to edge
        % distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure
        szd = aaa.theta - pi/2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TR6
        
        % sz_pos = (mod(szd + pi, 2*pi) - pi);
        szd = abs(mod(szd + pi, 2*pi) - pi);
        approach = (szd(1:end-1)>pi/3 & szd(2:end)<=pi/3);
        ai = find(approach);
        exits = (szd(1:end-1)<=pi/3 & szd(2:end)>pi/3);
        exi = find(exits);
        if exi(1)<ai(1)
            exi = exi(2:end);
        end
        if ai(end)>exi(end)
            ai = ai(1:end-1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TR6
        
        %%% remove avoid indices that are too close in time
        good_ai = ai;
        good_exi = exi;
        % j=0;
        % while any(diff(t_sec(good_ai)) < 5) && j<100
        %     close_ind = find(diff(t_sec(good_ai)) < 5, 1) + 1;
        %     good_ai = setdiff(good_ai, good_ai(close_ind));
        %     j=j+1;
        %     if j ==100
        %         disp('looping out')
        %     end
        % end
        
        for i = 1:length(ai)-1
            if t_sec(ai(i+1))- t_sec(ai(i)) < 5
                good_ai = setdiff(ai, ai(i+1));
                good_exi = setdiff(exi, exi(i+1));
            end
        end
        ai = good_ai;
        exi = good_exi;
        
        is_avoid = true(length(ai),1);
        closest_pos = NaN(length(ai),1);
        for i = 1:length(ai)
            if any(ei >= ai(i) & (ei <= exi(i)))
                is_avoid(i) = false;
            end
            closest_pos(i) = find(szd(ai(i):exi(i)) == min(szd(ai(i):exi(i))), 1) + ai(i) - 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TR6
        
        avoid_ind = closest_pos(is_avoid);
        shk_ind = closest_pos(~is_avoid);
        scatter(avoid_ind, avoid_ind*0, 'go')
        scatter(shk_ind,   shk_ind*0,   'mx')
        i = 20;
        avoid_peth = gb_PETH(ipos_mean, avoid_ind, i, 4);
        fail_peth = gb_PETH(ipos_mean, shk_ind, i, 4);
        figure(5);
        subplot(nsess,1, jj)
        hold on;
        plot(mean(avoid_peth'), 'g'); plot(mean(fail_peth'), 'r')
        plot([i+1 i+1], [-.05 .05], 'k')
        
        ttick = mod(0:24, 4);
        ttick = find(ttick==0);
        set(gca, 'XTick', ttick, 'XtickLabel', [-20:4:4]*ip.integration_val, 'YTick', [-.025 .0 .025])
        
    end
    % shk_ind = ei;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TR6
    
end
%%
figure(6); clf
hold on
xs = [-.1 0 .1];
markers = {'o' 'd' 's'};
% subplo
for jj = 1:3
    scatter(1+xs(jj), arena_perc(jj), 30, 'Marker', markers{jj}, 'MarkerFaceColor', [0 0 .8], 'MarkerEdgeColor', 'none')
    scatter(2+xs(jj), room_perc(jj), 30, 'Marker', markers{jj}, 'MarkerFaceColor', [.8 0 0], 'MarkerEdgeColor', 'none')
end
% boxplot([arena_perc' room_perc'])
% boxplot(room_perc)
axis([ .5 2.5 0 1])
set(gca, 'XTick', [1 2], 'XTickLabel', {'Arena', 'Room'}, 'YTick', [0:.4:.8])
ylabel('Proportion of time in the preferred frame')
axis square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TR6

%%
figure(1);
clf
hold on
i = 1:length(aaa.pred_real);
% i = 400:1450;
plot(i, aaa.pred_rand, 'b.')
plot(i, aaa.var_binned, 'k')

% % % figure(1);
% % % clf
% % % hold on
% % % i = 1:length(x);
% % % % i = 400:1450;
% % % plot3(i, x(i), y(i), 'k')
% % % plot3(i, same_x(i), same_y(i), 'b-o')
% % % plot3(i, diff_x(i), diff_y(i), 'r-o')

figure(2);
clf
% subplot(131)
hold on
decode_pos = aaa.bin_center(aaa.pred_real);
% plot(abs(aaa.pred_err), 'k')
plot(decode_pos, 'k.')
% yyaxis('right')
plot(theta, 'b')
plot(shk_ind, shk_ind*0, 'r*')
plot(avoid_ind, avoid_ind*0, 'g*')















