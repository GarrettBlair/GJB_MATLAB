%%
ddir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\';
animals = {'Acc20832', 'Acc19947', 'Hipp18240'};
sessions = cell(3,1);
sessions{1} = {'HAB2' 'TR0' 'TR3' 'TR6' 'RET0' 'NEW0' 'NEW4'};
sessions{2} = {'HAB2' 'TR0' 'TR3' 'TR6' 'RET0' 'NEW0' 'NEW4'};
sessions{3} = {'TR0' 'TR3' 'TR6' 'WTR10' 'TR11' 'TR17'};

datadir = [ddir 'ipos_decoding\'];
figdir  = [ddir 'ipos_decoding_figs\'];

ipos_int_time = .25;
x_fold_training = 5;
for a = 1:length(animals)
    for s = 1:length(sessions{a})
        sessname = sprintf('%s_%s', animals{a}, sessions{a}{s});
%         disp(sessname)
        ipos_file = sprintf('%s%s\\processed_files\\ipos_files\\%s_int_%1.2fsec.mat',...
            ddir, animals{a}, sessname, ipos_int_time);
        dataSaveName = sprintf('%s%s_decoding_%1.2fsec.mat', datadir, sessname, ipos_int_time);
        figSaveName = sprintf('%sa%d_s%d_%sdecoding_%1.2fsec.png', figdir, a, s, sessname, ipos_int_time);
        if ~isfile(dataSaveName)
            fprintf('\n\t\t~~~~~ EVALUATING %s', sessname)
            ipos_svm_decoding(ipos_file, dataSaveName, figSaveName, x_fold_training)
        else
            fprintf('\n\t\t~~~~~ SKIP %s', sessname)
        end
    end
end

%%
plot_only = false;
if plot_only == true
for a = 1:length(animals)
    for s = 1:length(sessions{a})
        sessname = sprintf('%s_%s', animals{a}, sessions{a}{s});
        ipos_file = sprintf('%s%s\\processed_files\\ipos_files\\%s_int_%1.2fsec.mat',...
            ddir, animals{a}, sessname, ipos_int_time);
        dataSaveName = sprintf('%s%s_decoding_%1.2fsec.mat', datadir, sessname, ipos_int_time);
        load(dataSaveName)
        figSaveName = sprintf('%sa%d_s%d_%sdecoding_%1.2fsec.png', figdir, a, s, sessname, ipos_int_time);
        %%
        fig = figure;
        set(fig, 'Position', [216 351 931 454])
        clrstyle = [.1 .1 .1; .5 .2 .8; .8 .2 .2; .2 .2 .7];
        subplot(2,2,1); cla; hold on;
        room_angle_means = summary_figs(room_angle_spks, room_angle_ipos, room_angle_iposr, room_angle_iposa, clrstyle);
        title('Room Frame Angle'); set(gca, 'XTickLabel', {'spks' 'ipos(rm-ar)' 'ipos(room)' 'ipos(arena)'})
        % subplot(2,1,2); hold on;
        % dist_means = summary_figs(room_dist_spks, room_dist_ipos, room_dist_iposr, room_dist_iposa);
        % title('Room Frame Shock Dist.'); set(gca, 'XTickLabel', {'spks' 'ipos(room)' 'ipos(arena)' 'ipos(room-arena)'})
        subplot(2,2,3); cla; hold on;
        room_pos_means = summary_figs(room_pos_spks, room_pos_ipos, room_pos_iposr, room_pos_iposa, clrstyle);
        ylim([-5 85])
        title('Room Frame Position'); set(gca, 'XTickLabel', {'spks' 'ipos(rm-ar)' 'ipos(room)' 'ipos(arena)'}, 'XTickLabelRotation', 0)
        
        subplot(2,2,2); cla; hold on;
        arena_angle_means = summary_figs(arena_angle_spks, arena_angle_ipos, arena_angle_iposr, arena_angle_iposa, clrstyle);
        ylim([-.2 3.5])
        title('Arena Frame Angle'); set(gca, 'XTickLabel', {'spks' 'ipos(rm-ar)' 'ipos(room)' 'ipos(arena)'})
        % subplot(2,1,2); hold on;
        % dist_means = summary_figs(arena_dist_spks, arena_dist_ipos, arena_dist_iposr, arena_dist_iposa);
        % title('Room Frame Shock Dist.'); set(gca, 'XTickLabel', {'spks' 'ipos(room)' 'ipos(arena)' 'ipos(room-arena)'})
        subplot(2,2,4); cla; hold on;
        arena_pos_means = summary_figs(arena_pos_spks, arena_pos_ipos, arena_pos_iposr, arena_pos_iposa, clrstyle);
        ylim([-5 85])
        title('Arena Frame Position'); set(gca, 'XTickLabel', {'spks' 'ipos(rm-ar)' 'ipos(room)' 'ipos(arena)'})
        
        saveas(fig, figSaveName)
    end
end
end

%%
function means = summary_figs(room_angle_spks, room_angle_ipos, room_angle_iposr, room_angle_iposa, clrstyle)
angdiff         = room_angle_spks.pred_diff';
angdiff_r       = room_angle_spks.rand_diff';
angdiff_ip      = room_angle_ipos.pred_diff';
angdiff_r_ip    = room_angle_ipos.rand_diff';
angdiff_ipr     = room_angle_iposr.pred_diff';
angdiff_r_ipr   = room_angle_iposr.rand_diff';
angdiff_ipa     = room_angle_iposa.pred_diff';
angdiff_r_ipa   = room_angle_iposa.rand_diff';
means = [];
means.diff_spks             = mean(angdiff);
means.diff_ipos             = mean(angdiff_ip);
means.diff_iposroom         = mean(angdiff_ipr);
means.diff_iposarena        = mean(angdiff_ipa);
means.rand_diff_spks        = mean(angdiff_r);
means.rand_diff_ipos        = mean(angdiff_r_ip);
means.rand_diff_iposroom    = mean(angdiff_r_ipr);
means.rand_diff_iposarena   = mean(angdiff_r_ipa);

% figure; hold on;
% plot(t, pos_angle)
% plot(t, room_angle_iposr.bin_center(room_angle_iposr.pred_real))
% plot(t, room_angle_spks.bin_center(room_angle_ipos.pred_real))

% figure(8); clf;
% hold on
% v = violinplot([angdiff angdiff_ip angdiff_ipr angdiff_ipa], [1 2 3 4]);%, 'Bandwidth', pi/12) ;
rms = .35;
plot([1+rms 1-rms], [median(angdiff_r)      median(angdiff_r)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([2+rms 2-rms], [median(angdiff_r_ip)     median(angdiff_r_ip)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([3+rms 3-rms], [median(angdiff_r_ipr)       median(angdiff_r_ipr)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([4+rms 4-rms], [median(angdiff_r_ipa)    median(angdiff_r_ipa)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
boxplot([angdiff angdiff_ip angdiff_ipr angdiff_ipa], 'Notch', 'off', 'Widths', .15, 'plotstyle', 'traditional', 'Jitter', .35, 'Symbol', '', 'Color', clrstyle)
drawnow
end