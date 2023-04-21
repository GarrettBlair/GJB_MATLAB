%%
ddir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\';
% animals = {'Acc20832', 'Acc19947', 'Hipp18240'};
animals = {'Acc20832', 'Acc19947', 'Hipp18240'};
% animals = {'Acc19947', 'Hipp18240'};
sessions = cell(length(animals),1);
% animals = {'Hipp18240'};
% sessions{1} = {'HAB2' 'TR0' 'TR3' 'TR6' 'RET0' 'NEW0' 'NEW4'};
% sessions{2} = {'HAB2' 'TR0' 'TR3' 'TR6' 'RET0' 'NEW0' 'NEW4'};
% sessions{3} = {'TR0' 'TR3' 'TR6' 'WTR10' 'TR11' 'TR17' 'DRK15'};
sessions{1} = {'TR0' 'TR3' 'TR6'};
sessions{2} = {'TR0' 'TR3' 'TR6'};
sessions{3} = {'TR0' 'TR3' 'TR6'};
% sessions{1} = {'TR3'};
% sessions{2} = {'TR3'};
% sessions{3} = {'TR3'};

% datadir = [ddir 'ipos_decoding\'];
% figdir  = [ddir 'ipos_decoding_figs\'];
% datadir = [ddir 'spk_decoding\'];
% figdir  = [ddir 'spk_decoding_figs\'];
datadir = [ddir 'spk2_decoding\'];
figdir  = [ddir 'spk2_decoding_figs\'];

ipos_int_time = .25;
x_fold_training = 10;
%%
tic
for a = 1:length(animals)
    for s = 1:length(sessions{a})
        for i = 1:length(ipos_int_time)
        %%
        sessname = sprintf('%s_%s', animals{a}, sessions{a}{s});
        %         disp(sessname)
        ipos_file = sprintf('%s%s\\processed_files\\ipos_files\\%s_int_%1.2fsec.mat',...
            ddir, animals{a}, sessname, ipos_int_time(i));
        if ~isfile(ipos_file)
            % use ipos_plotting.m to create files
            fprintf('\n\t\t~~~~~ MISSING %s', ipos_file)
        else
            dataSaveName = sprintf('%s%s_decoding_%1.2fsec.mat', datadir, sessname, ipos_int_time(i));
            figSaveName = sprintf('%sa%d_s%d_%sdecoding_%1.2fsec.png', figdir, a, s, sessname, ipos_int_time(i));
            if true%~isfile(dataSaveName)
                fprintf('\n\t\t~~~~~ EVALUATING %s\n', sessname)
                ipos_svm_decoding(ipos_file, dataSaveName, figSaveName, x_fold_training)
            else
                fprintf('\n\t\t~~~~~ SKIP %s', sessname)
            end
        end
        end
    end
end
toc
%%
plot_only = true;
vars = {'room_pos_spks' 'arena_pos_spks' 'room_pos_ipos' 'arena_pos_ipos' 'room_pos_ipos_mean' 'arena_pos_ipos_mean' 'room_pos_ipos_ensemble' 'arena_pos_ipos_ensemble'};
room_ipos.all = NaN(length(animals), length(sessions{1})); room_ipos.same = room_ipos.all; room_ipos.diff = room_ipos.all;
room_ipos.allr = room_ipos.all; room_ipos.samer = room_ipos.all; room_ipos.diffr = room_ipos.all;
room_spks = room_ipos; 
if plot_only == true
for a = 1:length(animals)
    for s = 1:length(sessions{a})
        sessname = sprintf('%s_%s', animals{a}, sessions{a}{s});
        ipos_file = sprintf('%s%s\\processed_files\\ipos_files\\%s_int_%1.2fsec.mat',...
            ddir, animals{a}, sessname, ipos_int_time);
        dataSaveName = sprintf('%s%s_decoding_%1.2fsec.mat', datadir, sessname, ipos_int_time);
        load(dataSaveName, vars{:})
        figSaveName = sprintf('%sa%d_s%d_%sdecoding_%1.2fsec.png', figdir, a, s, sessname, ipos_int_time);
        %%
        fig = figure;
        clf(fig)
        set(fig, 'Position', [56 351 931 454])
        clrstyle = [.1 .1 .1; .5 .2 .8; .8 .2 .2; .2 .2 .7];
        subplot(2,3,1); cla; hold on;
        room_all_means = summary_figs(room_pos_spks, room_pos_ipos, room_pos_ipos_mean, room_pos_ipos_ensemble, clrstyle, 'all');
        title('Room Frame Pos, all'); set(gca, 'XTickLabel', {'spks' 'ipos' 'mean' 'ensem'},'XTickLabelRotation', 0)
        ylim([-10 90])
        subplot(2,3,2); cla; hold on;
        room_same_means = summary_figs(room_pos_spks, room_pos_ipos, room_pos_ipos_mean, room_pos_ipos_ensemble, clrstyle, 'same');
        title('Room Frame Pos, same'); set(gca, 'XTickLabel', {'spks' 'ipos' 'mean' 'ensem'},'XTickLabelRotation', 0)
        ylim([-10 90])
        subplot(2,3,3); cla; hold on;
        room_diff_means = summary_figs(room_pos_spks, room_pos_ipos, room_pos_ipos_mean, room_pos_ipos_ensemble, clrstyle, 'diff');
        title('Room Frame Pos, diff'); set(gca, 'XTickLabel', {'spks' 'ipos' 'mean' 'ensem'},'XTickLabelRotation', 0)
        ylim([-10 90])
        
        subplot(2,3,4); cla; hold on;
        arena_all_means = summary_figs(arena_pos_spks, arena_pos_ipos, arena_pos_ipos_mean, arena_pos_ipos_ensemble, clrstyle, 'all');
        title('Arena Frame Pos, all'); set(gca, 'XTickLabel', {'spks' 'ipos' 'mean' 'ensem'},'XTickLabelRotation', 0)
        ylim([-10 90])
        subplot(2,3,5); cla; hold on;
        arena_same_means = summary_figs(arena_pos_spks, arena_pos_ipos, arena_pos_ipos_mean, arena_pos_ipos_ensemble, clrstyle, 'same');
        title('Arena Frame Pos, same'); set(gca, 'XTickLabel', {'spks' 'ipos' 'mean' 'ensem'},'XTickLabelRotation', 0)
        ylim([-10 90])
        subplot(2,3,6); cla; hold on;
        arena_diff_means = summary_figs(arena_pos_spks, arena_pos_ipos, arena_pos_ipos_mean, arena_pos_ipos_ensemble, clrstyle, 'diff');
        title('Arena Frame Pos, diff'); set(gca, 'XTickLabel', {'spks' 'ipos' 'mean' 'ensem'},'XTickLabelRotation', 0)
        ylim([-10 90])
        room_ipos.all(a,s) = room_all_means.arg2;
        room_ipos.same(a,s) = room_same_means.arg2;
        room_ipos.diff(a,s) = room_diff_means.arg2;
        room_ipos.allr(a,s) = room_all_means.arg2_rand;
        room_ipos.samer(a,s) = room_same_means.arg2_rand;
        room_ipos.diffr(a,s) = room_diff_means.arg2_rand;
                
        room_spks.all(a,s) = room_all_means.arg1;
        room_spks.same(a,s) = room_same_means.arg1;
        room_spks.diff(a,s) = room_diff_means.arg1;
        room_spks.allr(a,s) = room_all_means.arg1_rand;
        room_spks.samer(a,s) = room_same_means.arg1_rand;
        room_spks.diffr(a,s) = room_diff_means.arg1_rand;
        %         saveas(fig, figSaveName)
    end
end
end
figure(167); clf; hold on
plot(room_ipos.all', 'k')
plot(room_ipos.same, 'g')
plot(room_ipos.diff', 'm')

figure(168); clf; hold on
plot((room_ipos.allr-room_ipos.all)', 'k')
plot((room_ipos.samer-room_ipos.same)', 'g')
plot((room_ipos.diffr-room_ipos.diff)', 'm')

figure(169); clf; hold on
plot((room_spks.allr-room_spks.all)', 'k')
% plot((room_spks.samer-room_spks.same)', 'g')
% plot((room_spks.diffr-room_spks.diff)', 'm')
plot((room_spks.same)', 'g')
plot((room_spks.diff)', 'm')
%%
function means = summary_figs(arg1, arg2, arg3, arg4, clrstyle, refFrame)
switch refFrame
    case 'all'
        angdiff         = arg1.pred_err';
        angdiff_r       = arg1.rand_err';
        angdiff_ip      = arg2.pred_err';
        angdiff_r_ip    = arg2.rand_err';
        angdiff_ipr     = arg3.pred_err';
        angdiff_r_ipr   = arg3.rand_err';
        angdiff_ipa     = arg4.pred_err';
        angdiff_r_ipa   = arg4.rand_err';
    case 'same'
        angdiff         = arg1.pred_err_same';
        angdiff_r       = arg1.rand_err_same';
        angdiff_ip      = arg2.pred_err_same';
        angdiff_r_ip    = arg2.rand_err_same';
        angdiff_ipr     = arg3.pred_err_same';
        angdiff_r_ipr   = arg3.rand_err_same';
        angdiff_ipa     = arg4.pred_err_same';
        angdiff_r_ipa   = arg4.rand_err_same';
    case 'diff'
        angdiff         = arg1.pred_err_diff';
        angdiff_r       = arg1.rand_err_diff';
        angdiff_ip      = arg2.pred_err_diff';
        angdiff_r_ip    = arg2.rand_err_diff';
        angdiff_ipr     = arg3.pred_err_diff';
        angdiff_r_ipr   = arg3.rand_err_diff';
        angdiff_ipa     = arg4.pred_err_diff';
        angdiff_r_ipa   = arg4.rand_err_diff';
end
means = [];
means.arg1      = nanmedian(angdiff);
means.arg2      = nanmedian(angdiff_ip);
means.arg3      = nanmedian(angdiff_ipr);
means.arg4      = nanmedian(angdiff_ipa);
means.arg1_rand = nanmedian(angdiff_r);
means.arg2_rand = nanmedian(angdiff_r_ip);
means.arg3_rand = nanmedian(angdiff_r_ipr);
means.arg4_rand = nanmedian(angdiff_r_ipa);

% figure; hold on;
% plot(t, pos_angle)
% plot(t, room_angle_iposr.bin_center(room_angle_iposr.pred_real))
% plot(t, room_angle_spks.bin_center(room_angle_ipos.pred_real))

% figure(8); clf;
% hold on
% v = violinplot([angdiff angdiff_ip angdiff_ipr angdiff_ipa], [1 2 3 4]);%, 'Bandwidth', pi/12) ;
cla
v = violinplot([angdiff angdiff_ip angdiff_ipr angdiff_ipa], [1 2 3 4]);%, 'Bandwidth', pi/12) ;
for i = 1:4
v(i).EdgeColor = clrstyle(i,:)*0;
v(i).ShowData = 0;
v(i).ShowMean = 0;
v(i).ShowNotches = 0;
v(i).ViolinColor = clrstyle(i,:);
v(i).ViolinAlpha = .3;
rectangle('Position', [i-.45 0 .425 max(max([angdiff angdiff_ip angdiff_ipr angdiff_ipa]))], 'FaceColor', 'w', 'EdgeColor', 'w')
end
rms = .35;
plot([1+rms 1-rms], [nanmedian(angdiff_r)      nanmedian(angdiff_r)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([2+rms 2-rms], [nanmedian(angdiff_r_ip)     nanmedian(angdiff_r_ip)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([3+rms 3-rms], [nanmedian(angdiff_r_ipr)       nanmedian(angdiff_r_ipr)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([4+rms 4-rms], [nanmedian(angdiff_r_ipa)    nanmedian(angdiff_r_ipa)], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
boxplot([angdiff angdiff_ip angdiff_ipr angdiff_ipa], 'Notch', 'off', 'Widths', .2,...
    'plotstyle', 'traditional', 'Jitter', .35, 'Symbol', '', 'Color', clrstyle, 'Positions', [1 2 3 4]-.2)
axis([0.25 4.75 -.2 3.45 ])


drawnow
end
