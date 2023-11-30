% accfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_28_H17_31_09_CON20_@placecells_ACC_miniscope2.mat';
% hpcfile = 'D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_28_H17_31_09_CON20_@placecells_HPC_miniscope1.mat';
% comboim = "F:\GarrettBlair\APA\hpcacc24500_8x_con20.tif";
accfile = "D:\GarrettBlair\APA\Acc20832\processed_files\2023_01_26_H16_10_07_NEW5_@placecells.mat";
hpcfile = "D:\GarrettBlair\APA\HPCACC24502\processed_files\2023_06_29_H14_56_34_CON21_@placecells_HPC_miniscope1.mat";
comboim = "F:\GarrettBlair\APA\hpcacc24502_acc20832_4x.tif";
hpc_ms = load(hpcfile);
acc_ms = load(accfile);

r = hpc_ms.ms.room;
a = hpc_ms.ms.arena;
ds = 4;

% him = "F:\GarrettBlair\APA\hpc24500_con20_msCam_MC.tiff";
% aim = "F:\GarrettBlair\APA\acc24500_con20_msCam_MC.tiff";
kern = gausswin(15); kern(1:floor(length(kern)/2)) = 0; %kern(1:floor(length(kern)/2)).^5; 
% kern = sqrt(kern);
kern = kern./sum(kern);
% sp = hpc_ms.ms.room.svm_decoding.spks_bin(1:50,:)>0;
% hsp = hpc_ms.ms.neuron.C+hpc_ms.ms.neuron.YrA/2;
hsp = hpc_ms.ms.spks;
hsp = bin_spks_average(hsp, ds, false);
hsp = normalize_rows(hsp);
ht = hpc_ms.ms.timestamps./1000;
% asp = acc_ms.ms.neuron.C+acc_ms.ms.neuron.YrA/2;
asp = acc_ms.ms.spks;
asp = bin_spks_average(asp, ds, false);
asp = normalize_rows(asp);
for i = 1:size(hsp,1); hsp(i,:) = conv(hsp(i,:), kern, 'same'); end
for i = 1:size(asp,1); asp(i,:) = conv(asp(i,:), kern, 'same'); end
at = acc_ms.ms.timestamps./1000;
t = bin_spks_average(at', ds, false);
t = round(t-t(1));
%%
win = 90;
nsegs = 30;
spkscale = 3;

figure(12345); clf
set(gcf, 'Color', 'k', 'Position', [200 182 1400 650])
colormap magma
c = [.2 .9 1];
m = [1 .2 .9];
for i = win+1:2:2000-win
    %%
    clf
    fr = double(imread(comboim, i));
    fr = fr(:, 20:end-50);
    hspks = hsp(1:nsegs, i-win:i+win);
    aspks = asp(1:nsegs, i-win:i+win);
    b = NaN(1,size(aspks,2));
    subplot_tight(1,11,1:8)
    imagesc(fr, [0 255]); axis image off
    text(40+69, 20, 'CA1', 'Color', c, 'FontSize', 32, 'FontWeight', 'bold')
    text(420, 20, 'ACC', 'Color', m, 'FontSize', 32, 'FontWeight', 'bold')
    text(210, 280, '8x speed', 'Color', 'w', 'FontSize', 18, 'FontWeight', 'normal')
    text(390, 280, sprintf('%d sec', t(i)), 'Color', 'w', 'FontSize', 18, 'FontWeight', 'normal')
    subplot_tight(1,11,9:11)
    stacked_traces([aspks*NaN; b; hspks], spkscale, {'Color', c, 'LineWidth', 1.5}, c./4); axis tight 
    stacked_traces([aspks; b; hspks*NaN], spkscale, {'Color', m, 'LineWidth', 1.5}, m./4); axis tight 
    hold on
    plot([win+1 win+1] , [-2 nsegs*2+size(b,1) + 4], 'w:', 'LineWidth', 2)
    set(gca, 'YTick', [nsegs/2, nsegs*3/2 + size(b,1)], 'YTickLabel', {'ACC cells', 'CA1 cells'},...
        'YColor', 'w', 'Color', 'none', 'XColor', 'k', 'YTickLabelRotation', 90, 'TickLength', [0 0], 'FontSize', 18)
    xlim([ -2, size(aspks,2)])
    ylim([ -6, nsegs*2+size(b,1) + 6])
    text(win+1-5, -4, sprintf('%d sec', t(i)), 'Color', 'w', 'FontSize', 18, 'FontWeight', 'normal')
    drawnow
    temp = getframe(gcf);
    imwrite(temp.cdata, 'F:\GarrettBlair\APA\test2.tiff', 'WriteMode', 'append')
end
%%
% hbg = imread(him, 1);
% abg = imread(aim, 1);
% fs =1:100:length(at);
% hbg = zeros(size(hbg));
% abg = zeros(size(abg));
% for i = fs
%     hbg = hbg + double(imread(him, i))./length(fs);
%     abg = abg + double(imread(aim, 1))./length(fs);
% 
% end
% abg = conv2(abg, ones(5,5)./25, 'same');
% hbg = conv2(hbg, ones(5,5)./25, 'same');
% %%
% timeclock = [1:1:60];
% hsmall = strel('disk', 3);
% hs = hsmall.Neighborhood./sum(hsmall.Neighborhood(:));
% hbig = strel('square', 10);
% hb = hbig.Neighborhood./sum(hbig.Neighborhood(:));
% for i = 1:length(timeclock)
%     %%
%     clf
%     hframe = find(min(abs(ht-timeclock(i))) == abs(ht-timeclock(i)), 1);
%     aframe = find(min(abs(at-timeclock(i))) == abs(at-timeclock(i)), 1);
%     hfr = double(imread(him, hframe));
%     afr = double(imread(aim, aframe));
%     hf = conv2(hfr, hs, 'same');
%     af = conv2(afr, hs, 'same');
%     hbg = conv2(hfr, hb, 'same');
%     abg = conv2(afr, hb, 'same');
%     
% %     hf = imclose(hfr, hsmall);
% %     hbg = imopen(hfr, hbig);
% %     af = imclose(afr, hsmall);
% %     abg = imopen(afr, hbig);
%     im = hfr-hbg; im = im(25:end-25, :); im = im(:, 25:end-25);
%     subplot_tight(1,2,1, [.05 .1]);
% %     imagesc(hf, [40 150])
%     imagesc(im, [-10 10])
%     subplot_tight(1,2,2, [.05 .1]);
% %     imagesc(af, [40 180])
%     im = afr-abg; im = im(25:end-25, :); im = im(:, 25:end-25);
%     imagesc(im, [-10 10])
%     drawnow
% end
% %%
% figure(1005); clf; 
% subplot_tight(2,3,3, [.05 .1])
% stacked_traces(hsp(1:15,:), 2.5, {'Color', cyan/1.5, 'LineWidth', 2}); axis tight
% set(gca, 'XTick', [t1 t2], 'XTickLabel', [60 120], 'YTick', [2 16], 'YTickLabel', [1 15])
% subplot_tight(2,3,6, [.05 .1])
% stacked_traces(asp(1:15,:), 2.5, {'Color', magenta/1.5, 'LineWidth', 2}); axis tight
% set(gca, 'XTick', [t1 t2], 'XTickLabel', [60 120], 'YTick', [2 16], 'YTickLabel', [1 15])
% 
% dims = double([length(hpc_ms.ms.frameNum), hpc_ms.ms.neuron.dims])
% %%
% r.pcell_stats.coherence = pfield_coherence_calc(r.pfields, r.vmap);
% a.pcell_stats.coherence = pfield_coherence_calc(a.pfields, a.vmap);
% % pull out cell examples
% figure(1001); clf; hold on
% subplot(1,2,1); hold on
% shksx = interp1(hpc_ms.ms.timestamps, r.x, r.entranceTimes, 'linear');
% shksy = interp1(hpc_ms.ms.timestamps, r.y, r.entranceTimes, 'linear');
% p = patch([0, -60*cos(pi/3), 60*cos(pi/3)], [0 60*sin(pi/3) 60*sin(pi/3)], 'r'); p.FaceAlpha = .3; p.EdgeColor='none';
% plot(r.x,r.y, 'Color', [.4 .4 .4], 'LineWidth', .5); axis([-44 44 -44 44])
% scatter(shksx,shksy, 50, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 .5 .5])
% axis square
% set(gca, 'Color', 'none')
% 
% subplot(1,2,2); hold on
% shksx = interp1(hpc_ms.ms.timestamps, a.x, r.entranceTimes, 'linear');
% shksy = interp1(hpc_ms.ms.timestamps, a.y, r.entranceTimes, 'linear');
% plot(a.x,a.y, 'Color', [.4 .4 .4], 'LineWidth', .5); axis([-44 44 -44 44])
% scatter(shksx,shksy, 50, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 .5 .5])
% axis square
% set(gca, 'Color', 'none')
% 
% ncells = 4;
% 
% goodc = (r.pcell_stats.coherence >= quantile(r.pcell_stats.coherence, .75));
% goodi = (r.pcell_stats.infoPerSpike >= quantile(r.pcell_stats.infoPerSpike, .75));
% goodcr = (r.split_corr >= quantile(r.split_corr, .75));
% good_r = find(goodc &  goodi & goodcr);
% 
% goodc = (a.pcell_stats.coherence >= quantile(a.pcell_stats.coherence, .75));
% goodi = (a.pcell_stats.infoPerSpike >= quantile(a.pcell_stats.infoPerSpike, .75));
% goodcr = (a.split_corr >= quantile(a.split_corr, .75));
% good_a = find(goodc &  goodi & goodcr);
% %
% r_pf = [];
% a_pf = [];
% % rord = randperm(length(good_r));
% rord = [11    12     6     8     9    10     5     3]; % for repro
% % aord = randperm(length(good_a));
% aord = [4     8     2     7     6     5     9    10]; % for repro
% for i = 1:ncells
% % r_pf = [squeeze(r.pfields_smooth_split1(good_r(i), :,:)); squeeze(r.pfields_smooth_split2(good_r(i), :,:))]; % squeeze(r.pfields_smooth(good, :,:));
% % a_pf = [squeeze(a.pfields_smooth_split1(good_a(i), :,:)); squeeze(a.pfields_smooth_split2(good_a(i), :,:))]; % squeeze(a.pfields_smooth(good, :,:));
% temp1 = squeeze(r.pfields_smooth(good_r(rord(i)), :,:));
% temp2 = squeeze(a.pfields_smooth(good_r(rord(i)), :,:));
% temp = cat(2, temp1, NaN(size(temp1,1),2), temp2);
% temp = normalize_matrix(temp);
% r_pf = [r_pf; NaN(4, size(temp,2)); temp];
% 
% temp1 = squeeze(r.pfields_smooth(good_a(aord(i)), :,:));
% temp2 = squeeze(a.pfields_smooth(good_a(aord(i)), :,:));
% temp = cat(2, temp1, NaN(size(temp1,1),2), temp2);
% temp = normalize_matrix(temp)*1;
% a_pf = [a_pf; NaN(4, size(temp,2)); temp];
% 
% end
% a_pf = cat(2, a_pf, NaN(size(a_pf,1),1));
% r_pf = cat(2, r_pf, NaN(size(r_pf,1),1));
% 
% % subplot(9,2,[3:3+(ncells-1)*2]);
% figure(1002); clf
% p = pcolor([r_pf a_pf]);
% axis image off
% p.EdgeColor = 'none';
% set(gca, 'Color', 'none')
% colormap viridis
% %%
% hpc_ms.params.skip_ensemble = true;
% h_time = hpc_ms.ms.room.svm_decoding.t;
% a_time = acc_ms.ms.room.svm_decoding.t;
% 
% % [ang_dist, inside_ang_ref] = angular_distance(ang1, ang_ref, bleed_size);
% 
% h_ipos = abs(hpc_ms.ms.room.momentary_pos_info) - abs(hpc_ms.ms.arena.momentary_pos_info);
% a_ipos = abs(acc_ms.ms.room.momentary_pos_info) - abs(acc_ms.ms.arena.momentary_pos_info);
% h_ipos_infer = hpc_ms.ms.room.inferred_ipos - hpc_ms.ms.arena.inferred_ipos;
% a_ipos_infer = acc_ms.ms.room.inferred_ipos - acc_ms.ms.arena.inferred_ipos;
% h_spks = hpc_ms.ms.room.svm_decoding.spks_bin; 
% a_spks = acc_ms.ms.room.svm_decoding.spks_bin;
% goodh = ~any(isnan([h_ipos]),1);
% gooda = ~any(isnan([a_ipos]),1);
% goodhi = ~any(isnan([h_ipos_infer]),1);
% goodai = ~any(isnan([a_ipos_infer]),1);
% if length(h_time)<length(a_time)
%     t = a_time(goodh);
%     h_ipos_match = NaN(size(h_ipos, 1), sum(gooda));
%     h_ipos_match_infer = NaN(size(h_ipos_infer, 1), sum(goodai));
%     h_spk_match = NaN(size(h_ipos, 1), sum(gooda));
%     for i = 1:size(h_ipos)
%         h_ipos_match(i,:) = interp1(h_time(goodh), h_ipos(i, goodh), t, 'linear')';
%         h_ipos_match_infer(i,:) = interp1(h_time(goodhi), h_ipos_infer(i, goodhi), a_time(goodai), 'linear')';
%         h_spk_match(i,:)  = interp1(h_time(goodh), h_spks(i, goodh), t, 'linear')';
%     end
%     a_ipos_match = a_ipos(:,gooda);
%     a_ipos_match_infer = a_ipos_infer(:,goodai);
%     a_spk_match = a_spks(:,gooda);
%     
%     h_ipos_mean = interp1(h_time(goodh), nanmean(h_ipos(:, goodh),1), t, 'linear')';
%     a_ipos_mean = nanmean(a_ipos(:, gooda),1);
% else
%     t = h_time(goodh);
%     a_ipos_match = NaN(size(a_ipos, 1), sum(goodh));
%     a_ipos_match_infer = NaN(size(a_ipos_infer, 1), sum(goodhi));
%     a_spk_match = NaN(size(a_ipos, 1), sum(goodh));
%     for i = 1:size(a_ipos)
%         a_ipos_match(i,:) = interp1(a_time(gooda), a_ipos(i, gooda), t, 'linear')';
%         a_ipos_match_infer(i,:) = interp1(a_time(goodai), a_ipos_infer(i, goodai), h_time(goodhi), 'linear')';
%         a_spk_match(i,:)  = interp1(a_time(gooda), a_spks(i, gooda), t, 'linear')';
%     end
%     h_ipos_match = h_ipos(:,goodh);
%     h_ipos_match_infer = h_ipos_infer(:,goodhi);
%     h_spk_match = h_spks(:,goodh);
%     
%     a_ipos_mean = interp1(a_time(gooda), nanmean(a_ipos(:, gooda),1), h_time(goodh), 'linear')';
%     h_ipos_mean = nanmean(h_ipos(:, goodh),1);
% end
% % [h_ipos] = Fenton_ipos(hpc_ms.ms, .25, 'arena', hpc_ms.params);
% 
% goodinds = ~any(isnan([h_ipos_match; a_ipos_match]),1);
% goodindsi = ~any(isnan([h_ipos_match_infer; a_ipos_match_infer]),1);
% % h_ipos = zscore( h_ipos_match(:, goodinds) );
% % a_ipos = zscore( a_ipos_match(:, goodinds) );
% h_ipos = ( h_ipos_match(:, goodinds) );
% a_ipos = ( a_ipos_match(:, goodinds) );
% h_ipos_infer = ( h_ipos_match_infer(:, goodindsi) );
% a_ipos_infer = ( a_ipos_match_infer(:, goodindsi) );
% [hpc_acc_iposcorr, hpc_acc_iposcorrpval] = nancorr( nanmean(h_ipos,1), nanmean(a_ipos,1) );
% [hpc_acc_iposcorr_infer, hpc_acc_iposcorrpval_infer] = nancorr( nanmean(h_ipos_infer,1), nanmean(a_ipos_infer,1) );
% h_ipos = (bin_spks_average(h_ipos, 5, false));
% a_ipos = (bin_spks_average(a_ipos, 5, false));
% h_ipos_infer = (bin_spks_average(h_ipos_infer, 5, false));
% a_ipos_infer = (bin_spks_average(a_ipos_infer, 5, false));
% 
% h_spks = ( h_spk_match(:, goodinds) );
% a_spks = ( a_spk_match(:, goodinds) );
% 
% hsp = normalize_rows(bin_spks_average(h_spks, 5, false));
% asp = normalize_rows(bin_spks_average(a_spks, 5, false));