%%
% h = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\ipos_files\Hipp18240_TR3_int_0.25sec.mat");
% a = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\ipos_files\Acc19947_TR3_int_0.25sec.mat");
clear
h = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\ipos_files\Hipp18240_TR3_int_0.25sec.mat");
a = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\ipos_files\Acc19947_TR3_int_0.25sec.mat");

ad = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\ipos_decoding\Acc19947_TR3_decoding_0.25sec.mat");
hd = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\ipos_decoding\Hipp18240_TR3_decoding_0.25sec.mat");


%%
hinds = 3600:4300-1;
ainds = 400:400+length(hinds)-1;%3600:4200-1;
hsmooth = h.ipos_ra*NaN;
asmooth = a.ipos_ra*NaN;
for i = 1:size(h.ipos_ra,1)
hsmooth(i,:) = conv2(h.ipos_ra(i,:), ones(4,1)./4, 'same');
end
for i = 1:size(a.ipos_ra,1)
asmooth(i,:) = conv2(a.ipos_ra(i,:), ones(4,1)./4, 'same');
end
hsmooth = zscore(h.ipos_ra, 0, 2);
asmooth = zscore(a.ipos_ra, 0, 2);
hsmooth = h.ipos_ra;
asmooth = a.ipos_ra;
% hi = sum(h.ipos_ra,1);
% ai = sum(a.ipos_ra,1);
hi = nanmean(hsmooth,1);
ai = nanmean(asmooth,1);
hspks = normalize_rows(h.ms_room_temp.spks_bin);
aspks = normalize_rows(a.ms_room_temp.spks_bin);

figure(1); clf
subplot(6,1,1);
plot(nanmean(hspks(:,hinds), 1'), 'b'); axis off tight
subplot(6,1,2:3);
imagesc(hspks(:,hinds), [-1 1]); colormap viridis
% imagesc(hsmooth(:,hinds), [-10 10]); colormap viridis
axis tight image
set(gca, 'YTick', [1 size(h.ipos_ra,1)])
subplot(6,1,4);
plot(hi(hinds), 'b'); axis off tight
subplot(6,1,5:6);
imagesc(hsmooth(:,hinds), [-1 1]); colormap viridis
axis tight image
set(gca, 'YTick', [1 size(h.ipos_ra,1)])
% yyaxis('right')
% plot(hi, 'b'); % axis off
% ylim([ -1.5  .15]) 
% set(gca, 'YTick', [-max(hi) max(hi)])


figure(2); clf
subplot(6,1,1);
plot(nanmean(aspks(:,ainds), 1'), 'm'); axis off tight
subplot(6,1,2:3);
imagesc(aspks(:,ainds), [-1 1]); colormap viridis
% imagesc(hsmooth(:,hinds), [-10 10]); colormap viridis
axis tight image
set(gca, 'YTick', [1 size(a.ipos_ra,1)])
subplot(6,1,4);
plot(ai(ainds), 'm'); axis off tight
subplot(6,1,5:6);
imagesc(asmooth(:,ainds), [-1 1]); colormap plasma
axis tight image
set(gca, 'YTick', [1 size(a.ipos_ra,1)])

% % % figure(3); clf; 
% % % x = abs(diff([ 0 find(diff(hi>0)) length(hi) ])) ;
% % % posx = hi(1:end-1)>0 & hi(2:end)>0;
% % % histogram(x, [0:1:30], 'Normalization', 'probability')
% % % x = abs(diff([ 0 find(diff(ai>0)) length(ai) ])) ;
% % % posx = ai(1:end-1)>0 & ai(2:end)>0;
% % % hold on
% % % histogram(x, [0:1:30], 'Normalization', 'probability')

% figure(5); clf
% plot(ai(ainds), 'm'); % 
% hold on
% plot(hi(hinds), 'b'); % axis off
% plot([1 max(length(hinds), length(ainds))], [0 0], 'k'); % axis off

%%
figure(99); clf

temp = hd.room_angle_spks;
subplot(4,2,1); hold on
plot(hinds, temp.var_binned(hinds), 'k', 'LineWidth', 2)
scatter(hinds, temp.pred_real(hinds), 'o', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'r')
subplot(4,2,3); hold on
plot(hinds, temp.var_binned(hinds), 'k', 'LineWidth', 2)
temp = hd.room_angle_ipos;
scatter(hinds, temp.pred_real(hinds), 'o', 'MarkerFaceColor', 'm', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'm')

temp = ad.room_angle_spks;
subplot(4,2,2); hold on
plot(ainds, temp.var_binned(ainds), 'k', 'LineWidth', 2)
scatter(ainds, temp.pred_real(ainds), 20, 'o', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'b')
temp = ad.room_angle_ipos;
subplot(4,2,4); hold on
plot(ainds, temp.var_binned(ainds), 'k', 'LineWidth', 2)
scatter(ainds, temp.pred_real(ainds), 20, 'o', 'MarkerFaceColor', [0 .6 .8], 'MarkerFaceAlpha', .5, 'MarkerEdgeColor', [0 .4 .4])


subplot(2,2,3); hold on
title('0.25 sec integration')
hdecode_acc = [nanmedian(abs(hd.room_angle_spks.pred_err)), nanmedian(abs(hd.room_angle_ipos.pred_err)),...
              nanmedian(abs(hd.arena_angle_spks.pred_err)), nanmedian(abs(hd.arena_angle_ipos.pred_err))];
hdecode_rand = [nanmedian(abs(hd.room_angle_spks.rand_err)), nanmedian(abs(hd.room_angle_ipos.rand_err)),...
              nanmedian(abs(hd.arena_angle_spks.rand_err)), nanmedian(abs(hd.arena_angle_ipos.rand_err))];
plot(hdecode_acc', 'b-')
plot(hdecode_rand', 'b:', 'LineWidth', 3)
% axis([0 5 0 pi])
% set(gca, 'XTick', [1:4], 'XTickLabel', {'spks' 'ipos' 'spks' 'ipos'})

% subplot(2,2,4); hold on
adecode_acc = [nanmedian(abs(ad.room_angle_spks.pred_err)), nanmedian(abs(ad.room_angle_ipos.pred_err)),...
              nanmedian(abs(ad.arena_angle_spks.pred_err)), nanmedian(abs(ad.arena_angle_ipos.pred_err))];
adecode_rand = [nanmedian(abs(ad.room_angle_spks.rand_err)), nanmedian(abs(ad.room_angle_ipos.rand_err)),...
              nanmedian(abs(ad.arena_angle_spks.rand_err)), nanmedian(abs(ad.arena_angle_ipos.rand_err))];
plot(adecode_acc', 'm-')
plot(adecode_rand', 'm:', 'LineWidth', 3)
axis([0 5 0 pi])
set(gca, 'XTick', [1:4], 'XTickLabel', {'spks' 'ipos' 'spks' 'ipos'})
xlabel('Room frame                         Arena frame')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
h = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\ipos_files\Hipp18240_TR3_int_5.00sec.mat");
a = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\processed_files\ipos_files\Acc19947_TR3_int_5.00sec.mat");

ad = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\ipos_decoding\Acc19947_TR3_decoding_5.00sec.mat");
hd = load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\ipos_decoding\Hipp18240_TR3_decoding_5.00sec.mat");


%%
% hinds = 3600:4300-1;
% ainds = 400:400+length(hinds)-1;%3600:4200-1;
hinds = 3600:4300-1;
ainds = 400:400+length(hinds)-1;%3600:4200-1;
hinds = hinds*.25/5; hinds = unique(round(hinds));
ainds = ainds*.25/5; ainds = unique(round(ainds));
hsmooth = h.ipos_ra*NaN;
asmooth = a.ipos_ra*NaN;
for i = 1:size(h.ipos_ra,1)
hsmooth(i,:) = conv2(h.ipos_ra(i,:), ones(4,1)./4, 'same');
end
for i = 1:size(a.ipos_ra,1)
asmooth(i,:) = conv2(a.ipos_ra(i,:), ones(4,1)./4, 'same');
end
hsmooth = zscore(h.ipos_ra, 0, 2);
asmooth = zscore(a.ipos_ra, 0, 2);
hsmooth = h.ipos_ra;
asmooth = a.ipos_ra;
% hi = sum(h.ipos_ra,1);
% ai = sum(a.ipos_ra,1);
hi = nanmean(hsmooth,1);
ai = nanmean(asmooth,1);
hspks = normalize_rows(h.ms_room_temp.spks_bin);
aspks = normalize_rows(a.ms_room_temp.spks_bin);

% figure(1); clf
% subplot(6,1,1);
% plot(nanmean(hspks(:,hinds), 1'), 'b'); axis off tight
% subplot(6,1,2:3);
% imagesc(hspks(:,hinds), [-1 1]); colormap viridis
% % imagesc(hsmooth(:,hinds), [-10 10]); colormap viridis
% axis tight image
% set(gca, 'YTick', [1 size(h.ipos_ra,1)])
% subplot(6,1,4);
% plot(hi(hinds), 'b'); axis off tight
% subplot(6,1,5:6);
% imagesc(hsmooth(:,hinds), [-1 1]); colormap viridis
% axis tight image
% set(gca, 'YTick', [1 size(h.ipos_ra,1)])
% % yyaxis('right')
% % plot(hi, 'b'); % axis off
% % ylim([ -1.5  .15]) 
% % set(gca, 'YTick', [-max(hi) max(hi)])
% 
% 
% figure(2); clf
% subplot(6,1,1);
% plot(nanmean(aspks(:,ainds), 1'), 'm'); axis off tight
% subplot(6,1,2:3);
% imagesc(aspks(:,ainds), [-1 1]); colormap viridis
% % imagesc(hsmooth(:,hinds), [-10 10]); colormap viridis
% axis tight image
% set(gca, 'YTick', [1 size(a.ipos_ra,1)])
% subplot(6,1,4);
% plot(ai(ainds), 'm'); axis off tight
% subplot(6,1,5:6);
% imagesc(asmooth(:,ainds), [-1 1]); colormap plasma
% axis tight image
% set(gca, 'YTick', [1 size(a.ipos_ra,1)])

% % % figure(3); clf; 
% % % x = abs(diff([ 0 find(diff(hi>0)) length(hi) ])) ;
% % % posx = hi(1:end-1)>0 & hi(2:end)>0;
% % % histogram(x, [0:1:30], 'Normalization', 'probability')
% % % x = abs(diff([ 0 find(diff(ai>0)) length(ai) ])) ;
% % % posx = ai(1:end-1)>0 & ai(2:end)>0;
% % % hold on
% % % histogram(x, [0:1:30], 'Normalization', 'probability')

% figure(5); clf
% plot(ai(ainds), 'm'); % 
% hold on
% plot(hi(hinds), 'b'); % axis off
% plot([1 max(length(hinds), length(ainds))], [0 0], 'k'); % axis off

%%


% temp = hd.room_angle_spks;
% subplot(4,2,1); hold on
% plot(hinds, temp.var_binned(hinds), 'k', 'LineWidth', 2)
% scatter(hinds, temp.pred_real(hinds), 'o', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'r')
% subplot(4,2,3); hold on
% plot(hinds, temp.var_binned(hinds), 'k', 'LineWidth', 2)
% temp = hd.room_angle_ipos;
% scatter(hinds, temp.pred_real(hinds), 'o', 'MarkerFaceColor', 'm', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'm')
% 
% temp = ad.room_angle_spks;
% subplot(4,2,2); hold on
% plot(ainds, temp.var_binned(ainds), 'k', 'LineWidth', 2)
% scatter(ainds, temp.pred_real(ainds), 20, 'o', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'b')
% temp = ad.room_angle_ipos;
% subplot(4,2,4); hold on
% plot(ainds, temp.var_binned(ainds), 'k', 'LineWidth', 2)
% scatter(ainds, temp.pred_real(ainds), 20, 'o', 'MarkerFaceColor', [0 .6 .8], 'MarkerFaceAlpha', .5, 'MarkerEdgeColor', [0 .4 .4])


figure(99); 
subplot(2,2,4); hold on
title('5 sec integration')
hdecode_acc = [nanmedian(abs(hd.room_angle_spks.pred_err)), nanmedian(abs(hd.room_angle_ipos.pred_err)),...
              nanmedian(abs(hd.arena_angle_spks.pred_err)), nanmedian(abs(hd.arena_angle_ipos.pred_err))];
hdecode_rand = [nanmedian(abs(hd.room_angle_spks.rand_err)), nanmedian(abs(hd.room_angle_ipos.rand_err)),...
              nanmedian(abs(hd.arena_angle_spks.rand_err)), nanmedian(abs(hd.arena_angle_ipos.rand_err))];
plot(hdecode_acc', 'b-')
plot(hdecode_rand', 'b:', 'LineWidth', 3)
% axis([0 5 0 pi/2])
% set(gca, 'XTick', [1:4], 'XTickLabel', {'spks' 'ipos' 'spks' 'ipos'})

% subplot(2,2,4); hold on
adecode_acc = [nanmedian(abs(ad.room_angle_spks.pred_err)), nanmedian(abs(ad.room_angle_ipos.pred_err)),...
              nanmedian(abs(ad.arena_angle_spks.pred_err)), nanmedian(abs(ad.arena_angle_ipos.pred_err))];
adecode_rand = [nanmedian(abs(ad.room_angle_spks.rand_err)), nanmedian(abs(ad.room_angle_ipos.rand_err)),...
              nanmedian(abs(ad.arena_angle_spks.rand_err)), nanmedian(abs(ad.arena_angle_ipos.rand_err))];
plot(adecode_acc', 'm-')
plot(adecode_rand', 'm:', 'LineWidth', 3)
axis([0 5 0 pi])
set(gca, 'XTick', [1:4], 'XTickLabel', {'spks' 'ipos' 'spks' 'ipos'})
xlabel('Room frame                         Arena frame')







