%%
clear
load('D:\Sample Data\ensemble_prob_linear\standard\256_0.mat');
% load("D:\Sample Data\ensemble_prob_linear\prop_samples\60_0.mat")
% load('D:\Sample Data\ensemble_prob_linear\prop_cells\20_0.mat')
state_smoothing = 3 % using movmedian, size n prior + current
%%
dt = abs(diff(t));
dt = [dt(1) dt];

ms = [];
params = [];
comm_frame = []; rare_frame = [];
% bins = [-10:20/20:10];
params.pfield_kernel_radius = 5;
params.occupancy_thresh = 1;
s = switch_binned;
idx = 1:length(s);
% construct_place_maps_2D(struct_in, x, y, dt, spks, bins, params)
inds = idx>sum(s==1);%s==0;
ms1 = construct_place_maps_2D(ms, x1(inds), y1(inds), dt(inds), bin_spks(:, inds), xbins, ybins, params);
ms2 = construct_place_maps_2D(ms, x2(inds), y2(inds), dt(inds), bin_spks(:, inds), xbins, ybins, params);
comm_frame.corr_s0 = ms1.split_corr;
rare_frame.corr_s0 = ms2.split_corr;
figure(104); clf; 
subplot(2,2,1);
hold on; 
histogram(ms1.split_corr, [-.5:.05:1], 'FaceColor', 'b'); 
histogram(ms2.split_corr, [-.5:.05:1], 'FaceColor', 'r');
title('Agnostic to switch, session pfield corr')
ylabel('session pfield corr')
inds = s==1;
ms1 = construct_place_maps_2D(ms, x1(inds), y1(inds), dt(inds), bin_spks(:, inds), xbins, ybins, params);
ms2 = construct_place_maps_2D(ms, x2(inds), y2(inds), dt(inds), bin_spks(:, inds), xbins, ybins, params);
comm_frame.corr_s1 = ms1.split_corr;
rare_frame.corr_s1 = ms2.split_corr;

subplot(2,2,2);
hold on; 
histogram(ms1.split_corr, [-.5:.05:1], 'Normalization', 'probability', 'FaceColor', 'b'); 
histogram(ms2.split_corr, [-.5:.05:1], 'Normalization', 'probability', 'FaceColor', 'r');
title('Red ref only')

% figure(1); clf; plot(t,x1); hold on; plot(t2, xx1)
% e = abs(ensem_av1) - abs(ensem_av2);
e = abs(ensem_av1+ensem_min1) - abs(ensem_av2+ensem_min2);
% e = abs(ensem_min1) - abs(ensem_min2);
if state_smoothing >0
e = movmedian(e, [state_smoothing 0]);
else
e = e;
end
% e = ensem_av1+ensem_min1 - (ensem_av2+ensem_min2);
% e = ensem_av1;
e = e-min(e);
binsize = .005;
e = e./max(e)+ binsize;
thresh_e = nanmedian(e)+2*std(e);
%         edges_of_bins  = [-.0015:.001:.0015;
edges_of_bins  = [binsize:binsize:1+binsize];
[m] = gb_conjoint_ensemble_model_fitting(e, thresh_e, edges_of_bins, true);
m.indep_sum = m.upper*m.upper' + m.lower*m.lower';
% m.product/(sum(denom(:)))

v1 = m.bin_center( (max(m.lower)==m.lower) );
v2 = m.bin_center( (max(m.upper)==m.upper) );
s2 = abs(e-v1) > abs(e-v2);


inds = s2==0;
ms1 = construct_place_maps_2D(ms, x1(inds), y1(inds), dt(inds), bin_spks(:, inds), xbins, ybins, params);
ms2 = construct_place_maps_2D(ms, x2(inds), y2(inds), dt(inds), bin_spks(:, inds), xbins, ybins, params);
comm_frame.corr_sr0 = ms1.split_corr;
rare_frame.corr_sr0 = ms2.split_corr;

figure(104); 
subplot(2,2,3);
hold on; 
histogram(ms1.split_corr, [-.5:.05:1], 'Normalization', 'probability', 'FaceColor', 'b'); 
histogram(ms2.split_corr, [-.5:.05:1], 'Normalization', 'probability', 'FaceColor', 'r');
title('Inferred switch using blue frame only')


inds = s2==1;
ms1 = construct_place_maps_2D(ms, x1(inds), y1(inds), dt(inds), bin_spks(:, inds), xbins, ybins, params);
ms2 = construct_place_maps_2D(ms, x2(inds), y2(inds), dt(inds), bin_spks(:, inds), xbins, ybins, params);
comm_frame.corr_sr1 = ms1.split_corr;
rare_frame.corr_sr1 = ms2.split_corr;

subplot(2,2,4);
hold on; 
histogram(ms1.split_corr, [-.5:.05:1], 'Normalization', 'probability', 'FaceColor', 'b'); 
histogram(ms2.split_corr, [-.5:.05:1], 'Normalization', 'probability', 'FaceColor', 'r');
title('Inferred switch using red frame only')

% comm_frame.delta_corr = nanmedian(comm_frame.corr_sr0 - comm_frame.corr_s0);
% rare_frame.delta_corr = nanmedian(rare_frame.corr_sr1 - rare_frame.corr_s0);
% comm_frame.delta_corr = nansum(comm_frame.corr_sr0 < comm_frame.corr_s0)/length(comm_frame.corr_s0);
% rare_frame.delta_corr = nansum(rare_frame.corr_sr1 < rare_frame.corr_s0)/length(rare_frame.corr_s1);
comm_frame.delta_corr = nansum(comm_frame.corr_sr0 > comm_frame.corr_s0)/nansum(comm_frame.corr_sr0 < comm_frame.corr_s0);
rare_frame.delta_corr = nansum(rare_frame.corr_sr1 > rare_frame.corr_s0)/nansum(rare_frame.corr_sr1 < rare_frame.corr_s0);
figure(105); clf; 
subplot(2,1,1); hold on
% plot(([comm_frame.corr_s0, comm_frame.corr_sr0, comm_frame.corr_sr1]'), 'k')
plot([1 2], [comm_frame.corr_s0, comm_frame.corr_sr0]', 'k')
violinplot([comm_frame.corr_s0, comm_frame.corr_sr0]);
title(sprintf('Common ref change: %2.3f', comm_frame.delta_corr))
subplot(2,1,2); hold on
plot([1 2], [rare_frame.corr_s0, rare_frame.corr_sr1]', 'r')
violinplot([rare_frame.corr_s0, rare_frame.corr_sr0]);
title(sprintf('Common ref change: %2.3f', rare_frame.delta_corr))

mean(s==s2)

figure(101); clf;
j=1;
maxval = 1.1*(max(abs(ensem_av1+ensem_min1) - abs(ensem_av2+ensem_min2)));
subplot(5,1,j); hold on; j=j+1;
plot(ensem_av1, 'r'); ylim([-.005 maxval])
subplot(5,1,j); hold on; j=j+1;
plot(ensem_min1, 'm'); ylim([-.005 maxval])
subplot(5,1,j); hold on; j=j+1;
plot(ensem_av2, 'b'); ylim([-.005 maxval])
subplot(5,1,j); hold on; j=j+1;
plot(ensem_min2, 'c'); ylim([-.005 maxval])
subplot(5,1,j); hold on; j=j+1;
ee = abs(ensem_av1+ensem_min1) - abs(ensem_av2+ensem_min2);
plot(s*maxval, '-' ,'Color', [.9 .1 .1], 'LineWidth', 1);
plot(ee, 'k');
ee = movmedian(ee, [state_smoothing 0]);
plot(ee, '-' ,'Color', [.9 .7 .4], 'LineWidth', 2); ylim([-.005 maxval])

%%

































