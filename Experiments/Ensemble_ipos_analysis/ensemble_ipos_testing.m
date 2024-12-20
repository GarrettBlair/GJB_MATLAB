clear
load("D:\Sample Data\ensemble_prob\standard\ipos_ensemb128_0.mat");

%
[nsegs, nt] = size(bin_spks1);
nt = 48015;
% t = linspace(0, nt, nt)/30;

figure;

% xx1 = average_spks_time(x1(1:nt), .25, t, 0, 'mean');
% xx1 = xx1(1:end-1);
% yy1 = average_spks_time(y1(1:nt), .25, t, 0, 'mean');
% yy1 = yy1(1:end-1);
% xx2 = average_spks_time(x2(1:nt), .25, t, 0, 'mean');
% xx2 = xx2(1:end-1);
% yy2 = average_spks_time(y2(1:nt), .25, t, 0, 'mean');
% yy2 = yy2(1:end-1);

dt = abs(diff(t));
dt = [dt(1) dt];

ms = [];
params = [];
bins = [-10:20/20:10];
params.pfield_kernel_radius = 5;
params.occupancy_thresh = 0;

% construct_place_maps_2D(struct_in, x, y, dt, spks, bins, params)
ms1 = construct_place_maps_2D(ms, x1, y1, dt, bin_spks1, bins, params);
ms2 = construct_place_maps_2D(ms, x2, y2, dt, bin_spks2, bins, params);


% figure(1); clf; plot(t,x1); hold on; plot(t2, xx1)
figure; hold on; plot3(x1(1:300), y1(1:300), 1:300)
%%
ii = 1:length(x1);
ss = switch_binned==0;

idx = 1:length(s);
% construct_place_maps_2D(struct_in, x, y, dt, spks, bins, params)
inds = idx>sum(s==1);%s==0;
ms1 = construct_place_maps_2D(ms, x1(inds), y1(inds), dt(inds), bin_spks(:, inds), xbins, ybins, params);
ms2 = construct_place_maps_2D(ms, x2(inds), y2(inds), dt(inds), bin_spks(:, inds), xbins, ybins, params);
for seg = 1:100%nsegs
    %%
s = bin_spks(seg,:)>0;
figure(2);
clf; 
subplot(3,2,1)
imagesc([squeeze(pf1(:,:,seg)) squeeze(pf2(:,:,seg))])
title('input pfield')
subplot(3,2,3)
% imagesc([squeeze(ms2.pfields(seg,:,:)), squeeze(ms2.pfields_smooth(seg,:,:))])
imagesc([squeeze(ms1.pfields_smooth(seg,:,:)), squeeze(ms2.pfields_smooth(seg,:,:))])
title('pfield from activity')

subplot(3,1,3)
hold on
scatter3(ii(ss==1), x1(ss==1), y1(ss==1), 'r.', 'MarkerEdgeAlpha', .15)
scatter3(ii(ss==0), x2(ss==0), y2(ss==0), 'b.', 'MarkerEdgeAlpha', .15)
scatter3(ii(s&ss), x1(s&ss), y1(s&ss), 10, 'ro', 'MarkerEdgeAlpha', .4)
scatter3(ii(s&~ss), x2(s&~ss), y2(s&~ss), 10,'bo', 'MarkerEdgeAlpha', .4)
drawnow

subplot(3,2,[2,4])
hold on
scatter3(ii(ss==1), y1(ss==1), x1(ss==1), 'r.', 'MarkerEdgeAlpha', .15)
scatter3(ii(ss==0), y2(ss==0), x2(ss==0), 'b.', 'MarkerEdgeAlpha', .15)
scatter3(ii(s&ss), y1(s&ss), x1(s&ss), 10, 'ro', 'MarkerEdgeAlpha', .4)
scatter3(ii(s&~ss), y2(s&~ss), x2(s&~ss), 10,'bo', 'MarkerEdgeAlpha', .4)
drawnow
pause(.5)
% input('')
end


