function [outstr, randstr] = APA_GrangerCausality_HPCACC_epochs(hpc_fname, acc_fname)

twin = 10; % seconds, pre/post entrance or avoidance
peth_samples = 2*twin/.25 + 1;
% load("D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_17_H13_25_13_TR14_@placecells_ACC_miniscope2.mat")
load(acc_fname, 'ms', 'params')
ipos_timeres = params.ipos_int_time;% .1;
params.skip_ensemble = true;
[ms.room.momentary_pos_info,  acc_rtemp]    = Fenton_ipos(ms, ipos_timeres, 'room', params);
[ms.arena.momentary_pos_info, ~]            = Fenton_ipos(ms, ipos_timeres, 'arena', params);
a_svm = ms.room.svm_decoding;

ipos_acc = abs(ms.room.momentary_pos_info) - abs(ms.arena.momentary_pos_info);
ipos_a = nanmean(ipos_acc,1);

load(hpc_fname, 'ms')
params.skip_ensemble = true;
h_svm = ms.room.svm_decoding;
[ms.room.momentary_pos_info,  hpc_rtemp]    = Fenton_ipos(ms, ipos_timeres, 'room', params);
[ms.arena.momentary_pos_info, ~]            = Fenton_ipos(ms, ipos_timeres, 'arena', params);

ipos_hpc = abs(ms.room.momentary_pos_info) - abs(ms.arena.momentary_pos_info);
ipos_h = nanmean(ipos_hpc,1);

if length(ipos_a)<length(ipos_h) % interp to the smaller vector
    ipos_h = interp1(hpc_rtemp.t, ipos_h, acc_rtemp.t, 'linear')';
    decode_h = interp1(h_svm.t, h_svm.pred_err, a_svm.t, 'linear')';
    decode_a = a_svm.pred_err';
    t = acc_rtemp.t;
    room_struct = acc_rtemp;
elseif length(ipos_a)>length(ipos_h)
    ipos_a = interp1(acc_rtemp.t, ipos_a, hpc_rtemp.t, 'linear')';
    decode_a = interp1(a_svm.t, a_svm.pred_err, h_svm.t, 'linear')';
    decode_h = h_svm.pred_err';
    t = hpc_rtemp.t;
    room_struct = hpc_rtemp;
else 
    t = hpc_rtemp.t;
    room_struct = hpc_rtemp;
    decode_a = a_svm.pred_err';
    decode_h = h_svm.pred_err';
end
ipos_a(isnan(ipos_a)) = interp1(find(~isnan(ipos_a)), ipos_a(~isnan(ipos_a)), find(isnan(ipos_a)), 'linear');
ipos_h(isnan(ipos_h)) = interp1(find(~isnan(ipos_h)), ipos_h(~isnan(ipos_h)), find(isnan(ipos_h)), 'linear');
% decode_a(isnan(decode_a)) = interp1(find(~isnan(decode_a)), decode_a(~isnan(decode_a)), find(isnan(decode_a)), 'linear');
% decode_h(isnan(decode_h)) = interp1(find(~isnan(decode_h)), decode_h(~isnan(decode_h)), find(isnan(decode_h)), 'linear');

fig1 = figure; 
subplot(2,3,1:2)
plot(t, ipos_a); hold on; plot(t, ipos_h)
%%
% r = ms.room;


% t = rtemp.t;%ms.timestamps;
% bt = average_spks_time(t', params.ipos_int_time, t./1000, false, 'mean')
e = ms.room.entranceTimes./1000;



[th,rho] = cart2pol(room_struct.x, room_struct.y);
% approach_min_time_sec = 3;
if contains(hpc_fname, '_TR') && contains(acc_fname, '_TR')
    ZONE_CENTER = pi/2;
elseif contains(hpc_fname, '_CON') && contains(acc_fname, '_CON')
    ZONE_CENTER = 3*pi/2;
end
ZONE_SIZE = pi/6;% distance from zone center to be counted appr_start in the zone
ZONE_NEARTHRESH = ZONE_SIZE + (pi-ZONE_SIZE)/2;% distance from zone center to be counted appr_start an approach
[ang_dist, ~] = angular_distance(th, ZONE_CENTER, 0);
in_zone = ang_dist<=ZONE_SIZE;
near    = ang_dist<=ZONE_NEARTHRESH & ~in_zone;
far     = ang_dist> ZONE_NEARTHRESH & ~in_zone;

appr = [];
entr = [];

appr.start = find(far(1:end-1)==1 & far(2:end)==0); % find the transitions between far and near
appr.end   = find(far(1:end-1)==0 & far(2:end)==1); % find the transitions between far and near
% appr.start = find(far(1:end-1)==1 & far(2:end)==0); % find the transitions between far and near
% appr.end = find(near(1:end-1)==1 & far(2:end)==1); % find the transitions between near and far
appr.end = appr.end(appr.end>appr.start(1)); % ensure the first end comes after a start 
appr.start = appr.start(appr.start<appr.end(end)); % ensure the last start has an end

if length(appr.end) ~= length(appr.end)
    error('check inds here')
end
appr.length = t(appr.end) - t(appr.start);
appr.interT = [appr.start(1)-t(1);...
    t(appr.start(2:end)) - t(appr.end(1:end-1))];


appr.mindist     = NaN(size(appr.start)); % minimum distance during approach
appr.mindistind  = NaN(size(appr.start)); % index of that dist
appr.mindisttime = NaN(size(appr.start)); % time of that ind
appr.petha       = NaN(size(appr.start)); % time of that ind
appr.pethb       = NaN(size(appr.start)); % time of that ind

entr = appr; % shocked enraces

appr.isentrance  = false(size(appr.start)); % minimum distance during entrance
entr.entranceTimes = e; % minimum distance during entrance
entr.entrInd       = NaN(size(appr.start)); % index of the shock for entrance

plotting = true;
% figure(89); clf; 
subplot(2,3,5)
hold on
axis([-40 40 -40 40])
ipos_a_sub = [];
ipos_h_sub = [];
ipos_a_entr = [];
ipos_h_entr = [];
ipos_a_peth = [];
ipos_h_peth = [];
ipos_a_entrpeth = [];
ipos_h_entrpeth = [];

decode_a_peth = [];
decode_h_peth = [];

skip_peth_count = 0;
for i = 1:length(appr.start)
    md = min(ang_dist(appr.start(i):appr.end(i)));
    td = ismember(t, t(appr.start(i):appr.end(i)));
    isentr = e>=t(appr.start(i)) & e<=t(appr.end(i));
    if ~any(isentr)
        ind = find(min(abs(t-md)) == abs(t-md), 1, 'first');
        appr.mindist(i) = md;
        appr.mindistind(i) = find(ang_dist==md & td==true);
        appr.mindisttime(i) = t(ang_dist==md & td==true);
        
        
        a = appr.mindisttime(i) - twin;
        b = appr.mindisttime(i) + twin;
        appr.petha(i) = find(min(abs(t-a)) == abs(t-a));
        appr.pethb(i) = find(min(abs(t-b)) == abs(t-b));
        inds = appr.petha(i):appr.pethb(i);
        if length(inds) == peth_samples
        ipos_h_sub = cat(2, ipos_h_sub, ipos_h(inds));
        ipos_a_sub = cat(2, ipos_a_sub, ipos_a(inds));
        ipos_h_peth = cat(1, ipos_h_peth, ipos_h(inds));
        ipos_a_peth = cat(1, ipos_a_peth, ipos_a(inds));
        
        
        decode_a_peth = cat(1, decode_a_peth, decode_a(inds));
        decode_h_peth = cat(1, decode_h_peth, decode_h(inds));
            if plotting ==true
                plot3(room_struct.x(inds), room_struct.y(inds), t(inds), 'k')
            end
        else
            skip_peth_count = skip_peth_count+1;
            if length(appr.start)*.25 < skip_peth_count
                warning('Many approaches have uneven lengths')
            end
        end
    elseif any(isentr)
        isentr = find(isentr==1,1);
        ind = find(min(abs(t-e(isentr))) == abs(t-e(isentr)), 1, 'first');
        entr.entrInd(i) = ind;
        appr.isentrance(i) = true;
        entr.mindist(i) = md;
        entr.mindistind(i) = find(ang_dist==md & td==true);
        entr.mindisttime(i) = t(ang_dist==md & td==true);
        
        a = e(isentr) - twin;
        b = e(isentr) + twin;
        entr.petha(i) = find(min(abs(t-a)) == abs(t-a));
        entr.pethb(i) = find(min(abs(t-b)) == abs(t-b));
        
        inds = entr.petha(i):entr.pethb(i);
        if length(inds) == peth_samples
            ipos_h_entr = cat(2, ipos_h_entr, ipos_h(inds));
            ipos_a_entr = cat(2, ipos_a_entr, ipos_a(inds));
            ipos_h_entrpeth = cat(1, ipos_h_entrpeth, ipos_h(inds));
            ipos_a_entrpeth = cat(1, ipos_a_entrpeth, ipos_a(inds));
            if plotting == true
                plot3(room_struct.x(inds), room_struct.y(inds), t(inds), 'r.-')
            end
        else
            skip_peth_count = skip_peth_count+1;
            if length(appr.start)*.25 < skip_peth_count
                warning('Many approaches have uneven lengths')
            end
        end
    end
end

valid = false(length(appr.mindistind),1);
valid(1) = true;
for i = 2:length(appr.mindistind)
    last_valid = find(valid==true, 1, 'last');
    if appr.mindistind(i) - appr.mindistind(last_valid) < peth_samples
        valid(i) = false;
    else
        valid(i) = true;
    end
end
ipos_a_peth     = ipos_a_peth(valid, :);
ipos_a_sub      = ipos_a_sub(valid, :);
decode_a_peth   = decode_a_peth(valid, :);
ipos_h_peth     = ipos_h_peth(valid, :);
ipos_h_sub      = ipos_h_sub(valid, :);
decode_h_peth   = decode_h_peth(valid, :);

% appr.valid = appr.length>=approach_min_time_sec & ~appr.isentrance;
appr.valid = ~appr.isentrance && valid;
% appr.start = appr.start(valid);
% appr.end = appr.end(valid);
% appr.length = appr.length(valid);
% appr.mindist = appr.mindist(valid);
% appr.mindistind  = appr.mindistind(valid); 
% appr.mindisttime = appr.mindisttime(valid); 


entr.mindist     = entr.mindist(~isnan(entr.mindisttime));
entr.mindistind  = entr.mindistind(~isnan(entr.mindisttime));
entr.petha       = entr.petha(~isnan(entr.mindisttime));
entr.pethb       = entr.pethb(~isnan(entr.mindisttime));
entr.mindisttime = entr.mindisttime(~isnan(entr.mindisttime));

% subplot(122);
% hold on; 
% shadedErrorBar([], nanmean(ipos_a_entrpeth,1), nanstd(ipos_a_entrpeth,[],1)./(size(ipos_a_entrpeth,1)-1),'lineProps', 'm.-'); 
% shadedErrorBar([], nanmean(ipos_h_entrpeth,1), nanstd(ipos_h_entrpeth,[],1)./(size(ipos_h_entrpeth,1)-1),'lineProps', 'c.-'); 
% plot(ipos_a_peth')
%%
% X = [ipos_a_sub; ipos_h_sub];
fs = .25;
source_labels = {'ACC' 'HPC'};
% nsamples = size(ipos_h_peth,2);
nsamples = [];
nrand = 1000;
verbose = false;
plotting = true;
nt = .5 + peth_samples/2;

% X = NaN(2, size(ipos_a_peth,2), size(ipos_a_peth,1)); % source x samples x trials
X = NaN(2, nt-1, size(ipos_a_peth,1)); % source x samples x trials
X(1,:,:) = ipos_a_peth(:, 1:nt-1)';
X(2,:,:) = ipos_h_peth(:, 1:nt-1)';
[outstr_pre, randstr_pre] = granger_causality_Seth2015(X, fs, source_labels, nsamples, nrand, verbose, plotting);

% X = NaN(2, size(ipos_a_peth,2), size(ipos_a_peth,1)); % source x samples x trials
X = NaN(2, nt-1, size(ipos_a_peth,1)); % source x samples x trials
X(1,:,:) = ipos_a_peth(:, nt+1:end)';
X(2,:,:) = ipos_h_peth(:, nt+1:end)';
[outstr_post, randstr_post] = granger_causality_Seth2015(X, fs, source_labels, nsamples, nrand, verbose, plotting);

% X = NaN(2, size(ipos_h_entrpeth,2), size(ipos_h_entrpeth,1)); % source x samples x trials
% X(1,:,:) = ipos_a_entrpeth';
% X(2,:,:) = ipos_h_entrpeth';
% [entr_str, entr_rand] = granger_causality_Seth2015(X, fs, source_labels, nsamples, nrand, verbose, plotting);
%%
% subplot(121);
subplot(2,3,4); cla
hold on; 

plot([nt nt], [ 2*min([nanmean(ipos_a_peth,1), nanmean(ipos_h_peth,1)])  2*max([nanmean(ipos_a_peth,1), nanmean(ipos_h_peth,1)])], ...
    'k:', 'LineWidth', 2)
np = sqrt(size(ipos_a_peth,1));
shadedErrorBar([], nanmean(ipos_a_peth,1), nanstd(ipos_a_peth,[],1)./np,'lineProps', 'm-');
shadedErrorBar([], nanmean(ipos_h_peth,1), nanstd(ipos_h_peth,[],1)./np,'lineProps', 'c-');
axis tight
%
subplot(2,3,3); cla
i=1; j=2;
d = randstr_pre.F(:,i,j);
[c, bins] = histcounts(d, 100, 'Normalization', 'probability');
hold on;
bins = bins(1:end-1) + abs(mean(diff(bins)));
plot(bins, c, 'k-')
plot([outstr_pre.F(i,j) outstr_pre.F(i,j)], [0 1.2*max(c)], 'g-')
xlim([  -.001 1.2*max(bins) ])
ylim([  -.01 1.2*max(c) ])
title(sprintf('ACC->HPC GC pre=%2.3f\nGC post== %2.3f', outstr_pre.pval(i,j), outstr_post.pval(i,j)))

d = randstr_post.F(:,i,j);
[c, bins] = histcounts(d, 100, 'Normalization', 'probability');
bins = bins(1:end-1) + abs(mean(diff(bins)));
plot(bins, c, 'r-')
plot([outstr_post.F(i,j) outstr_post.F(i,j)], [0 1.2*max(c)], 'm-')



subplot(2,3,6); cla
i=2; j=1;

d = randstr_pre.F(:,i,j);
[c, bins] = histcounts(d, 100, 'Normalization', 'probability');
hold on;
bins = bins(1:end-1) + abs(mean(diff(bins)));
plot(bins, c, 'k-')
plot([outstr_pre.F(i,j) outstr_pre.F(i,j)], [0 1.2*max(c)], 'g-')
xlim([  -.001 1.2*max(bins) ])
ylim([  -.01 1.2*max(c) ])
title(sprintf('HPC->ACC GC pre=%2.3f\nGC post== %2.3f', outstr_pre.pval(i,j), outstr_post.pval(i,j)))

d = randstr_post.F(:,i,j);
[c, bins] = histcounts(d, 100, 'Normalization', 'probability');
bins = bins(1:end-1) + abs(mean(diff(bins)));
plot(bins, c, 'r-')
plot([outstr_post.F(i,j) outstr_post.F(i,j)], [0 1.2*max(c)], 'm-')


fig2 = figure; clf
subplot(221); hold on
imagesc(decode_a_peth)
plot([nt nt], [0 size(decode_h_peth,1)],'r:', 'LineWidth', 2)
axis tight

subplot(222); hold on
imagesc(decode_h_peth)
plot([nt nt], [0 size(decode_h_peth,1)], 'r:', 'LineWidth', 2)
axis tight

subplot(2,2,3:4); hold on
plot([nt nt], [10 40], 'k:', 'LineWidth', 2)
shadedErrorBar([], nanmean(decode_a_peth,1), nanstd(decode_a_peth,[],1)./np,'lineProps', 'm-');
shadedErrorBar([], nanmean(decode_h_peth,1), nanstd(decode_h_peth,[],1)./np,'lineProps', 'c-');
axis tight
% plot(nanmean(decode_a_peth,1));
% hold on
% plot(nanmean(decode_h_peth,1))

outstr = outstr_pre;
randstr = randstr_pre;

drawnow()



 