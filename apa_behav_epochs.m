function [outstr, randstr] = APA_GrangerCausality_HPCACC_epochs(hpc_fname, acc_fname)

twin = 5; % seconds, pre/post entrance or avoidance

% load("D:\GarrettBlair\APA\HPCACC24500\processed_files\2023_06_17_H13_25_13_TR14_@placecells_ACC_miniscope2.mat")
load(acc_fname)
params.skip_ensemble = true;
[ms.room.momentary_pos_info,  rtemp]  = Fenton_ipos(ms, params.ipos_int_time, 'room', params);
[ms.arena.momentary_pos_info, ~]  = Fenton_ipos(ms, params.ipos_int_time, 'arena', params);
a_time = rtemp.t;

ipos_a = abs(ms.room.momentary_pos_info) - abs(ms.arena.momentary_pos_info);
ipos_a = nanmean(ipos_a,1);

load(hpc_fname)
params.skip_ensemble = true;
[ms.room.momentary_pos_info,  rtemp]  = Fenton_ipos(ms, params.ipos_int_time, 'room', params);
[ms.arena.momentary_pos_info, atemp]  = Fenton_ipos(ms, params.ipos_int_time, 'arena', params);
h_time = rtemp.t;

ipos_h = abs(ms.room.momentary_pos_info) - abs(ms.arena.momentary_pos_info);
ipos_h = nanmean(ipos_h,1);

if length(ipos_a)<length(ipos_h) % interp to the smaller vector
    ipos_h = interp1(h_time, ipos_h, a_time, 'linear');
elseif length(ipos_a)>length(ipos_h)
    ipos_a = interp1(a_time, ipos_a, h_time, 'linear');
end
ipos_a(isnan(ipos_a)) = interp1(find(~isnan(ipos_a)), ipos_a(~isnan(ipos_a)), find(isnan(ipos_a)), 'linear');
ipos_h(isnan(ipos_h)) = interp1(find(~isnan(ipos_h)), ipos_h(~isnan(ipos_h)), find(isnan(ipos_h)), 'linear');

figure; plot(ipos_a); hold on; plot(ipos_h)
%%
% r = ms.room;
r = rtemp;


t = rtemp.t;%ms.timestamps;
% bt = average_spks_time(t', params.ipos_int_time, t./1000, false, 'mean')
e = ms.room.entranceTimes./1000;



[th,rho] = cart2pol(rtemp.x,rtemp.y);
% approach_min_time_sec = 3;
ZONE_CENTER = pi/2;
ZONE_SIZE = pi/6;% distance from zone center to be counted appr_start in the zone
ZONE_NEARTHRESH = ZONE_SIZE + (pi-ZONE_SIZE)/2;% distance from zone center to be counted appr_start an approach
[ang_dist, ~] = angular_distance(th, ZONE_CENTER, 0);
in_zone = ang_dist<=ZONE_SIZE;
near    = ang_dist<=ZONE_NEARTHRESH & ~in_zone;
far     = ang_dist> ZONE_NEARTHRESH & ~in_zone;

appr = [];
entr = [];

appr.start = find(far(1:end-1)==1 & near(2:end)==1); % find the transitions between far and near
appr.end = find(near(1:end-1)==1 & far(2:end)==1); % find the transitions between near and far
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
figure(89); clf; hold on
axis([-40 40 -40 40])
ipos_a_sub = [];
ipos_h_sub = [];
ipos_a_entr = [];
ipos_h_entr = [];
ipos_a_peth = [];
ipos_h_peth = [];
ipos_a_entrpeth = [];
ipos_h_entrpeth = [];
for i = 1:length(appr.start)
    md = min(ang_dist(appr.start(i):appr.end(i)));
    td = ismember(t, t(appr.start(i):appr.end(i)));
    isentr = (e>=t(appr.start(i)) & e<=t(appr.end(i)));
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
        ipos_h_sub = cat(2, ipos_h_sub, ipos_h(inds));
        ipos_a_sub = cat(2, ipos_a_sub, ipos_a(inds));
        ipos_h_peth = cat(1, ipos_h_peth, ipos_h(inds));
        ipos_a_peth = cat(1, ipos_a_peth, ipos_a(inds));

    if plotting ==true
        plot3(r.x(inds), r.y(inds), inds, 'k')
    end
    elseif any(isentr)
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
        ipos_h_entr = cat(2, ipos_h_entr, ipos_h(inds));
        ipos_a_entr = cat(2, ipos_a_entr, ipos_a(inds));
        ipos_h_entrpeth = cat(1, ipos_h_entrpeth, ipos_h(inds));
        ipos_a_entrpeth = cat(1, ipos_a_entrpeth, ipos_a(inds));
    if plotting ==true
        plot3(r.x(inds), r.y(inds), inds, 'r')
    end
    end
end


% appr.valid = appr.length>=approach_min_time_sec & ~appr.isentrance;
appr.valid = ~appr.isentrance;
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
X = NaN(2, size(ipos_a_peth,2), size(ipos_a_peth,1)); % source x samples x trials
X(1,:,:) = ipos_a_peth';
X(2,:,:) = ipos_h_peth';
fs = .25;
source_labels = {'ACC' 'HPC'};
% nsamples = size(ipos_h_peth,2);
nsamples = [];
nrand = 1000;
verbose = false;
plotting = true;

[outstr, randstr] = granger_causality_Seth2015(X, fs, source_labels, nsamples, nrand, verbose, plotting);
figure;

figure; 
% subplot(121);
hold on; 
np = sqrt(size(ipos_a_peth,1));
shadedErrorBar([], nanmean(ipos_a_peth,1), nanstd(ipos_a_peth,[],1)./np,'lineProps', 'm-'); 
shadedErrorBar([], nanmean(ipos_h_peth,1), nanstd(ipos_h_peth,[],1)./np,'lineProps', 'c-'); 
drawnow()



 