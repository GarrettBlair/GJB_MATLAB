function [avoid_ind, is_entrance, avoid_prom, avoid_width, d, shk_ind, entr_ind, approach_ind] = shock_zone_dist_peakfinding(refStruct, t, shock_zone_center, shock_zone_size, distance_entrance_size, plotting);
% clear
% load('2023_01_13_H19_26_00_TR7_@placecells.mat')
%%
% shock_zone_center = pi/2; % typical room shock configuration
% shock_zone_size = pi/6; % size in rad from center to edge
% distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure

x = refStruct.x; 
y = refStruct.y; 
if median(abs(diff(t)))>1
    warning('time vec should be in seconds, check')
end
% t = ms.timestamps./1000;

[th, rth] = cart2pol(x,y);

% [edge_dist, inside_SZ] = angular_distance(th, shock_zone_center, shock_zone_size);
% [centerdist, ~] = angular_distance(th, shock_zone_center, 0);


d = shock_zone_center - th;
d = abs(mod(d + pi, 2*pi) - pi);
dd = d*-1;
peaksparams = {'MinPeakProminence',pi/8, 'MinPeakWidth', 1, 'MinPeakHeight', -1*distance_entrance_size};
if plotting==true
% figure(1); clf; polarhistogram(th,24);
figure(289); clf; set(gcf, 'Position', [ 50 600 1200 350])
% findpeaks(dd,t,peaksparams{:}, 'Annotate', 'extents');
subplot(2,1,1);
hold on
findpeaks(dd,t,peaksparams{:});
plot([t(1) t(end)], [shock_zone_size shock_zone_size], 'r:')
plot([t(1) t(end)], -1*[shock_zone_size shock_zone_size], 'r:')
plot([t(1) t(end)], [distance_entrance_size distance_entrance_size], 'g:')
plot(refStruct.shockTimes./1000, refStruct.shockTimes./1000*0 + 1, 'rx')
plot(refStruct.entranceTimes./1000, refStruct.entranceTimes./1000*0 + 1, 'md')
ylim([-4 1.5])
legend off
end
[pv, peak_time, peak_width, peak_prom] = findpeaks(dd,t, peaksparams{:});

valid_avoid = pv<shock_zone_size*0;%-.9
pv_shk = pv(~valid_avoid);
peak_time_shk = peak_time(~valid_avoid);
peak_width_shk = peak_width(~valid_avoid);
peak_prom_shk = peak_prom(~valid_avoid);

pv = pv(valid_avoid);
peak_time = peak_time(valid_avoid);
peak_width = peak_width(valid_avoid);
peak_prom = peak_prom(valid_avoid);
is_entrance = false(length(pv),1);

e = refStruct.shockTimes./1000;
shk_ind = NaN(size(e));
for i = 1:length(e)
    ei = e(i);
    shk_ind(i) = find(min(abs(ei-t)) == (abs(ei-t)), 1, 'first'); % get the index of shock times
    td = abs(ei-peak_time);
    if min(td)<3
%         tdind = find( min(td) == td , 1);
        is_entrance(td<3) = true;
%         peak_time(tdind) = ei;
    else
%         warning('no avoid near entr')
    end
    
end
e = refStruct.entranceTimes./1000;
entr_ind = NaN(size(e));
for i = 1:length(e)
    entr_ind(i) = find(min(abs(ei-t)) == (abs(ei-t)), 1, 'first'); % get the index of entrance times
end
avoid_pv = pv(~is_entrance);
avoid_time = peak_time(~is_entrance);
avoid_width = peak_width(~is_entrance);
avoid_prom = peak_prom(~is_entrance);

avoid_ind = find(ismember(t, avoid_time));
% shk_ind = find(ismember(t, e));
approach_ind = find(ismember(t, peak_time));
if plotting==true
% figure(1); clf; polarhistogram(th,24);
% scatter(peak_time(~is_entrance), pv(~is_entrance)+.1, 150, 'g.')
% scatter(peak_time(is_entrance),  pv(is_entrance)+.1,  150, 'r.')

subplot(2,1,2); cla
hold on
plot(t, -1*d, '-k')
plot([t(1) t(end)], [shock_zone_size shock_zone_size], 'r:')
plot([t(1) t(end)], -1*[shock_zone_size shock_zone_size], 'r:')
plot([t(1) t(end)], [0 0], 'k:')
plot(refStruct.shockTimes./1000, refStruct.shockTimes./1000*0 + 1, 'rx')
plot(refStruct.entranceTimes./1000, refStruct.entranceTimes./1000*0 + 1, 'md')
legend off
axis tight
scatter(peak_time, pv+.1, 20, 'ro')
scatter(avoid_time, avoid_pv+.1, 150, 'g.')
ylim([-4 1.5])
% scatter(peak_time_shk, pv_shk+.1, 150, 'r.')
end
%%
