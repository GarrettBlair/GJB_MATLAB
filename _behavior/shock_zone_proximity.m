function [avoid_ind, shks_ind] = shock_zone_proximity(ms, shock_zone_center, shock_zone_size, distance_entrance_size, plotting);
% clear
% load('2023_01_13_H19_26_00_TR7_@placecells.mat')
%%
% shock_zone_center = pi/2; % typical room shock configuration
% shock_zone_size = pi/6; % size in rad from center to edge
% distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure

x = ms.room.x; 
y = ms.room.y; 
t = ms.timestamps./1000;

[th, rth] = cart2pol(x,y);

[edge_dist, inside_SZ] = angular_distance(th, shock_zone_center, shock_zone_size);
[centerdist, ~] = angular_distance(th, shock_zone_center, 0);


d = shock_zone_center - th;
d = abs(mod(d + pi, 2*pi) - pi);
enter_ind = find( d(1:end-1)>distance_entrance_size & d(2:end)<=distance_entrance_size );
exit_ind = enter_ind*0;
ent_nframes  = enter_ind*0;
ent_dur  = enter_ind*0;
mindist_ind  = enter_ind*0;
was_shkd  = enter_ind*0==1;
exit_temp = find( d(1:end-1)<=distance_entrance_size & d(2:end)>distance_entrance_size );
if plotting==true
% figure(1); clf; polarhistogram(th,24);
figure; 
plot3(t,x,y); hold on;
plot3(t(inside_SZ), x(inside_SZ), y(inside_SZ), 'r.')
end
for j = 1:length(enter_ind)
    temp = exit_temp(exit_temp>enter_ind(j));
    
    if ~isempty(temp)
        exit_ind(j) = temp(1);
        inds = enter_ind(j):exit_ind(j);
%         dd = abs(centerdist(inds));
%         dd(find(inside_SZ(inds)==1,1):end) = NaN;
        dd = abs(edge_dist(inds));
        dd(find(inside_SZ(inds)==1,1):end) = NaN;
%         cd = centerdist(inds) <= shock_zone_size*2;
%         dd(cd==true) = NaN;
        mindist_ind(j) = enter_ind(j) + find(min(dd) == dd);
        ent_nframes(j) = exit_ind(j) - enter_ind(j);
        ent_dur(j) = t(exit_ind(j)) - t(enter_ind(j));
        was_shkd(j) = any(ms.room.entranceTimes./1000 >= t(enter_ind(j)) & ms.room.entranceTimes./1000 <= t(exit_ind(j)));
    else
        exit_ind(j) = NaN;
        ent_nframes(j) = NaN;
        ent_dur(j) = NaN;
        was_shkd(j) = 0;
    end
    if ent_dur(j)>1
%         plot3(t(inds), x(inds), y(inds), 'm-')
%         figure(2); clf;
%         hold on
%         plot(t, d, 'b')
%         plot([t(1) t(end)], [shock_zone_size shock_zone_size], 'r:')
%         plot([t(1) t(end)], -1*[shock_zone_size shock_zone_size], 'r:')
%         plot([t(1) t(end)], [distance_entrance_size distance_entrance_size], 'g:')
%         plot(ms.room.entranceTimes./1000, ms.room.entranceTimes./1000*0, 'md')
%         avoids = ~was_shkd & ent_dur>2 & ent_nframes > 30;
%         avoid_ind = mindist_ind(avoids);
%         shks = was_shkd & ent_dur>2 & ent_nframes > 30;
%         shks_ind = mindist_ind(shks);
%         plot(t(avoid_ind), d(avoid_ind), 'g.')
%         plot(t(shks_ind),  d(shks_ind), 'r.')
%         drawnow
    end
end
avoids = ~was_shkd & ent_dur>2 & ent_nframes > 30;
avoid_ind = mindist_ind(avoids);
shks = was_shkd & ent_dur>2 & ent_nframes > 30;
shks_ind = mindist_ind(shks);

% figure; 
% plot3(t,x,y); hold on;
% plot3(t(avoid_ind), x(avoid_ind), y(avoid_ind), 'go')
% plot3(t(shks_ind), x(shks_ind), y(shks_ind), 'ro')

% d1 = pi/3 - th;
% d1 = abs(mod(d1 + pi, 2*pi) - pi);
% d2 = 2*pi/3 - th;
% d2 = abs(mod(d2 + pi, 2*pi) - pi);
% d = min(d1, d2);
% d = d./(pi-pi/6);
%
if plotting == true
figure(2); clf;
% subplot(211); hold on; 
% plot(t,   th, 'k');
% plot([t(1) t(end)], shock_zone_center+shock_zone_size*[1 1], 'r:')
% plot([t(1) t(end)], shock_zone_center-shock_zone_size*[1 1], 'r:')
% subplot(212); cla; 
hold on
plot(t, d, 'b')
plot([t(1) t(end)], [shock_zone_size shock_zone_size], 'r:')
plot([t(1) t(end)], -1*[shock_zone_size shock_zone_size], 'r:')
plot([t(1) t(end)], [distance_entrance_size distance_entrance_size], 'g:')
plot(ms.room.entranceTimes./1000, ms.room.entranceTimes./1000*0, 'md')
% plot(ms.room.shockTimes./1000, ms.room.shockTimes./1000*0, 'ro')

plot(t(avoid_ind), d(avoid_ind), 'g.')
plot(t(shks_ind),  d(shks_ind), 'r.')
% plot(t(enter_ind), d(enter_ind), 'g.')
% plot(t(exit_ind),  d(enter_ind), 'r.')
% plot(ms.room.shockTimes./1000, ms.room.shockTimes./1000*0, 'ro')

% ylim([-.1 3.25])
set(gca, 'YDir', 'reverse')
end
