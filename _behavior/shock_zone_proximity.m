% clear
% load('2023_01_13_H19_26_00_TR7_@placecells.mat')
%%
shock_zone_center = pi/2; % typical room shock configuration
shock_zone_size = pi/6; % size in rad from center to edge
distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure

x = ms.room.x; y = ms.room.y; t = ms.timestamps./1000;
[th, rth] = cart2pol(x,y);
figure(1); clf; polarhistogram(th,24);

d = shock_zone_center - th;
d = abs(mod(d + pi, 2*pi) - pi);
enter_ind = find( d(1:end-1)>distance_entrance_size & d(2:end)<=distance_entrance_size );
exit_ind = enter_ind*0;
ent_dur  = enter_ind*0;
was_shkd  = enter_ind*0==1;
exit_temp = find( d(1:end-1)<=distance_entrance_size & d(2:end)>distance_entrance_size );
for j = 1:length(enter_ind)
    temp = exit_temp(exit_temp>enter_ind(j));
    exit_ind(j) = temp(1);
    ent_dur(j) = exit_ind(j) - enter_ind(j);
    was_shkd(j) = any(ms.room.entranceTimes./1000 >= t(enter_ind(j)) & ms.room.entranceTimes./1000 <= t(exit_ind(j)));
end
% d1 = pi/3 - th;
% d1 = abs(mod(d1 + pi, 2*pi) - pi);
% d2 = 2*pi/3 - th;
% d2 = abs(mod(d2 + pi, 2*pi) - pi);
% d = min(d1, d2);
% d = d./(pi-pi/6);

%%
figure(2); clf;
subplot(211); hold on; 
plot(t,   th, 'k');
plot([t(1) t(end)], shock_zone_center+shock_zone_size*[1 1], 'r:')
plot([t(1) t(end)], shock_zone_center-shock_zone_size*[1 1], 'r:')
subplot(212); cla; hold on
plot(t, d, 'b')
plot([t(1) t(end)], [shock_zone_size shock_zone_size], 'r:')
plot([t(1) t(end)], -1*[shock_zone_size shock_zone_size], 'r:')
plot([t(1) t(end)], [distance_entrance_size distance_entrance_size], 'g:')
plot(ms.room.entranceTimes./1000, ms.room.entranceTimes./1000*0, 'md')
plot(ms.room.shockTimes./1000, ms.room.shockTimes./1000*0, 'ro')

plot(t(enter_ind), d(enter_ind), 'g.')
plot(t(exit_ind),  d(enter_ind), 'r.')
plot(ms.room.shockTimes./1000, ms.room.shockTimes./1000*0, 'ro')

% ylim([-.1 3.25])
set(gca, 'YDir', 'reverse')
%%
dub = 0*rth + 20;
% d = sqrt(rth.^2 + dub.^2 - 2.*rth.*dub.*cos(th- pi/2));

% d = pi/2 - th;
% d = abs(mod(d + pi, 2*pi) - pi)./pi;



figure(3); clf
inds = 1:1500;
for i = 1:2:length(inds)%length(x)
    %%
    figure(3); clf;
    hold on
    plot(x(inds(1):inds(i)), y(inds(1):inds(i)), 'k-')
    plot(x(inds(i)), y(inds(i)), 'ro')
    axis([-42 42 -42 42])
    title(sprintf('%3.1f sec', t(inds(i))))
    figure(4); clf;
    polarplot(d(inds(1):inds(i)), rth(inds(1):inds(i)), 'k-')
    rlim([0 42])
drawnow
pause(.01)
end