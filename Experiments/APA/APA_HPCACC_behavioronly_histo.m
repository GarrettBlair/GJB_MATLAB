clear
%%
ddir = 'D:\GarrettBlair\APA\TF_behavonly\';
day0 = datetime(2024,3,26);
aname = 'TF';
anums = 1:5;

max_tr = 15;
params = [];
APA_rat_imaging_params_current;

n = length(anums);
num_entr    = NaN(n, max_tr);
first_entr  = NaN(n, max_tr);
mta         = NaN(n, max_tr);
dist        = NaN(n, max_tr);
x           = NaN(n, max_tr);
trn         = NaN(n, max_tr);
sessname    = cell(n, max_tr);

for a = 1:n
    for s = 0:max_tr
        r_fn = [ddir aname num2str(anums(a)) '_TR' num2str(s) '_Room.dat'];
        a_fn = [ddir aname num2str(anums(a)) '_TR' num2str(s) '_Arena.dat'];
        r_ret = [ddir aname num2str(anums(a)) '_RET' num2str(s) '_Room.dat'];
        a_ret = [ddir aname num2str(anums(a)) '_RET' num2str(s) '_Arena.dat'];
        sname = [];
        if isfile(r_ret) && isfile(a_ret)
            r_fn = r_ret;
            a_fn = a_ret;
            sname = ['RET' num2str(s)];
        elseif isfile(r_fn) && isfile(a_fn)
            sname = ['TR' num2str(s)];
        end
        sessname{a, s+1}    = sname;
        if ~isempty(sname)
            disp(r_fn)
            [room, arena, params_sub] = behavior_DAT_tracking_eval(r_fn, a_fn, params);
            x(a, s+1) = floor(days( room.DAT_fileinfo.datetime - day0 ));
            trn(a, s+1) = s;
            dist(a, s+1) = sum(arena.speed.*arena.dt)/100; % total meters
            sessmins = (room.timestamps(end)/60000);
            if room.num_entrances == 0
                num_entr(a, s+1) = room.num_entrances/sessmins;
                max_time = room.timestamps(end)/1000 + 30;
                first_entr(a, s+1) = max_time;
                mta(a, s+1) = max_time;
            else
                num_entr(a, s+1) = room.num_entrances/sessmins;
                first_entr(a, s+1) = room.first_entrance;
                et = room.timestamps(room.entrance_start_idx);
                m = max(abs(diff([0; et; room.timestamps(end)])))/1000;
                mta(a, s+1) = m;
            end
        end
    end
end
%%
figure(1); clf; hold on
max_days = max(max(floor(x)+1));
% TF1 TF4 TF5 will be remote
% TF2 TF3     will be recent
isremote = [1 0 0 1 1]; 
ismale   = [1 0 1 0 1];
s_style = {};
remote_color = [.1,1,.7; .3,.7,1]*.8;
linestyle = {'-', '--'};
ssize = 50;
lw = 1.5;
time_scale = 'log';
% l_style = {'LineWidth', .25, 'LineColor'}
first_epm = NaN(n, max_days);
first_mta = NaN(n, max_days);
first_fe  = NaN(n, max_days);
first_x  = NaN(n, max_days);

for a = 1:n
d = floor(x(a,:));
ddiff = diff([-1, d])>0;

subplot(3,2,1); 
hold on
plot(trn(a,:), num_entr(a,:), 'Color', remote_color(isremote(a)+1,:), 'LineWidth', lw, 'LineStyle', linestyle{ismale(a)+1})

subplot(3,2,3); 
hold on
plot(trn(a,:), mta(a,:), 'Color', remote_color(isremote(a)+1,:), 'LineWidth', lw, 'LineStyle', linestyle{ismale(a)+1})

subplot(3,2,5); 
hold on
plot(trn(a,:), first_entr(a,:), 'Color', remote_color(isremote(a)+1,:), 'LineWidth', lw, 'LineStyle', linestyle{ismale(a)+1})

subplot(3,2,2); 
hold on
plot(x(a,ddiff), num_entr(a,ddiff), 'Color', remote_color(isremote(a)+1,:), 'LineWidth', lw, 'LineStyle', linestyle{ismale(a)+1})
scatter(x(a,ddiff), num_entr(a,ddiff), ssize, 'k.')

subplot(3,2,4); 
hold on
plot(x(a,ddiff), mta(a,ddiff), 'Color', remote_color(isremote(a)+1,:), 'LineWidth', lw, 'LineStyle', linestyle{ismale(a)+1})
scatter(x(a,ddiff), mta(a,ddiff), ssize, 'k.')

subplot(3,2,6); 
hold on
plot(x(a,ddiff), first_entr(a,ddiff), 'Color', remote_color(isremote(a)+1,:), 'LineWidth', lw, 'LineStyle', linestyle{ismale(a)+1})
scatter(x(a,ddiff), first_entr(a,ddiff), ssize, 'k.')

first_epm(a,1:sum(ddiff)) = num_entr(a,ddiff);
first_mta(a,1:sum(ddiff)) = mta(a,ddiff);
first_fe(a,1:sum(ddiff))  = first_entr(a,ddiff);
first_x(a,1:sum(ddiff))   = x(a,ddiff);
end

subplot(3,2,1); 
hold on
scatter(trn(:), num_entr(:), ssize, 'k.')
plot(nanmean(trn,1), nanmean(num_entr,1), 'k-', 'LineWidth', lw*2)
axis([-1 max_tr+1 -.5 4])

subplot(3,2,3); 
hold on
scatter(trn(:), mta(:), ssize, 'k.')
plot(nanmean(trn,1), nanmean(mta,1), 'k-', 'LineWidth', lw*2)
axis([-1 max_tr+1 0 1600])
set(gca, 'YScale', time_scale)

subplot(3,2,5); 
hold on
scatter(trn(:), first_entr(:), ssize, 'k.')
plot(nanmean(trn,1), nanmean(first_entr,1), 'k-', 'LineWidth', lw*2)
axis([-1 max_tr+1 0 1600])
set(gca, 'YScale', time_scale)

subplot(3,2,2); 
plot(nanmean(first_x,1), nanmean(first_epm,1), 'k-', 'LineWidth', lw*2)
axis([-1 max_days+1 -.5 4])

subplot(3,2,4); 
plot(nanmean(first_x,1), nanmean(first_mta,1), 'k-', 'LineWidth', lw*2)
axis([-1 max_days+1 0 1600])
set(gca, 'YScale', time_scale)

subplot(3,2,6); 
plot(nanmean(first_x,1), nanmean(first_fe,1), 'k-', 'LineWidth', lw*2)
axis([-1 max_days+1 0 1600])
set(gca, 'YScale', time_scale)

%
figure(2); clf
first_dist = NaN(n, max_days);
for a = 1:n
d = floor(x(a,:));
ddiff = diff([-1, d])>0;

subplot(1,2,1); 
hold on
plot(trn(a,:), dist(a,:), 'Color', remote_color(isremote(a)+1,:), 'LineWidth', lw, 'LineStyle', linestyle{ismale(a)+1})

subplot(1,2,2); 
hold on
plot(x(a,ddiff), dist(a,ddiff), 'Color', remote_color(isremote(a)+1,:), 'LineWidth', lw, 'LineStyle', linestyle{ismale(a)+1})
scatter(x(a,ddiff), dist(a,ddiff), ssize, 'k.')
first_dist(a,1:sum(ddiff))  = dist(a,ddiff);
end

subplot(1,2,1); 
hold on
scatter(trn(:), dist(:), ssize, 'k.')
plot(nanmean(trn,1), nanmean(dist,1), 'k-', 'LineWidth', lw*2)
axis([-1 max_tr+1 -1 150])

subplot(1,2,2); 
plot(nanmean(first_x,1), nanmean(first_dist,1), 'k-', 'LineWidth', lw*2)
axis([-1 max_days+1 -1 150])
