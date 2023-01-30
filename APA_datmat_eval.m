% % % function [entrances, numShocks, path_length] = APA_datmat_eval(datmatFolder, animals, params)
%%
datmatFolder = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\DAT_files\matlab files';
% fnames = dir(datmatFolder);
% fnames = fnames(3:end); % remove the '.' and '..' entries
% sess = {'HAB0', 'HAB1', 'HAB2', 'HAB3', 'TR0', 'TR1', 'TR2', 'TR3', 'X0', 'X1'};
sess = {'TR0' 'TR1' 'TR2' 'TR3' 'X0' 'X1' 'X2' 'X3' 'TR4' 'TR5' 'TR6' 'TR7' 'RET0' 'NEW0' 'NEW1' 'NEW2' 'NEW3' 'NEW4'};
animals = {'Acc20832', 'Acc19947'};

nsess = length(sess);
nsubs = length(animals);
entrances = NaN(nsubs, nsess);
numShocks = NaN(nsubs, nsess);
shockOcc = NaN(nsubs, nsess);
path_length = NaN(nsubs, nsess);
entrances_m = NaN(nsubs, nsess);
numShocks_m = NaN(nsubs, nsess);
path_length_m = NaN(nsubs, nsess);

params.arena_radius             = 40; % in cm
params.arena_center             = [ 127.5, 127.5]; % pixel center of behav cam, [x,y]
params.pixpercm                 = 3.1220; % pixel radius of behav cam env
params.behav_fps                = 30;
params.behav_smoothing_interval = .25; % in seconds, length of smoothing kernel
params.nan_interp               = true;
params.correct_dt               = true; % correct for large jumps in timestamp file when constructing vmap

figure(94); clf; figure(95); clf
for a = 1:nsubs
    for s = 1:nsess
        %%
        arena_tracking_fname = [datmatFolder '\' animals{a} '_' sess{s} '_Arena.mat'];
        room_tracking_fname  = [datmatFolder '\' animals{a} '_' sess{s} '_Room.mat'];
        if isfile(arena_tracking_fname) && isfile(room_tracking_fname)
        [room, arena, ~] = behavior_DAT_tracking_eval(room_tracking_fname, arena_tracking_fname, params);
%         [arena, ~] = read_APA_csv(arena_tracking_fname, []);
%         [room , ~] = read_APA_csv(room_tracking_fname, []);
        aaa = room.pol_theta >= pi/3 & room.pol_theta < 2*pi/3;
        bbb = sum(aaa)/length(aaa);

        path_length(a,s) = sum(arena.speed.*arena.dt)./100; % path length in meters
        entrances(a,s) = room.num_entrances + arena.num_entrances;
        numShocks(a,s) = room.num_shocks + arena.num_shocks;
        shockOcc(a,s) = bbb;
        mins = room.timestamps(end)/(1000*60);
        path_length_m(a,s) = path_length(a,s)/mins; % path length in meters per min
        entrances_m(a,s) = entrances(a,s)/mins; % per min
        numShocks_m(a,s) = numShocks(a,s)/mins; % per min
%         path_length_m(a,s) = path_length(a,s)/path_length(a,s);% pre meter    mins; % path length in meters per min
%         entrances_m(a,s) = entrances(a,s)/path_length(a,s);% per meter    mins; % per min
%         numShocks_m(a,s) = numShocks(a,s)/path_length(a,s);% per meter    mins; % per min

%         figure(94);
%         subplot_tight(nsubs, nsess, nsess*(a-1) + s)
%         polarhistogram(room.pol_theta, [0:pi/12:2*pi],'FaceColor','red','FaceAlpha',.3)
%         title(sprintf('%s\n%s Room', animals{a}, sess{s}))
%         axis off
% 
%         figure(95);
%         subplot_tight(nsubs, nsess, nsess*(a-1) + s)
%         polarhistogram(arena.pol_theta, [0:pi/12:2*pi],'FaceColor','blue','FaceAlpha',.3)
%         title(sprintf('%s\n%s Arena', animals{a}, sess{s}))
%         axis off
        else
            fprintf('\nFile not found: %s', room_tracking_fname)
            fprintf('\nFile not found: %s', arena_tracking_fname)
            fprintf('\n')
        end
        
    end
end

%%

figure(1); clf;
subplot(1,3,1);
hold on;
plot(1:nsess, entrances_m')
% axis([0 nsess+1 0 .75])
% plot(1:nsess, entrances')
% axis([0 nsess+1 0 90])
title('Entrances / min')
set(gca, 'XTick', 1:nsess, 'XTickLabels', sess, 'XTickLabelRotation', 90, 'YAxisLocation', 'right')

subplot(1,3,2);
hold on;
plot(1:nsess, numShocks_m')
% axis([0 nsess+1 0 1.5])
% plot(1:nsess, numShocks')
% axis([0 nsess+1 0 150])
title('Shock / min')
set(gca, 'XTick', 1:nsess, 'XTickLabels', sess, 'XTickLabelRotation', 90, 'YAxisLocation', 'right')

subplot(1,3,3);
hold on;
% plot(1:nsess, shockOcc')
plot(1:nsess, path_length')
% axis([0 nsess+1 0 250])
title('Meters ran')
set(gca, 'XTick', 1:nsess, 'XTickLabels', sess, 'XTickLabelRotation', 90, 'YAxisLocation', 'right')
legend(animals)


% % % end
