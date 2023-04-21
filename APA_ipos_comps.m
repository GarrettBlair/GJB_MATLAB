dataDir = 'D:\APA recordings\';
animals = {'Acc20832' 'Acc19947' 'Hipp18240'};% 'Hipp16942'};
animals = {'Acc20832' 'Acc19947' 'Hipp18240' 'Hipp16942'};
if false % adding ipos data to files
for aLoop = 1:length(animals)
    cd([dataDir animals{aLoop} '\processed_files\']);
    temp = dir('*mat');
    nFiles = length(temp);
    fprintf('\n\tANIMAL: %s', animals{aLoop})
    for i = 1:nFiles
        tempdata = load([temp(i).name]);
        if ~isfield(tempdata, 'iposData')
            fprintf('\nIPOS: %s', [temp(i).name])
            [room_ipos, sub_room, room_ensemble]    = Fenton_ipos(tempdata.ms, .25, 'room', []);
            [arena_ipos, sub_arena, arena_ensemble] = Fenton_ipos(tempdata.ms, .25, 'arena', []);
            iposData.time_res = .25;
            iposData.room = sub_room;
            iposData.arena = sub_arena;
            iposData.room_info = room_ipos;
            iposData.arena_info = arena_ipos;
            iposData.room_ensemble = room_ensemble;
            iposData.arena_ensemble = arena_ensemble;
            iposData.ipos = abs(room_ipos) - abs(arena_ipos);
            save([temp(i).name], 'iposData', '-append');
        else
            fprintf('\nSKIP: %s', [temp(i).name])
        end
    end
end
end
if false % adding shock and entrance times to ms
animals = {'Hipp18240'};
for aLoop = 1:length(animals)
    cd([dataDir animals{aLoop} '\processed_files\']);
    temp = dir('*mat');
    nFiles = length(temp);
    fprintf('\n\tANIMAL: %s', animals{aLoop})
    for i = 1:nFiles
        tempdata = load([temp(i).name]);
        ms = tempdata.ms;
        if ~isfield(ms.room, 'entranceTimes')
            sess = tempdata.AnimalDir.SessType{i};
            Arena_files = dir(['D:\APA recordings\DAT_files\Hipp18240_' sess '_Arena.dat']);
            arena_tracking_fname = [Arena_files.folder '\' Arena_files.name];
            Room_files  = dir(['D:\APA recordings\DAT_files\Hipp18240_' sess '_Room.dat']);
            room_tracking_fname = [Room_files.folder '\' Room_files.name];
            [ms, room, arena, params_sub] = behavior_DAT_tracking_file(ms, room_tracking_fname, arena_tracking_fname, []);
            save([temp(i).name], 'ms', '-append');
        else
            fprintf('\nSKIP: %s', [temp(i).name])
        end
    end
end
end
%% COLLECT DATA
dataDir = 'D:\APA recordings\';
animals = {'Acc20832' 'Acc19947' 'Hipp18240'};% 'Hipp16942'};
sfep = NaN(30, 3);
n_e = NaN(30, 3);
ipos_dist_corr = NaN(30, 3);
iposE_dist_corr = NaN(30, 3);
ensemble_sfep = NaN(30, 3);
corrbins = [-.2:.025:1]; corrcenter = corrbins(1:end-1) + mean(abs(diff(corrbins)));
ipos_corr_hists = zeros(length(corrcenter), 3);
for aLoop = 1:length(animals)
    cd([dataDir animals{aLoop} '\processed_files\']);
    temp = dir('*mat');
    nFiles = length(temp);
    fprintf('\n\tANIMAL: %s', animals{aLoop})
    for i = 1:nFiles
        tempdata = load([temp(i).name]);
        if isfield(tempdata, 'iposData')
            ip = tempdata.iposData;
            [a_theta, a_rho] = cart2pol(ip.room.x, ip.room.y);
            [a_dist] = abs(angular_distance(a_theta, pi/2));
            
            n_e(i,aLoop) = length(tempdata.ms.room.entranceTimes);
            d_ipos = nanmean(ip.ipos,1);
            valid = ~isnan(d_ipos);
            [seg_iposcor, p] = corr(d_ipos(valid)', ip.ipos(:, valid)');
            ipos_corr_hists(:, aLoop) = ipos_corr_hists(:, aLoop) + histcounts(seg_iposcor, corrbins)';
            
            nanind = isnan(d_ipos.*a_dist);
            ipos_dist_corr(i,aLoop) = corr(d_ipos(~nanind)', a_dist(~nanind)');
            sfep(i,aLoop) = sum(d_ipos>0) / (sum(d_ipos>0)+ sum(d_ipos<0));
            
            ensemble_ipos = (abs(ip.room_ensemble) - abs(ip.arena_ensemble))';
            ensemble_sfep(i,aLoop) = sum(ensemble_ipos>0) / (sum(ensemble_ipos>0)+ sum(ensemble_ipos<0));
            nanind = isnan(ensemble_ipos.*a_dist);
            iposE_dist_corr(i,aLoop) = corr(ensemble_ipos(~nanind)', a_dist(~nanind)');
%             fprintf('\nsfep %2.4f  -  %s', sfep, [temp(i).name])
            
        else
            fprintf('\nNO IPOS: %s', [temp(i).name])
        end
    end
end
%% PLOTTING
figure(4); clf
figure(5); clf
figure(7); clf
for aLoop = 1:3
    cd([dataDir animals{aLoop} '\processed_files\']);
    temp = dir('*mat');
    nFiles = length(temp);
    nam = {};
    for i = 1:nFiles
        fname = temp(i).name;
        a = strfind(fname, '_');
        nam{i} = fname(a(6)+1:a(7)-1);
    end
    fprintf('\n\tANIMAL: %s', animals{aLoop})
    tempdata = load([temp(1).name]);
    figure(4); subplot(3,1,aLoop)
    hold on
    plot([1 nFiles], [.5 .5], 'k:')
    plot(1:nFiles, sfep(1:nFiles, aLoop), 'k-o')
    ylabel('SFEP')
    set(gca, 'XTick', 1:nFiles, 'XTickLabel', nam, 'XTickLabelRotation', -90)
    ylim([.2 .7])
    xlim([0 nFiles+1])
    %     ylim([0 1])
    yyaxis('right')
    set(gca, 'YColor', 'r')
    ylabel('Entrances')
    plot(1:nFiles, n_e(1:nFiles, aLoop), 'r.-')
    ylim([-10 80])
    drawnow

    figure(5); subplot(3,1,aLoop)
    hold on
    plot(1:nFiles, ipos_dist_corr(1:nFiles, aLoop), 'k.-')
    ylabel('Ipos corr avoidance')
    set(gca, 'XTick', 1:nFiles, 'XTickLabel', nam, 'XTickLabelRotation', -90)
    ylim([-.7 .7])
    xlim([0 nFiles+1])
    %     ylim([0 1])
    yyaxis('right')
    set(gca, 'YColor', 'r')
    ylabel('Entrances')
    plot(1:nFiles, n_e(1:nFiles, aLoop), 'r.-')
    ylim([-10 80])
    drawnow
    
    figure(7);
    hold on
    counttrace = ipos_corr_hists(:,aLoop)./sum(ipos_corr_hists(:,aLoop));
    plot(corrcenter, cumsum(counttrace), 'LineWidth', 2)
end
    figure(7);
legend({'ACC1' 'ACC2' 'HIPP1'})
entr1 = n_e(:);
valid = entr1>0; entr1 = entr1(valid);
sfep1 = sfep(:); sfep1 = sfep1(valid);
ipos_dist_corr1 = ipos_dist_corr(:); ipos_dist_corr1 = ipos_dist_corr1(valid);
iposE_dist_corr1 = iposE_dist_corr(:); iposE_dist_corr1 = iposE_dist_corr1(valid);

figure(6); clf; hold on;
scatter(n_e(:), sfep(:), 'ko')
scatter(n_e(:), ipos_dist_corr(:), 'go')
scatter(n_e(:), iposE_dist_corr(:), 'bo')
axis([ 0 100 -1 1])
[h,p_ipos] = corr(entr1(~isnan(ipos_dist_corr1)), ipos_dist_corr1(~isnan(ipos_dist_corr1)));
[h,p_ipose] = corr(entr1(~isnan(iposE_dist_corr1)), iposE_dist_corr1(~isnan(iposE_dist_corr1)));
[h,p_sfep] = corr(entr1(~isnan(sfep1)), sfep1(~isnan(sfep1)));
text(85, .5, sprintf('sfep\n%0.5f', p_sfep), 'Color', 'k')
text(85, .2, sprintf('ensem ipos\n%0.5f', p_ipose), 'Color', 'b')
text(85, -.2, sprintf('ipos corr\n%0.5f', p_ipos), 'Color', 'g')
% scatter(n_e(:), ensemble_sfep(:), 'ro')
% scatter(sfep(:), ensemble_sfep(:), 'bo')
%%
testset = load("D:\APA recordings\Acc20832\processed_files\2023_01_13_H18_48_32_TR6_@placecells.mat");
int_time = .25; % seconds
[a_theta, a_rho] = cart2pol(sub_room.x, sub_room.y);
[a_dist] = abs(angular_distance(a_theta, pi/2));
[room_pos, sub_room, ensemble_room] = Fenton_ipos(testset.ms, int_time, 'room', []);
[arena_pos, sub_arena, ensemble_arena] = Fenton_ipos(testset.ms, int_time, 'arena', []);
e_pos = abs(ensemble_room)-abs(ensemble_arena);
i_pos = nanmean(abs(room_pos)-abs(arena_pos), 1);
% create a data set to test with varying number of outliers and outlier
% distance
figure; plot(abs(ensemble_room)); hold on; plot(abs(ensemble_arena))
figure(27); clf; hold on
plot(e_pos, 'b')
% ylim([-8 8])
yyaxis('right')
plot(i_pos, 'k')
% ylim([-4 4])

bins = [-45, -36:4:36, 45];
[emap]   = make_occupancymap_2D(sub_room.x, sub_room.y, e_pos, bins, bins);
[imap]   = make_occupancymap_2D(sub_room.x, sub_room.y, i_pos, bins, bins);


