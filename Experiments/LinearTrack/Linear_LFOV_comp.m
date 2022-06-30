clear
%%
ddir = 'C:\Users\gjb326\Desktop\Sample Data\Linear Maze Data';
sess2use = 1:10;
temp = dir(ddir);
%
fprintf('\n\t~~~DATA:~~~\n')
filecount = 0;
allfiles = cell(10,1);
sessnum = cell(10,1);
allnames = cell(10,1);
isLFOV = cell(10,1);
num_animals = 0;
a_list = {};
sess_list = {};
for i = 1:size(temp,1)
    aname = temp(i).name;
    %     d = strcat(temp(i).folder, '\', aname, '\Session2\', aname, '_Sess2.mat');
    for j = 1:length(sess2use)
        d = sprintf('%s\\%s\\Session%d\\%s_Sess%d.mat', temp(i).folder, aname, sess2use(j), aname, sess2use(j));
        %     d_alt = dir(strcat(temp(i).folder, '\', aname, '\Session2_H*'));
        d_alt = dir(sprintf('%s\\%s\\Session%d_H*', temp(i).folder, aname, sess2use(j)));
        if ~isempty(d_alt)
            d = sprintf('%s\\%s\\%s_Sess%d.mat', d_alt.folder, d_alt.name, aname, sess2use(j));
        end
        if exist(d, 'file')
            filecount = filecount +1;
            fprintf('%s\n', d)
            allfiles{filecount} = d;
            aa = strfind(d, '_Sess')+5;
            bb = strfind(d, '.mat')-1;
            sessnum{filecount} = str2double(d(aa:bb));

            allnames{filecount} = aname;
            isLFOV{filecount} = false;
            
            if num_animals>0
                animal_ind = strcmp(a_list, aname);
                    if any(animal_ind)==0
                        num_animals = num_animals+1;
                        a_list{num_animals} = aname;
                    else
                        try
                            sess_list{animal_ind} = cat(1, sess_list{animal_ind}, sess2use(j));
                        catch
                            sess_list{num_animals} = sess2use(j);
                        end
                    end
            else
                num_animals = num_animals+1;
                a_list{1} = aname;
            end

        end
    end
end
%%
ddir2 = 'C:\Users\gjb326\Desktop\Sample Data\MiniData Hipp Blair';
hippnums = [6 7 8 9 12 13 15 18 30 32 34 35 36];
end_nonlfov = filecount;

for i = 1:length(hippnums)
    aname = strcat('Hipp', num2str(hippnums(i)));
    %     d = strcat(ddir2, '\', aname, '\', aname, '_p\', aname, '_linear2.mat');
    for j = 1:length(sess2use)
        d = sprintf('%s\\%s\\%s_p\\%s_linear%d.mat', ddir2, aname, aname, aname, sess2use(j));
        if exist(d, 'file')
            filecount = filecount +1;
            fprintf('%s\n', d)
            allfiles{filecount} = d;
            aa = strfind(d, '_linear')+7;
            bb = strfind(d, '.mat')-1;
            sessnum{filecount} = str2double(d(aa:bb));
            allnames{filecount} = aname;
            isLFOV{filecount} = true;
        end
        
        if num_animals>0
            animal_ind = strcmp(a_list, aname);
            if any(animal_ind)==0
                num_animals = num_animals+1;
                a_list{num_animals} = aname;
                sess_list{num_animals} = sess2use(j);
            else
%                 try
                    sess_list{animal_ind} = cat(1, sess_list{animal_ind}, sess2use(j));
%                 catch
%                     sess_list{num_animals} = sess2use(j);
%                 end
            end
        else
            num_animals = num_animals+1;
            a_list{1} = aname;
        end
    end
end
%%
% A = struct();
for i = 1:filecount
    fprintf('%s... ', allfiles{i})
    temp = load(allfiles{i});
        tempname = []; % sprintf('%s.sess%d', allnames{i}, sessnum{i});
    if isLFOV{i}
        tempname.aname = allnames{i};
        tempname.sessnum = sessnum{i};
        tempname.x = temp.ms.x;
        tempname.y = temp.ms.y;
        tempname.t = (temp.ms.time - temp.ms.time(1))*1000;
        sr = length(temp.ms.SR.is_bee) + length(temp.ms.SL.is_bee);
        lr = length(temp.ms.LR.is_bee) + length(temp.ms.LL.is_bee);
        tempname.numr = sr+lr;
        temp.ms.timestamp_file;
        tempname.dual = true;
        tempname.isLFOV = true;
        tempname.sr = sr;
        tempname.lr = lr;
    else
        ind = find(~isnan(temp.xpos), 1, 'last');
        tempname.aname = allnames{i};
        tempname.x = temp.xpos(1:ind);
        tempname.y = temp.ypos(1:ind);
        tempname.t = temp.VT_time(1:ind) - temp.VT_time(1);
        ind = find(~isnan(temp.unique_zone), 1, 'last');
        tempname.unique_zone = temp.unique_zone(1:ind);
        ind = find(~isnan(temp.zone_seq), 1, 'last');
        tempname.zone_seq = temp.zone_seq(1:ind);
        tempname.isLFOV = false;
        try
            sr = temp.short_rewards;
            lr = temp.long_rewards;
            tempname.sr = sr;
            tempname.lr = lr;
            tempname.numr = sr+lr;
            tempname.dual = true;
        catch
            tempname.numr = temp.num_rewards;
            tempname.sr = tempname.numr;
            tempname.lr = 0;
            tempname.dual = false;
        end
    end
    eval(sprintf('A.animal%s.sess%d = tempname;', allnames{i}, sessnum{i}));
    fprintf('Done!\n')
end
save('C:\Users\gjb326\Desktop\Sample Data\Linear Maze Data\AllData.mat', 'A');
%%
controls = {'A80','A81','A82','A83','A84','A85','A86','A87','A88','A89','A90','A91','A92','FC1','FC2','FC5','FC7','FC9','FC10','FC11'};
implants ={'Hipp6'};
for i = 2:length(hippnums)
    implants{end+1} = sprintf('Hipp%d', hippnums(i));
end
figure(9); clf; 
rr = NaN(filecount,1);
switch sess2use
    case 1
behav = [29:39, 42:48];%1:55; % [10, 12, 14:25];38:50;% 
lfov = [56:67];
    case 2
behav = [26:50];%1:55; % [10, 12, 14:25];38:50;% 
lfov = [57:69];
    case 3
behav = 38:44;%26:44; % [10, 12, 14:25];38:50;% 
lfov = [50:62];
end
% behav = 1:55; % [10, 12, 14:25];38:50;% 
% lfov = [57:69];
nsess = length(controls) + length(controls))
all_an = [behav, lfov];
for i = 1:length([behav, lfov])
    %
    ind = all_an(i);
    a = A{ind}.numr;
    b = (A{ind}.t(end)/10^6)/60;
    c = 'k';
    lr = A{ind}.lr;
    if isLFOV{ind}
        c = 'b';
    end

%     scatter(i, A{ind}.lr, 15, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .5)
    if b<=21 && b>=7 && A{ind}.lr>0 % A{ind}.dual == true
    rr(ind) = a/b; % rewards per min
%     subplot_tight(1, 3, 1:2, [.2,.05])
    subplot(1, 3, 1:2)
    hold on
%     scatter(i, lr, 15, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .5)
    scatter(i, rr(ind), 15, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .5)
%     h = text(ind, rr(ind)*1.2, A{ind}.aname, 'Interpreter', 'latex');
    end
end
set(gca, 'XTick', 1:length([behav, lfov]), 'XTickLabel', allnames(all_an), 'XTickLabelRotation', 90)
axis([0 length([behav, lfov])+1 0 4])
% rr_null = rr(6:end_nonlfov); % first 5 are the AQ animals that have a lot of reward experience so have very high reward rates
% rr_null = rr_null(~isnan(rr_null));

rr_null = rr(behav); % first 5 are the AQ animals that have a lot of reward experience so have very high reward rates
rr_null = rr_null(~isnan(rr_null));

rr_lfov = rr(lfov); 

xr = gb_rand_jitter(rr_null, 2);
xl = gb_rand_jitter(rr_lfov, 2);
subplot(1,3,3)
hold on
bar(1, median(rr_null), 'FaceColor', [.8 .8 .8])
c = 'k';
scatter(1+xr, rr_null, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .75)
c = 'b';
bar(2, median(rr_lfov), 'FaceColor', 'b', 'FaceAlpha', .25)
scatter(2+xl, rr_lfov, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .5)
a = ceil(1.3*max([rr_null; rr_lfov]));
axis([0 3 0 a])
ylabel('Rewards/min')
title(sprintf('Session %d', sess2use))
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Control', 'MiniLFOV'}, 'XTickLabelRotation', 0)

[p, h] = ranksum(rr_null, rr_lfov);
a = 1.25*max([rr_null; rr_lfov]);
plot([1, 2], [a, a], 'k', 'LineWidth', 2)
plot([1, 1], [a, a*.95], 'k', 'LineWidth', 2)
plot([2, 2], [a, a*.95], 'k', 'LineWidth', 2)
if p<=.05
text(1.5, a, sprintf('* p=%1.2e', p), 'FontSize', 10)
else
text(1.5*.95, 1.05*a, 'ns', 'FontSize', 10)
end


