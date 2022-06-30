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
% save('C:\Users\gjb326\Desktop\Sample Data\Linear Maze Data\AllData.mat', 'A');
%%
% controls = {'A80','A81','A82','A83','A84','A85','A86','A87','A88','A89','A90','A91','FC1','FC2','FC5','FC7','FC9','FC10','FC11'};
controls = {'A80','A81','A82','A83','A84','A85','A86','A87','A88','A89','A90','A91','FC2','FC5','FC9','FC10','FC11'}; % removed the animals that never reached criterion
ncontrols = length(controls);
implants ={'Hipp6'};
for i = 1:length(hippnums)
    controls{end+1} = sprintf('Hipp%d', hippnums(i));
end
rewardrate = NaN(length(controls),length(sess2use));
shortpref = NaN(length(controls),length(sess2use));
% allan = 
% nsess = (length(controls) + length(implants))*sess2use;
iscontrol = []
for i = 1:length(controls)
    for j = 1:length(sess2use)
        
        fieldss = sprintf('A.animal%s.sess%d;', controls{i}, sess2use(j));
        try
            temp = eval(fieldss);
            
%             ind = all_an(i);
            a = temp.numr;
            b = (temp.t(end)/10^6)/60;
            c = 'k';
            lr = temp.lr;
            if temp.isLFOV
                c = 'b';
            end
            
            %     scatter(i, temp.lr, 15, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .5)
            if b<=21 && b>=7 && temp.lr>0 % temp.dual == true
                if i < ncontrols
                    iscontrol(i) = true;
                else
                    iscontrol(i)=false;
                end
                rewardrate(i,j) = a/b; % rewards per min
                shortpref(i,j) = temp.sr/temp.lr; % rewards per min
                %     subplot_tight(1, 3, 1:2, [.2,.05])
                subplot(1, 3, 1:2)
                hold on
                %     scatter(i, lr, 15, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .5)
                scatter3(i, j, rewardrate(i,j), 15, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .5)
                %     h = text(ind, rr(ind)*1.2, temp.aname, 'Interpreter', 'latex');
            end
        catch
            fprintf('\n%s', fieldss)
        end
        
    end
    
end
%
rateCrit = NaN(length(controls),1);
choiceCrit = NaN(length(controls),1);
bothCrit = NaN(length(controls),1);
for i = 1:size(rewardrate,1)
    ind = find(rewardrate(i,:)>2 & shortpref(i,:)>2);
    try
        bothCrit(i) = sess2use(ind(1));
    end
    try
        rateCrit(i) = sess2use(find(rewardrate(i,:)>2, 1, 'first'));
    end
    try
        choiceCrit(i) = sess2use(find(shortpref(i,:)>2, 1, 'first'));
    end
end
goodinds = ~isnan(bothCrit);
iscontrol = iscontrol(goodinds);
bothCrit = bothCrit(goodinds);
rateCrit = rateCrit(goodinds);
choiceCrit = choiceCrit(goodinds);
% rateCrit(isnan(rateCrit)) = 12;
% choiceCrit(isnan(choiceCrit)) = 12;

null_both = bothCrit(iscontrol==true);
lfov_both = bothCrit(iscontrol==false);
[p, h] = ranksum(null_both, lfov_both)

null_rate = rateCrit(iscontrol==true);
lfov_rate = rateCrit(iscontrol==false);
[p, h] = ranksum(null_rate, lfov_rate)

null_choice = choiceCrit(iscontrol==true);
lfov_choice = choiceCrit(iscontrol==false);
[p, h] = ranksum(null_choice, lfov_choice)

%%
critSpd = NaN(30,1);
critRR = NaN(30,1);
day1Spd = NaN(30,1);
day1RR = NaN(30,1);
for i = 1:size(bothCrit,1)
    %%
    critSess = bothCrit(i);
    critDay = sprintf('A.animal%s.sess%d;', controls{i}, critSess);
    c = eval(critDay);
    if i < 18
        badinds = c.x==0 & c.y==0;
        c.x(find(badinds)) = interp1(find(~badinds), c.x(~badinds), find(badinds), 'linear');
        c.y(find(badinds)) = interp1(find(~badinds), c.y(~badinds), find(badinds), 'linear');
        
        xo = c.x; yo = c.y;
        xcenter = 350; % center x,y point on the short path
        ycenter = 150; % center x,y point on the short path
        ymin = 190-ycenter;
        yextreme = 346-ycenter;
        xmin = 120-xcenter;
        xextreme = 220-xcenter;
        xo = xo-xcenter;
        yo = yo-ycenter;
        xs = (xextreme - xmin)/(xextreme-xcenter);
        ys = (yo-ymin)/(yextreme-ymin);
        xo = xo +  xo*xs.*ys;
        
        % Shorten y and x dims
        x = 250*xo/440;
        y = 125*yo/(149-10);
        cwin = gausswin(11); cwin = cwin./sum(cwin);
        x = conv(x, cwin, 'same');
        y = conv(y, cwin, 'same');
        
        
        dt = [median(diff(c.t')); diff([c.t'])]./1000000;
        spd = sqrt(diff([x(1); x']).^2 + diff([y(1); y']).^2);
        spd = spd./dt;
    else
        x = c.x;
        y = c.y;
        dt = [median(diff(c.t)); diff([c.t])]./1000000;
        median(dt)
        spd = sqrt(diff([x(1); x]).^2 + diff([y(1); y]).^2)./dt;
    end
    center = x<=100 & x>=-100 & y <= 50;
    spd(~center) = NaN;
    
    figure(1); %clf;
    subplot(6,6,i)
    histogram(spd, [0:5:250])
%     subplot(1,3,2)
%     plot(y, 'r'); hold on; plot(c.y, 'k')
%     subplot(1,3,3)
%     hold on
% %     subplot(1,2,1)
%     plot3(x, y, spd)
%     plot3(c.x, c.y, spd)
    critSpd(i) = nanmedian(spd);
    critRR(i) = rewardrate(i, critSess);
end
for i = 1:size(bothCrit,1)
    %%
    try
    critSess = 1;%bothCrit(i);
    critDay = sprintf('A.animal%s.sess%d;', controls{i}, critSess);
    c = eval(critDay);
    catch
    critSess = 2;%bothCrit(i);
    critDay = sprintf('A.animal%s.sess%d;', controls{i}, critSess)
    c = eval(critDay);
    end
    if i < 18
        badinds = c.x==0 & c.y==0;
        c.x(find(badinds)) = interp1(find(~badinds), c.x(~badinds), find(badinds), 'linear');
        c.y(find(badinds)) = interp1(find(~badinds), c.y(~badinds), find(badinds), 'linear');
        
        xo = c.x; yo = c.y;
        xcenter = 350; % center x,y point on the short path
        ycenter = 150; % center x,y point on the short path
        ymin = 190-ycenter;
        yextreme = 346-ycenter;
        xmin = 120-xcenter;
        xextreme = 220-xcenter;
        xo = xo-xcenter;
        yo = yo-ycenter;
        xs = (xextreme - xmin)/(xextreme-xcenter);
        ys = (yo-ymin)/(yextreme-ymin);
        xo = xo +  xo*xs.*ys;
        
        % Shorten y and x dims
        x = 250*xo/440;
        y = 125*yo/(149-10);
        cwin = gausswin(11); cwin = cwin./sum(cwin);
        x = conv(x, cwin, 'same');
        y = conv(y, cwin, 'same');
        
%         center = x<=120 & x>=-120;
        dt = [median(diff(c.t')); diff([c.t'])]./1000000;
        median(dt)
        spd = sqrt(diff([x(1); x']).^2 + diff([y(1); y']).^2);
        spd = spd./dt;
    else
        x = c.x;
        y = c.y;
        dt = [median(diff(c.t)); diff([c.t])]./1000000;
        median(dt)
        spd = sqrt(diff([x(1); x]).^2 + diff([y(1); y]).^2)./dt;
    end
    center = x<=100 & x>=-100 & y <= 50;
    spd(~center) = NaN;

    day1Spd(i) = nanmedian(spd);
    day1RR(i) = rewardrate(i, critSess);
end
day1Spd_C = day1Spd(iscontrol==true);
day1Spd_L = day1Spd(iscontrol==false);
[ps1, h] = ranksum(day1Spd_C, day1Spd_L)

critSpd_C = critSpd(iscontrol==true);
critSpd_L = critSpd(iscontrol==false);
[ps2, h] = ranksum(critSpd_C, critSpd_L)

day1RR_C = day1RR(iscontrol==true);
day1RR_L = day1RR(iscontrol==false);
[pr1, h] = ranksum(day1RR_C, day1RR_L)

critRR_C = critRR(iscontrol==true);
critRR_L = critRR(iscontrol==false);
[pr2, h] = ranksum(critRR_C, critRR_L)

%%
figure(3); clf; 
subplot(1,2,1);
hold on
sub_barscatter(day1Spd_C, .5, [.9 .9 .9])
sub_barscatter(day1Spd_L, 1.5, [.2 .5 .8])

sub_barscatter(critSpd_C, 3.5, [.7 .7 .7])
sub_barscatter(critSpd_L, 4.5, [.2 .5 .8]*.75)
ylabel('Running speed (cm/sec)')
xlabel(sprintf('\nSession 1            Criterion Sess'))
set(gca, 'XTick', [.5 1.5 3.5, 4.5], 'XTickLabel', {'Con' 'LFOV' 'Con' 'LFOV'})
ylim([0 200])
sub_sigbar([.5 1.5], [155 170], ps1)
sub_sigbar([3.5 4.5], [155 170], ps2)


subplot(1,2,2);
hold on
sub_barscatter(day1RR_C, .5, [.9 .9 .9])
sub_barscatter(day1RR_L, 1.5, [.2 .5 .8])

sub_barscatter(critRR_C, 3.5, [.7 .7 .7])
sub_barscatter(critRR_L, 4.5, [.2 .5 .8]*.75)
ylabel('Rewards/min')
xlabel(sprintf('\nSession 1            Criterion Sess'))
set(gca, 'XTick', [.5 1.5 3.5, 4.5], 'XTickLabel', {'Con' 'LFOV' 'Con' 'LFOV'})
ylim([0 8])
sub_sigbar([.5 1.5], [5 5.55], pr1)
sub_sigbar([3.5 4.5], [5 5.55], pr2)

%%
figure(99); clf; hold on;
% c = 'k';
% scatter(null_rate, null_choice, 15, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .75)
% c = 'b';
% scatter(lfov_rate, lfov_choice, 15, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .5)
% axis([0 11 0 11])

X{1} = null_rate;
X{2} = null_choice;
X{3} = null_both;
Z{3} = lfov_both;
Z{2} = lfov_choice;
Z{1} = lfov_rate;
for i = 1:3
    
    rr_null = X{i};
    rr_lfov = Z{i};
    xr = gb_rand_jitter(rr_null, 2);
    xl = gb_rand_jitter(rr_lfov, 2);
    subplot(1,3,i)
    hold on
    bar(1, median(rr_null), 'FaceColor', [.8 .8 .8])
    c = 'k';
    scatter(1+xr, rr_null, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .75)
    c = 'b';
    bar(2, median(rr_lfov), 'FaceColor', 'b', 'FaceAlpha', .25)
    scatter(2+xl, rr_lfov, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerFaceAlpha', .5)
    a = 15;% ceil(1.3*max([rr_null; rr_lfov]));
    axis([0 3 0 a])
    ylabel('Rewards/min')
    % title(sprintf('Session %d', sess2use))
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'Control', 'MiniLFOV'}, 'XTickLabelRotation', 0)
    
    [p, h] = ranksum(rr_null, rr_lfov);
    a = .90*a;%max([rr_null; rr_lfov]);
    plot([1, 2], [a, a], 'k', 'LineWidth', 2)
    plot([1, 1], [a, a*.95], 'k', 'LineWidth', 2)
    plot([2, 2], [a, a*.95], 'k', 'LineWidth', 2)
    if p<=.05
        text(1.5*.65, 1.08*a, sprintf('* p=%1.2e', p), 'FontSize', 10)
    else
        text(1.5*.95, 1.08*a, 'ns', 'FontSize', 10)
    end
    switch i
        case 1
            ylabel('Sessions to meet rate criterion')
        case 2
            ylabel('Sessions to meet choice criterion')
        case 3
            ylabel('Sessions to meet both criterion')
    end
    
end


function sub_barscatter(ydata, xpos, colors)
    xr = gb_rand_jitter(ydata, 2);
    c = colors;
    bar(xpos, median(ydata), 'FaceColor', c.*.75)
    scatter(xpos+xr, ydata, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .75)
end
function sub_sigbar(xpos, ypos, p)
%     a = .90*a;%max([rr_null; rr_lfov]);
    plot(xpos, [max(ypos) max(ypos)], 'k', 'LineWidth', 2)
    plot([xpos(1), xpos(1)], [ypos], 'k', 'LineWidth', 2)
    plot([xpos(2), xpos(2)], [ypos], 'k', 'LineWidth', 2)
    if p<=.001
        text(mean(xpos) - abs(diff(xpos)), max(ypos) + abs(diff(ypos)), sprintf('***p=%1.2e', p), 'FontSize', 10)
    elseif p<=.01
        text(mean(xpos) - abs(diff(xpos)), max(ypos) + abs(diff(ypos)), sprintf('**p=%1.2e', p), 'FontSize', 10)
    elseif p<=.05
        text(mean(xpos) - abs(diff(xpos)), max(ypos) + abs(diff(ypos)), sprintf('*p=%1.2e', p), 'FontSize', 10)
    else
        text(mean(xpos) - abs(diff(xpos)), max(ypos) + abs(diff(ypos)), sprintf('ns, p=%1.2e', p), 'FontSize', 10)
    end
end
