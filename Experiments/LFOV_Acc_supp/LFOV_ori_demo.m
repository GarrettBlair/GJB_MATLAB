% ori_ts = ORI_Data.TimeStamp_ms_(shared_ts);
% w = ORI_Data.qw(shared_ts);
% x = ORI_Data.qx(shared_ts);
% y = ORI_Data.qy(shared_ts);
% z = ORI_Data.qz(shared_ts);
% [roll2, pitch2, yaw2] = quaternion_conversion_LFOV(w, x ,y, z);


% FIGDIR = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LFOV_ori_supp'

ddir = 'D:\Sample Data\Linear Maze Data\Hipp_processed\';
cd(ddir)
fnames = dir('*linear*');

% a = [6 7 8 9 12 13 15 18 35 37];
% aname = {'6' '7' '8' '9' '12' '13' '15' '18' '35' '37'};
a = [6 8 9];
aname = {'6' '8' '9'};
% s = [3:5];
s = [1:10];

        n = 24;
hmap(1:n,1) = linspace(0,1,n);
hmap(:,[2 3]) = 1; %brightness
huemap = hsv2rgb(hmap);
% colormap(huemap)
        mcolors = huemap;%jet(n);
% xbins = [-200:100:200];
% ybins = [-50 50 100 150];
xbins = [-200 -100:20:100, 200];
ybins = [-50 -20:20:120 160];%50 100 150];
yawbins = [-pi:2*pi/n:pi];
% binstep = 2*pi/n
% yawbins = [-pi-binstep/2:binstep:pi+binstep/2];
% polarmap = cell(length(xbins), length(ybins));
polarmap = zeros(length(xbins)-1, length(ybins)-1, length(yawbins)-1);
cornermap = [];
cornermap.y = zeros(length(yawbins)-1, 1);
cornermap.r = zeros(length(yawbins)-1, 1);
cornermap.p = zeros(length(yawbins)-1, 1);
cornermap.std = [];
cornermap.mean = [];
pathmap = [];
pathmap.y = zeros(length(yawbins)-1, 1);
pathmap.r = zeros(length(yawbins)-1, 1);
pathmap.p = zeros(length(yawbins)-1, 1);
beemap = [];
beemap.y_sl = zeros(length(yawbins)-1, 1);
beemap.y_sr = zeros(length(yawbins)-1, 1);
beemap.r_sl = zeros(length(yawbins)-1, 1);
beemap.r_sr = zeros(length(yawbins)-1, 1);
beemap.p_sl = zeros(length(yawbins)-1, 1);
beemap.p_sr = zeros(length(yawbins)-1, 1);

x_q = -100:10:100; % interp values for evaluating ori sensor consistency
min_bee_num = 8;
%
sr_bee_mean = [];
sr_bee_std = [];
sr_bee_animal = [];
sl_bee_mean = [];
sl_bee_std = [];
sl_bee_animal = [];

sr_bee_indx = 0;
sl_bee_indx = 0;

plotpaths = false;
for aLoop = 1:length(a)
    for sLoop = 1:length(s)
        %%
        fname = sprintf('%sHipp%d_linear%d.mat', ddir, a(aLoop), s(sLoop));
        if isfile(fname)
            temp = load(fname);
            x = temp.ms.x; y = temp.ms.y; t = temp.ms.time;
            q = temp.ms.ori; %[ORI_Data.qw, ORI_Data.qx, ORI_Data.qy, ORI_Data.qz];
            shared_ts = ismember(q(:,1), temp.ms.time);
            q = q(shared_ts, :);
            %         [~, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
            [roll, pitch, yaw] = quaternion_conversion_LFOV(q(:,2), q(:,3), q(:,4), q(:,5));
            pitch = real(pitch);
            %        if std(yaw) >=.03
            %%
            if plotpaths==true
                figure(a(aLoop)*100 + s(sLoop)); clf
                hold on
                plot3(temp.ms.x, temp.ms.y, temp.ms.time./1000,'k');
            end
            %         plot3(xb, yb, temp.ms.time./1000,'k');
            yawb = discretize(yaw, yawbins);
            [xycounts, ~, ~, xb, yb] = histcounts2(x, y, xbins, ybins);
            %         yb = discretize(y, [-150:25:150]);
            if plotpaths==true
                
                for i = 1:max(yawb)
                    ind = yawb==i;
                    if any(ind)
                        scatter3(temp.ms.x(ind), temp.ms.y(ind), temp.ms.time(ind)./1000, 20,'o',...
                            'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8);
                        %                 scatter3(xb(ind), yb(ind), temp.ms.time(ind)./1000, 20,'o',...
                        %                     'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .7);
                    end
                end
                set(gca, 'View', [-143.4612   21.8522])
                drawnow
                
            end
            long_end1 = x>100 & y>100;
            long_end2 = x<-100 & y>100;
            short_end1 = x>100 & y<25;
            short_end2 = x<-100 & y<25;
            northsouth= x<=100 | x>=-100;
            eastwest= ~northsouth;
            in_corner = long_end1 | long_end2 | short_end1 | short_end2;
            %         in_NS_path = (~in_corner & x<=100 ) | (~in_corner & x>=-100)  ;
            %         in_EW_path = ~in_corner & eastwest== 1;
            cornermap.y = cornermap.y + histcounts(yaw(in_corner), yawbins)';
            cornermap.r = cornermap.r + histcounts(roll(in_corner), yawbins)';
            cornermap.p = cornermap.p + histcounts(pitch(in_corner), yawbins)';
            pathmap.y = pathmap.y + histcounts(yaw(~in_corner), yawbins)';
            pathmap.r = pathmap.r + histcounts(roll(~in_corner), yawbins)';
            pathmap.p = pathmap.p + histcounts(pitch(~in_corner), yawbins)';
            
            %         bee_ind = temp.ms.short_bees;
            sl_bee_ind = temp.ms.SL.bee;
            sr_bee_ind = temp.ms.SR.bee;
            beemap.y_sl = beemap.y_sl + histcounts(yaw(sl_bee_ind), yawbins)';
            beemap.r_sl = beemap.r_sl + histcounts(roll(sl_bee_ind), yawbins)';
            beemap.p_sl = beemap.p_sl + histcounts(pitch(sl_bee_ind), yawbins)';
            beemap.y_sr = beemap.y_sr + histcounts(yaw(sr_bee_ind), yawbins)';
            beemap.r_sr = beemap.r_sr + histcounts(roll(sr_bee_ind), yawbins)';
            beemap.p_sr = beemap.p_sr + histcounts(pitch(sr_bee_ind), yawbins)';
            
            
            sr_bee_ind = temp.ms.SR.is_bee;
            %         figure; hold on
            sl_bee_yaw = [];
            for i = 1:size(temp.ms.SL.pass_ind,1)
                if temp.ms.SL.is_bee(i)
                    xx= x(temp.ms.SL.pass_ind(i,1):temp.ms.SL.pass_ind(i,2));
                    yy= yaw(temp.ms.SL.pass_ind(i,1):temp.ms.SL.pass_ind(i,2));
                    y_q = interp1(xx,yy, x_q);
                    %                 plot(xx,yy, 'k')
                    %                 plot(x_q, y_q, 'r')
                    sl_bee_yaw = cat(1, sl_bee_yaw, y_q);
                    
                end
            end
            %         figure; hold on
            sr_bee_yaw = [];
            for i = 1:size(temp.ms.SR.pass_ind,1)
                if temp.ms.SR.is_bee(i)
                    xx= x(temp.ms.SR.pass_ind(i,1):temp.ms.SR.pass_ind(i,2));
                    yy= yaw(temp.ms.SR.pass_ind(i,1):temp.ms.SR.pass_ind(i,2));
                    y_q = interp1(xx,yy, x_q);
                    %                 plot(xx,yy, 'b')
                    %                 plot(x_q, y_q, 'm')
                    sr_bee_yaw = cat(1, sr_bee_yaw, y_q);
                    
                end
            end
            nsr = size(sr_bee_yaw,1);
            nsl = size(sl_bee_yaw,1);
            if nsr>min_bee_num
                sr_bee_indx = sr_bee_indx+1;
                sr_bee_mean = cat(1, sr_bee_mean, mean(sr_bee_yaw,1));
                %         sr_bee_std(bee_indx, :) = std(sr_bee_yaw,[],1)./(sqrt(nsr-1));
                sr_bee_std = cat(1, sr_bee_std, std(sr_bee_yaw,[],1));
                sr_bee_animal = cat(1, sr_bee_animal, a(aLoop));
            end
            if nsl>min_bee_num
                sl_bee_indx = sl_bee_indx+1;
                sl_bee_mean = cat(1, sl_bee_mean, mean(sl_bee_yaw,1));
                %         sr_bee_std(bee_indx, :) = std(sr_bee_yaw,[],1)./(sqrt(nsr-1));
                sl_bee_std = cat(1, sl_bee_std, std(sl_bee_yaw,[],1));
                sl_bee_animal = cat(1, sl_bee_animal, a(aLoop));
            end
            if nsl>min_bee_num  || nsr>min_bee_num
                cornermap.std = cat(1, cornermap.std, nanstd(yaw(in_corner)));
                cornermap.mean = cat(1, cornermap.mean, nanmean(yaw(in_corner)));
            end
            
            %         NS_pathmap = NS_pathmap + histcounts(yaw(in_NS_path), yawbins)';
            %         EW_pathmap = EW_pathmap + histcounts(yaw(in_EW_path), yawbins)';
            for xl = 1:length(xbins)-1
                for yl = 1:length(ybins)-1
                    ind = xb==xl & yb==yl;
                    if any(ind)
                        
                        polarmap(xl, yl, :) = squeeze(polarmap(xl, yl, :)) + histcounts(yaw(ind), yawbins)';
                    end
                    
                end
            end
            
            %        else
            %           fprintf('\nBad %s', fname)
            %        end
        else
            
            fprintf('%s not found', fname)
        end
        
        
    end
end

%%
% figure(23); clf;
xb = yawbins(2:end)+mean(abs(diff(yawbins)))/2; 
yaw_xb = circshift(xb, -4); % circshift the bins to match changs convention
yaw_xb = [yaw_xb, yaw_xb(1)]; % just to make it circular

roll_xb = circshift(xb, -10); % circshift the bins to match changs convention
roll_xb = [roll_xb, roll_xb(1)]; % just to make it circular

pitch_xb = circshift(xb, -1); % circshift the bins to match changs convention
pitch_xb = [pitch_xb, pitch_xb(1)]; % just to make it circular

% subplot(3,2,1)
% polarplot(yaw_xb, [cornermap.y; cornermap.y(1)]./max(cornermap.y))
% set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
% subplot(3,2,2)
% polarplot(yaw_xb, [pathmap.y; pathmap.y(1)]./max(pathmap.y))
% set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
% subplot(3,2,3)
% polarplot(roll_xb, [cornermap.r; cornermap.r(1)]./max(cornermap.r))
% set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
% subplot(3,2,4)
% polarplot(roll_xb, [pathmap.r; pathmap.r(1)]./max(pathmap.r))
% set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
% subplot(3,2,5)
% polarplot(pitch_xb, [cornermap.p; cornermap.p(1)]./max(cornermap.p))
% set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
% subplot(3,2,6)
% polarplot(pitch_xb, [pathmap.p; pathmap.p(1)]./max(pathmap.p))
% set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')

figure(24); clf
subplot(3,3,1); v = cornermap.y;
polarplot(yaw_xb, [v; v(1)]./max(v))
title('yaw, corner')
set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
subplot(3,3,2); v = beemap.y_sl;
polarplot(yaw_xb, [v; v(1)]./max(v))
set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
subplot(3,3,3); v = beemap.y_sr;
polarplot(yaw_xb, [v; v(1)]./max(v))
set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')

subplot(3,3,4); v = cornermap.r;
polarplot(roll_xb, [v; v(1)]./max(v))
title('roll, corner')
set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
subplot(3,3,5); v = beemap.r_sl;
polarplot(roll_xb, [v; v(1)]./max(v))
set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
subplot(3,3,6); v = beemap.r_sr;
polarplot(roll_xb, [v; v(1)]./max(v))
set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')

subplot(3,3,7); v = cornermap.p;
polarplot(pitch_xb, [v; v(1)]./max(v))
title('pitch, corner')
set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
subplot(3,3,8); v = beemap.p_sl;
polarplot(pitch_xb, [v; v(1)]./max(v))
set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
subplot(3,3,9); v = beemap.p_sr;
polarplot(pitch_xb, [v; v(1)]./max(v))
set(gca, 'RTick', [.2 .4 .6 .8 1], 'RTickLabel', '', 'ThetaTick', [0:30:360], 'ThetaTickLabel', '')
%%
figure(35); clf; 
set(gcf, 'Position', [192   404   954   381])
subplot(1,2,1)
hold on
cm = jet(length(a)*3).*.7;
cm1 = light_colormap(cm(1:3:9,:), 2);
cm2 = light_colormap(cm(1:3:9,:), -2);

cornerx = linspace(-180, -110, length(cornermap.mean));
marker = {'.' 'x' '*'};
for i = 1:length(cornermap.mean)
 scatter(cornerx(i), cornermap.mean(i), 'Marker', marker{sr_bee_animal(i)==a}, 'MarkerEdgeColor', cm1(sr_bee_animal(i)==a, :))
 plot([cornerx(i) cornerx(i)], [cornermap.mean(i)-cornermap.std cornermap.mean(i)+cornermap.std],...
     'Color', cm2(sr_bee_animal(i)==a,:)/2)
end
for i = 1:size(sr_bee_mean,1)
    plot(x_q, sr_bee_mean(i,:), '-', 'Color', cm1(sr_bee_animal(i)==a,:))
end
for i = 1:size(sl_bee_mean,1)
    plot(x_q, sl_bee_mean(i,:), '-', 'Color', cm2(sr_bee_animal(i)==a,:))
end
set(gca, 'XTick', [-100:50:100], 'YTick', [-pi:pi/2:pi], 'YTickLabel', [270 0 90 180 270])
axis([-200 120 -3.5 3.5])
xlabel('xpos')
ylabel('Mean Yaw')

subplot(1,2,2)
hold on
shadedErrorBar([-150 150], [mean(cornermap.std) mean(cornermap.std)], [std(cornermap.std) std(cornermap.std)])
for i = 1:size(sr_bee_mean,1)
    plot(x_q, sr_bee_std(i,:), '-', 'Color', cm1(sr_bee_animal(i)==a,:))
end
for i = 1:size(sl_bee_mean,1)
    plot(x_q, sl_bee_std(i,:), '-', 'Color', cm2(sr_bee_animal(i)==a,:))
end
axis([-120 120 -.1 2.2])
xlabel('xpos')
ylabel('Yaw Standard dev')
%%
% subplot(3,2,6)
% polarplot(xb, [pathmap.p; pathmap.p(1)]./max(pathmap.p))

% subplot(133)
% polarplot(xb, [EW_pathmap; EW_pathmap(1)])
%%
aLoop = 4;
sLoop = 6;
fname = sprintf('%sHipp%d_linear%d.mat', ddir, a(aLoop), s(sLoop));
        figure(99); clf
subplot(121);
colorwheel;
set(gca, 'YDir', 'normal', 'XDir', 'reverse')
if isfile(fname)
    temp = load(fname);
    sub = 4351:4864;
%     sub = 1:2000;
    x = temp.ms.x(sub); y = temp.ms.y(sub);
    q = temp.ms.ori; %[ORI_Data.qw, ORI_Data.qx, ORI_Data.qy, ORI_Data.qz];
    shared_ts = ismember(q(:,1), temp.ms.time);
    q = q(shared_ts, :);
    %         [~, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
    [roll, pitch, yaw] = quaternion_conversion_LFOV(q(:,2), q(:,3), q(:,4), q(:,5));
    pitch = real(pitch);
    roll = roll(sub); yaw = yaw(sub); pich = pitch(sub);
        %%
        figure(a(aLoop));
        subplot(122); cla
        hold on
%         plot3(temp.ms.x, temp.ms.y, temp.ms.time./1000,'k');
        ticks = 1:length(x);
        plot(x, y,'k');
%         plot3(x, y, ticks,'k');
        %         plot3(xb, yb, temp.ms.time./1000,'k');
        yawb = discretize(yaw, yawbins);
        [xycounts, ~, ~, xb, yb] = histcounts2(x, y, xbins, ybins);
        %         yb = discretize(y, [-150:25:150]);
        for i = 1:max(yawb)
            ind = yawb==i;
            if any(ind)
%                 scatter3(x(ind), y(ind), temp.ms.time(ind)./1000, 20,'o',...
%                     'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8);
%                 scatter3(x(ind), y(ind), ticks(ind), 20,'o',...
%                     'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8);
                scatter(x(ind), y(ind), 20,'o',...
                    'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8);
                %                 scatter3(xb(ind), yb(ind), temp.ms.time(ind)./1000, 20,'o',...
                %                     'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .7);
            end
        end  
        scatter(x(1), y(1), 30, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
        scatter(x(end), y(end), 30, 'd', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
        axis([-200 200 -50 200])
%         axis equal
        set(gca, 'View', [-90 90])
        xlabel('X Position (cm)')
        ylabel('Y Position (cm)')
        drawnow
else
    
    fprintf('%s not found', fname)
end
%%
aLoop = 4;
sLoop = 6;
fname = sprintf('%sHipp%d_linear%d.mat', ddir, a(aLoop), s(sLoop));
        figure(99); clf
% subplot(121);
% colorwheel;
% set(gca, 'YDir', 'normal', 'XDir', 'reverse')
if isfile(fname)
%     temp = load(fname);
%     sub = 4351:4864;
    sub = 1:2800;
    x = temp.ms.x(sub); y = temp.ms.y(sub); t = temp.ms.time(sub)./1000;
    q = temp.ms.ori; %[ORI_Data.qw, ORI_Data.qx, ORI_Data.qy, ORI_Data.qz];
    shared_ts = ismember(q(:,1), temp.ms.time);
    q = q(shared_ts, :);
    %         [~, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
    [roll, pitch, yaw] = quaternion_conversion_LFOV(q(:,2), q(:,3), q(:,4), q(:,5));
    pitch = real(pitch);
    roll = roll(sub); yaw = yaw(sub); pich = pitch(sub);
        %
        figure(a(aLoop)); clf
%         subplot(122); cla
        hold on
        ticks = 1:length(x);
%         plot(x, y,'k');
        plot3(x, y, t,'k');
        %         plot3(xb, yb, temp.ms.time./1000,'k');
        yawb = discretize(yaw, yawbins);
        [xycounts, ~, ~, xb, yb] = histcounts2(x, y, xbins, ybins);
        %         yb = discretize(y, [-150:25:150]);
        for i = 1:max(yawb)
            ind = yawb==i;
            if any(ind)
%                 scatter3(x(ind), y(ind), temp.ms.time(ind)./1000, 20,'o',...
%                     'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8);
                scatter3(x(ind), y(ind), t(ind), 20,'o',...
                    'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8);
%                 scatter(x(ind), y(ind), 20,'o',...
%                     'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8);
                %                 scatter3(xb(ind), yb(ind), temp.ms.time(ind)./1000, 20,'o',...
                %                     'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .7);
            end
        end  
        scatter(x(1), y(1), 30, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
        scatter(x(end), y(end), 30, 'd', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
        axis([-200 200 -50 200])
%         axis equal
        set(gca, 'View', [-90 90])
        xlabel('X Position (cm)')
        ylabel('Y Position (cm)')
        zlabel('Time (sec')
        drawnow
else
    
    fprintf('%s not found', fname)
end



%%
q = ms.ori; %[ORI_Data.qw, ORI_Data.qx, ORI_Data.qy, ORI_Data.qz];
% q = q(shared_ts, :);
[~, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
yaw = yaw(1:2:end);
% ms.ori = ORI_Data(shared_ts, :);
% % [ori_str] = LFOV_ori_extract(q, [0:30:360]);
% 
% ms.ori.time = ori_ts; 
% ms.ori.roll = roll; 
% ms.ori.pitch = pitch; 
% ms.ori.yaw = yaw;

figure(9); clf
hold on
plot3(ms.x, ms.y, ms.time./1000,'K');
n =16;
yb = discretize(yaw, [-pi/2:pi/n:pi/2]);
mcolors = jet(n);
for i = 1:max(yb)
ind = yb==i;
if any(ind)
scatter3(ms.x(ind), ms.y(ind), ms.time(ind)./1000, 20,'o',...
    'MarkerFaceColor', mcolors(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .7);
end
end

%%
figure(1); clf; plot3(x,y,t);


figure(2); 
clf; hold on
plot(yaw, 'k');
plot(dx, 'r');

    figure(3); clf
    hold on
    plot3([x], [y], [t], 'k.-')
%     plot([x(i), x(i+1)], [y(i), y(i+1)], 'k.-')
    xd = cos(yaw+yaw_offset);
    yd = sin(yaw+yaw_offset);
    xd2 = cos(dx);
    yd2 = sin(dx);
    for i = 4:length(x)-1
    for j = 1%2:length(xd)-1
    plot3([x(i+j-1) x(i+j-1)+xd(j)], [y(i+j-1) y(i+j-1)+yd(j)], [t(i+j-1) t(i+j-1)], 'r-')
%     plot3([x(i+j-1) x(i+j-1)+xd2(j)], [y(i+j-1) y(i+j-1)+yd2(j)], [t(i+j-1) t(i+j-1)], 'm-')
    end
    end
    axis([-40 40 -40 40 0 100])
    drawnow

%%
win = 1;

yawms = ms.ori.yaw(m_inds);
yaw = interp1(ms.timestamps(m_inds), yawms, bt);
yaw_offset = pi;%-pi/2;'

dx = atan(diff(x)./diff(y));


imf = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\Hipp16942\2022_06_28_15_58_07\behavCam.tiff';
for i = 1:1:500%length(x)-100
    %%
    figure(1); clf
%     bi = find(min(abs(bt./1000-t(i))) == abs(bt./1000-t(i)))    +  b_inds(1);
    im = imread(imf, b_inds(i));
    imagesc(im)
    hold on
%     plot3([x(i:i+win)], [y(i:i+win)], [t(i:i+win)], 'k.-')
    plot([x(i:i+win)], [y(i:i+win)], 'k.-')
%     plot([x(i), x(i+1)], [y(i), y(i+1)], 'k.-')
    xd = 4*cos(yaw(i:i+win) + yaw_offset);
    yd = 4*sin(yaw(i:i+win) + yaw_offset);
    xd2 = 4*cos(dx(i:i+win));
    yd2 = 4*sin(dx(i:i+win));
    for j = 1:length(xd)
%     plot3([x(i+j-1) x(i+j-1)+xd(j)], [y(i+j-1) y(i+j-1)+yd(j)], [t(i+j-1) t(i+j-1)], 'r-')
    plot([x(i+j-1) x(i+j-1)+xd(j)], [y(i+j-1) y(i+j-1)+yd(j)], 'r-')
%     plot3([x(i+j-1) x(i+j-1)+xd2(j)], [y(i+j-1) y(i+j-1)+yd2(j)], [t(i+j-1) t(i+j-1)], 'm-')
    end
% %     axis([-40 40 -40 40 0 1000])
    axis([130 230 200 300])
    drawnow
    pause(.05)
%     plot([x(i), x(i)+xd], [y(i), y(i)+yd], 'r-')
end
