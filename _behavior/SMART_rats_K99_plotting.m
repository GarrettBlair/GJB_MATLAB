clear
outdir = "E:\SuaKim\SMART task data\behavior output\";
% anames = {  'SUA34991', 'SUA34992', 'SUA34993', 'SUA34994', 'SUA34995', 'SUA34996',...
%     'SUA35008', 'SUA35009', 'SUA35010', 'SUA35011', 'SUA35012', 'SUA35013'};
metricsFileName = sprintf('%s%s', outdir, '\SUA_behavior_metrics.mat');
load(metricsFileName)
% SMART_rats_sua_plotting2.m
%%

metersran_per_min = ((distRan./(sessmins*60))*100);
% vars2plot = {'entrpermin', 'metersran_per_min', 'right_preference_room', 'probability_correct', 'average_streak'}; 
% names2plot = {'rewards/min', 'av speed (cm/sec)', 'reinf. zone pref.', 'prob. correct', 'av. streak'}; 
% ylims = {[0, 12], [0, 50], [.2, 1.0], [.2, 1.0], [0, 10]};

% vars2plot = {'metersran_per_min', 'probability_correct', 'swap_performance_win_win', 'right_preference_room'}; 
% names2plot = {'av speed (cm/sec)', 'prob. correct', 'w-w', 'rpref'}; 
% ylims = {[-10, 10], [-.4, .4], [-.4, .4], [-.4, .4]};
strat = (swap_performance_win_win);% - swap_performance_lose_win) ./(swap_performance_win_win + swap_performance_lose_win) 
vars2plot = {'swap_performance_win_win', 'metersran_per_min', 'probability_correct'}; 
names2plot = {'win-win', 'av speed (cm/sec)', 'prob. correct'}; 
ylims = {[-.45, .45], [-15, 15], [-.3, .3]};
nv = ceil(sqrt(length(vars2plot)));
x = trainNum;
is_expt = [0 0 0 0 0 1 0 1 1 1 0 ]==1;
is_fem  = [1 1 1 1 1 1 0 0 0 0 0 ]==1;
% clrs = [.7 .7 .7; 1, .3, .3];
marks = {'d', 'o'};

% c = [.7 .7 .7]; e = [1 .4 0];
% rats_cmap = cat(1, c, c, c, c, c, e, c, e, e, e, c, e);
% sess2plot = [14:18];
% sess2plot = [17:21];
% sess2plot = [26:27];
sess2plot = [15 20];
%%%%%% only include rats above chance

rats2include = probability_correct(trainNum==sess2plot(1))>.7;
rats2include = rats2include(1:end-1);

figidx = 0;
% 000
figure(1000*min(sess2plot) + figidx); clf; 
% set(gcf, 'Position', [50+figidx*420 800-30*min(sess2plot) 560 420]); figidx=figidx+1;
set(gcf, 'Position', [50   572   740   198]); figidx=figidx+1;

% sess2plot = [15, 17];
expt_offset_scale  = length(sess2plot)+1;
expt_offset  = is_expt*expt_offset_scale;
xticknames = x;% {'-1', '+1'};
ii = 0;
for i = 1:2
    for j = 1:3
        ii=ii+1;
        if ii<=length(vars2plot)
            %%
%             subplot(2, 3, ii); hold on
            subplot(1,3, ii); hold on
            idx_all = [];
            for aLoop = 1:length(anames)
                idx = ismember(x(aLoop,:), sess2plot);
                if sum(idx)==length(sess2plot)
                    idx_all = cat(1, idx_all, idx);
                    y = eval(sprintf('%s(aLoop,:);', vars2plot{ii}));
%                     plot(expt_offset(aLoop)+x(aLoop,idx), y(idx), 'Color', clrs(is_expt(aLoop)+1, :))
                else
                    idx_all = cat(1, idx_all, idx*0);
                end
                %             scatter(x(aLoop,idx), y(idx), 'Marker', marks{is_expt(aLoop)+1})
            end
            xlim([min(sess2plot)-1, max(sess2plot+max(expt_offset))+1])
            ylim(ylims{ii})
            ylabel(sprintf('%s', names2plot{ii}), 'Interpreter', 'none')
%             set(gca, 'XTick', [sess2plot sess2plot+expt_offset_scale], 'XTickLabel', [sess2plot sess2plot])% 'XTickLabel', xticknames)
            set(gca, 'XTick', [mean(sess2plot) mean(sess2plot)+expt_offset_scale], 'XTickLabel', [sess2plot(1)-sess2plot(1) sess2plot(2)-sess2plot(1)])% 'XTickLabel', xticknames)
            idx_rat = all(idx_all(:, sess2plot+1),2);
            y = eval(sprintf('%s(idx_rat,:);', vars2plot{ii}));
            idx_sess = all(idx_all(idx_rat,:),1);
            y = y(:,idx_sess);
            %%%%%%%%%%%%%%%%%%
%             y = y(rats2include,:);
            %%%%%%%%%%%%%%%%%%
            if size(sess2plot,2) == 2
                yd = y(:,2) - y(:,1);
                y_sal = yd(~is_expt & rats2include');
                y_ibo = yd(is_expt & rats2include');
                plot([-100 100], [0 0], 'Color', 'k', 'LineWidth', 1)
                plot([mean(sess2plot-.5) mean(sess2plot)+.5], [mean(y_sal) mean(y_sal)], 'Color', clrs(1, :)/2, 'LineWidth', 2)
                plot([mean(sess2plot-.5) mean(sess2plot)+.5]+expt_offset_scale, [mean(y_ibo) mean(y_ibo)], 'Color', clrs(2, :)/2, 'LineWidth', 2)
                jitter = (rand(length(y_sal),1)-.5)/4;
                scatter(jitter + mean(sess2plot), y_sal, 30, 'MarkerFaceColor', clrs(1, :)/2, 'MarkerEdgeColor', 'k')
                jitter = (rand(length(y_ibo),1)-.5)/4;
                scatter(jitter + mean(sess2plot)+expt_offset_scale, y_ibo, 30, 'MarkerFaceColor', clrs(2, :)/2, 'MarkerEdgeColor', 'k')
                
                for loglog = [0,1]
                    p = ranksum(y(is_expt==loglog & rats2include',1), y(loglog==is_expt & rats2include',2));
                    %             [h,p] = ttest2(y(is_expt==loglog,1), y(loglog==is_expt,2));
                    [p,tbl,stats] = kruskalwallis([y_sal; y_ibo], [y_sal*0; y_ibo*0+1], 'off');
                    plotp = true;%%%
                    if p<.001
                        str = '***';
                    elseif p<.01
                        str = '**';
                    elseif p<.05
                        str = '*';
                    else
                        str = '~';
                        plotp = true;
                    end
                    if plotp == true
                        text(mean(sess2plot)+mean(mean(expt_offset(is_expt==loglog))), ylims{ii}(2)*.9, sprintf('%s', str))
                        text(mean(sess2plot)+mean(mean(expt_offset(is_expt==loglog))), ylims{ii}(2)*.85, sprintf('%0.3f', p));end
                end
            else
%                 plot(sess2plot, mean(y(~is_expt & rats2include',:),1), 'Color', clrs(1, :)/2, 'LineWidth', 2)
%                 plot(sess2plot+expt_offset_scale, mean(y(is_expt & rats2include',:),1), 'Color', clrs(2, :)/2, 'LineWidth', 2)
                
                for loglog = [0,1]
                    p = ranksum(y(is_expt==loglog & rats2include',1), y(loglog==is_expt & rats2include',2));
                    %             [h,p] = ttest2(y(is_expt==loglog,1), y(loglog==is_expt,2));
                    plotp = true;%%%
                    if p<.001
                        str = '***';
                    elseif p<.01
                        str = '**';
                    elseif p<.05
                        str = '*';
                    else
                        str = '~';
                        plotp = false;
                    end
                    if plotp == true
                        text(mean(sess2plot)+mean(mean(expt_offset(is_expt==loglog))), ylims{ii}(2)*.9, sprintf('%s', str))
                        text(mean(sess2plot)+mean(mean(expt_offset(is_expt==loglog))), ylims{ii}(2)*.85, sprintf('%0.3f', p));end
                end
                
            end
            set(gca, 'Box', 'off')
        end
    end
end
% ii=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 001
figure(1000*min(sess2plot) + figidx); clf; set(gcf, 'Position', [50+figidx*420 800-30*min(sess2plot) 560 420]); figidx=figidx+1;
% sess2plot = [15, 17];
abs_zone_sampling_angdiff = abs(zone_sampling_angdiff);
probability_correct_diff = probability_correct - probability_correct(:,16);
vars2plot = {'abs_zone_sampling_angdiff', 'zone_sampling_selectivity', 'probability_correct'}; 
names2plot = {'angdiff', 'select.', 'prob. correct'}; 

nv = (length(vars2plot));
% x = trainNum;
x = trainNum;

ylims = {[0, 200], [-.1 1], [.2, 1.0]};

% sess2plot = [15, 17:22];
% sess2plot = [22:26];
expt_offset_scale  = length(sess2plot)+1;
expt_offset  = is_expt*expt_offset_scale;
xticknames = x; % {'-1', '+1'};
ii = 0;
for i = 1:nv
    ii=ii+1;
    if ii<=length(vars2plot)
        subplot(nv, 1, i); hold on
        idx_all = [];
        for aLoop = 1:length(anames)
            idx = ismember(x(aLoop,:), sess2plot);
            if sum(idx)==length(sess2plot)
                idx_all = cat(1, idx_all, idx);
                y = eval(sprintf('%s(aLoop,:);', vars2plot{ii}));
                plot(expt_offset(aLoop)+x(aLoop,idx), y(idx), 'Color', clrs(is_expt(aLoop)+1, :))
            else
                idx_all = cat(1, idx_all, idx*0);
            end
            %             scatter(x(aLoop,idx), y(idx), 'Marker', marks{is_expt(aLoop)+1})
        end
        xlim([min(sess2plot)-.5, max(sess2plot+max(expt_offset))+.5])
        ylim(ylims{ii})
        ylabel(sprintf('%s', names2plot{ii}), 'Interpreter', 'none')
        set(gca, 'XTick', [sess2plot sess2plot+expt_offset_scale], 'XTickLabel', [sess2plot sess2plot])
        idx_rat = all(idx_all(:, sess2plot+1),2);
        y = eval(sprintf('%s(idx_rat,:);', vars2plot{ii}));
        idx_sess = all(idx_all(idx_rat,:),1);
        y = y(:,idx_sess);
        plot(sess2plot, mean(y(~is_expt,:),1), 'Color', clrs(1, :)/2, 'LineWidth', 2)
        plot(sess2plot+expt_offset_scale, mean(y(is_expt,:),1), 'Color', clrs(2, :)/2, 'LineWidth', 2)
        
        for loglog = [0,1]
            %             p = ranksum(y(is_expt==loglog,1), y(loglog==is_expt,2));
            [h,p] = ttest2(y(is_expt==loglog,1), y(loglog==is_expt,2));
            plotp = false;%%%
            if p<.001
                str = '***';
            elseif p<.01
                str = '**';
            elseif p<.05
                str = '*';
            else
                str = '~';
                plotp = false;
            end
            if plotp;
                text(mean(sess2plot)+mean(mean(expt_offset(is_expt==loglog))), ylims{ii}(2)*.9, sprintf('%s', str))
                text(mean(sess2plot)+mean(mean(expt_offset(is_expt==loglog))), ylims{ii}(2)*.85, sprintf('%0.3f', p));end
        end
    end
end

figure(1000*min(sess2plot) + figidx); clf; set(gcf, 'Position', [50+figidx*420 800-30*min(sess2plot) 560 420]); figidx=figidx+1;
subplot(3,1,1);
hold on
d = zone_sampling_angdiff;% - zone_sampling_angdiff(:, 16);
s = zone_sampling_selectivity;% - zone_sampling_selectivity(:, 16);
p = probability_correct;% - probability_correct(:, 16);
x = trainNum;
plot(x(~is_expt,:)', abs(d(~is_expt,:))', 'k')
plot(x(is_expt,:)', abs(d(is_expt,:))', 'r')
% xlim([nanmin(x(:))-.5, nanmax(x(:))+.5 ])
xlim([nanmin(sess2plot)-.5, nanmax(sess2plot)+.5 ])

subplot(3,1,2);
hold on
plot(x(~is_expt,:)', (s(~is_expt,:))', 'k')
plot(x(is_expt,:)', (s(is_expt,:))', 'r')
xlim([nanmin(sess2plot)-.5, nanmax(sess2plot)+.5 ])

subplot(3,1,3);
hold on
plot(x(~is_expt,:)', (p(~is_expt,:))', 'k')
plot(x(is_expt,:)', (p(is_expt,:))', 'r')
xlim([nanmin(sess2plot)-.5, nanmax(sess2plot)+.5 ])

%%%%%%%%%%%%%%%%%%%%%%%%%%

% sess2plot = [22:27];
x = trainNum;
figure(1000*min(sess2plot) + figidx); clf; set(gcf, 'Position', [50+figidx*420 800-30*min(sess2plot) 560 420]); figidx=figidx+1;
pxy_all = [];
for sn = 1:length(sess2plot)
%     set(gcf, 'Position', [100+a*100   442   850   536], 'Name', anames{a})
    hold on
    y_exp = [];
    y_con = [];
        s = sess2plot(sn);
        pxy = [];
        for aLoop = 1:numan-1
        fname = filenames(aLoop,s+1);
        if isfile(fname)
            subplot(2,length(sess2plot),sn)
            d=load(fname);
            [d.zonestruct] = DUAL8MAZE_performance_eval(d.room, d.arena, params, false);
            binedge = d.zonestruct.metrics.sampling_room_arena_bins;%    linspace(-pi, pi, 37);
            bincenters = binedge(2:end) - mean(abs(diff(binedge)))/2;
            y = d.zonestruct.metrics.sampling_room; 
%             bincenters = [bincenters, bincenters(1)];
%             y = [y, y(1)];
            if is_expt(aLoop)
            y_exp = cat(1, y_exp, y);
            else
            y_con = cat(1, y_con, y);
            end
%             polarplot(bincenters, y, 'Color', clrs(is_expt(aLoop)+1, :))
            plot(rad2deg(bincenters), y, 'Color', clrs(is_expt(aLoop)+1, :))
            hold on;
            pxy = cat(2, pxy, d.zonestruct.metrics.sampling_room_xy);
            
%             figure;
%             set(gcf, 'Position', [680   802   641   176])
%             im = [a_pxy(:,1:4)*NaN, r_pxy*-1, a_pxy(:,1:4)*NaN, a_pxy, a_pxy(:,1:4)*NaN];
%             im = [im(1:4, :)*NaN; im; im(1:4, :)*NaN];
%             subplot(1,2,1);
%             h = imagesc(im, [-.8 .8]); colorbar
%             axis image
%             colormap(lbmap(100,'RedBlue'));
%             set(h, 'AlphaData', ~isnan(im))
%             subplot(1,2,2); hold on;
%             % a_pt = [a_pt(end) a_pt];
%             % r_pt = [r_pt(end) r_pt];
%             % xs = [xs(end) xs];
%             plot(xs, r_pt, 'r');
%             plot(xs, a_pt, 'b')
            set(gca, 'XTick', [-180, -90 0 90 180])
            ylim([0 1])
            drawnow
%             pause(1)
        else
            pxy = cat(2, pxy, d.zonestruct.metrics.sampling_room_xy);
        end
        end
        pxy_all = cat(1, pxy_all, pxy);
%     shadedErrorBar(rad2deg(bincenters), median(y_con,1), std(y_con,1), 'LineProps', {'-', 'Color', clrs(1, :)/2, 'LineWidth', 2})
%     shadedErrorBar(rad2deg(bincenters), median(y_exp,1), std(y_exp,1), 'LineProps', {'-', 'Color', clrs(2, :)/2, 'LineWidth', 2})
    plot(rad2deg(bincenters), median(y_exp,1), '-', 'Color', clrs(2, :)/2, 'LineWidth', 2)
    plot(rad2deg(bincenters), median(y_con,1), '-', 'Color', clrs(1, :)/2, 'LineWidth', 2)
%     polarplot(bincenters, mean(y_exp,1), '-', 'Color', clrs(2, :)/2, 'LineWidth', 2)
%     polarplot(bincenters, mean(y_con,1), '-', 'Color', clrs(1, :)/2, 'LineWidth', 2)
%     rlim([0, 1]);
%     set(gca, 'RTick', [.25:.25:1], 'ThetaTick', [0:30:330])
    title(sprintf('%s_%d', sesstype(1, sess2plot(sn)+1), sess2plot(sn)), 'Interpreter', 'none')

end
% subplot(2,length(sess2plot),1)
% title('last sess, 30 min')
% subplot(2,length(sess2plot),2)
% title('retrieval, 10 min')
% subplot(2,length(sess2plot),3)
% title('next sess, 30 min')
% subplot(2,length(sess2plot),3)
% title('+2 sess, 30 min')

subplot(2,1,2); hold on
    y_exp = [];
    y_con = [];
for aLoop = 1:length(anames)-1
    idx = ismember(x(aLoop,:), sess2plot);
%     y = right_preference_room(aLoop,:); %eval(sprintf('%s(aLoop,:);', vars2plot{ii}));
    y = probability_correct(aLoop,:); %eval(sprintf('%s(aLoop,:);', vars2plot{ii}));
    plot(x(aLoop,idx), y(idx), 'Color', clrs(is_expt(aLoop)+1, :))
            if is_expt(aLoop)
            y_exp = cat(1, y_exp, y(idx));
            else
            y_con = cat(1, y_con, y(idx));
            end
%             scatter(x(aLoop,idx), y(idx), 'Marker', marks{is_expt(aLoop)+1})
end
xlim([min(sess2plot)-.5, max(sess2plot+.5)])
ylim([.2 1])
    plot(sess2plot, mean(y_con,1), 'Color', clrs(1, :)/2, 'LineWidth',2)
    plot(sess2plot, mean(y_exp,1), 'Color', clrs(2, :)/2, 'LineWidth',2)
ylabel(sprintf('%s', names2plot{3}), 'Interpreter', 'none')
% set(gca, 'XTick', sess2plot, 'XTickLabel', {'-1' 'RET' '+1'})
title('Reinforced area sampling liklihood')


%%
figure(3); clf;
x = trainNum;
% distbins = linspace(0, pi, 5);
distbins_center = rad2deg( distbins(1:end-1) + mean(abs(diff(distbins))) );
% sess2plot = [15, 17 18];
for sLoop = 1:length(sess2plot)
    y_exp = [];
    y_con = [];
    subplot(1,length(sess2plot),sLoop); hold on
    for aLoop = 1:length(anames)-1
        idx = ismember(x(aLoop,:), sess2plot(sLoop));
        y = squeeze(correct_by_dist(aLoop,idx ,:))'; %eval(sprintf('%s(aLoop,:);', vars2plot{ii}));
%         y=y./sum(y);
        plot(distbins_center, y, 'Color', clrs(is_expt(aLoop)+1, :))
        if is_expt(aLoop)
            y_exp = cat(1, y_exp, y);
        else
            y_con = cat(1, y_con, y);
        end
    end
xlim([0 225])
ylim([-.1 1.1])
set(gca, 'XTick', distbins_center, 'YTick', [0:.25:1])

plot(distbins_center, mean(y_con,1), 'Color', clrs(1, :)/2, 'LineWidth',2)
plot(distbins_center, mean(y_exp,1), 'Color', clrs(2, :)/2, 'LineWidth',2)
end

%%
% sess2plot = [15, 17];
sess2plot = [22 23];
figure(5); clf
for sn = 1:length(sess2plot)
%     set(gcf, 'Position', [100+a*100   442   850   536], 'Name', anames{a})
    hold on
    y_exp = [];
    y_con = [];
        s = sess2plot(sn);
        for aLoop = 1:numan-1
        fname = filenames(aLoop,s+1);
        if isfile(fname)
            subplot(1,length(sess2plot),sn)
            d=load(fname);
            prop=1/2;
            [outvars] = SMART_useful_figs(d.zonestruct, 1/prop, 300, true, false, pi/8);
%             y1 = outvars.correct_runav_time;
%             y2 = outvars.left_runav_time;
%             x  = outvars.correct_runav_time_xv/60;
            y1 = outvars.correct_runav_prop;
            y2 = outvars.left_runav_prop;
%             x  = outvars.correct_runav_prop_xv;
            x = prop:prop:1;
            binedge = d.zonestruct.metrics.sampling_room_arena_bins;%    linspace(-pi, pi, 37);
            bincenters = binedge(2:end) - mean(abs(diff(binedge)));
            bincenters = [bincenters, bincenters(1)];
            y = d.zonestruct.metrics.sampling_room(1,:); y = [y, y(1)];
%             x = bincenters;
%             y1 = y';
            if is_expt(aLoop)
            y_exp = cat(1, y_exp, y1');
            else
            y_con = cat(1, y_con, y1');
            end
            plot(x, y1, 'Color', clrs(is_expt(aLoop)+1, :))
            hold on;
            
        end
    end
    shadedErrorBar(x, mean(y_exp,1), std(y_exp,1), 'LineProps', {'Color', clrs(2, :)/2, 'LineWidth', 2})
    shadedErrorBar(x, mean(y_con,1), std(y_con,1), 'LineProps', {'Color', clrs(1, :)/2, 'LineWidth', 2})
%     plot(x, mean(y_exp,1), '-', 'Color', clrs(2, :)/2, 'LineWidth', 2)
%     plot(x, mean(y_con,1), '-', 'Color', clrs(1, :)/2, 'LineWidth', 2)
%     axis([-1 41 .4 1.1])
    axis([min(x)-mean(abs(diff(x)))/2 max(x)+mean(abs(diff(x)))/2 0 1.1])
    set(gca, 'XTick', unique(x), 'YTick', [0:.25:1])
end












































%%
figure; 
subplot(2,2,1)
hold on; plot(dayNum', probability_correct', 'k')
plot(dayNum([6 8 9 10],:)', probability_correct([6 8 9 10],:)', 'r')
ylabel('probability correct')

subplot(2,2,2) 
hold on; plot(dayNum', average_streak', 'k'); % right_preference_room
plot(dayNum([6 8 9 10],:)', average_streak([6 8 9 10],:)', 'r')
ylabel('av streak')

subplot(2,2,3)
hold on; plot(dayNum', entrpermin', 'k')
plot(dayNum([6 8 9 10],:)', entrpermin([6 8 9 10],:)', 'r')
ylabel('entries / min')

subplot(2,2,4)
hold on; plot(dayNum', lapsPerMin', 'k')
plot(dayNum([6 8 9 10],:)', lapsPerMin([6 8 9 10],:)', 'r')
ylabel('laps / min')

figure; 
subplot(1,2,1)
hold on; plot(dayNum', swap_performance_lose_win', 'k')
plot(dayNum([6 8 9 10],:)', swap_performance_lose_win([6 8 9 10],:)', 'r')
ylabel('lose win')
subplot(1,2,2)
hold on; plot(dayNum', swap_performance_win_win', 'k')
plot(dayNum([6 8 9 10],:)', swap_performance_win_win([6 8 9 10],:)', 'r')
ylabel('win win')

swap_performance_lose_lose = squeeze(swap_performance_all(:,:,2,2));
swap_performance_win_lose  = squeeze(swap_performance_all(:,:,1,2));
%%
% figure; plot(trainNum', probability_correct', 'b'); hold on; plot(trainNum', 1-probability_left', 'k')
% figure; plot(trainNum', swap_performance_win_win', 'b'); hold on; plot(trainNum', swap_performance_lose_win', 'k')
% figure; plot(trainNum', average_streak'); % hold on; plot(average_err); 

a = correct_by_dist(correct_by_dist(:,1)>.6,2:end);
b = correct_by_dist(correct_by_dist(:,1)<=.6,2:end);
figure; hold on;
plot(a', 'b'); % hold on; plot(average_err); 
plot(b', 'k'); % hold on; plot(average_err); 
plot(nanmean(a, 1), 'b', 'LineWidth' ,2)
plot(nanmean(b, 1), 'k', 'LineWidth' ,2)
drawnow

fprintf('\n\n');
n=3; aa = NaN(numan, n); bb = NaN(numan, n);
for i = 1:numan
    a = probability_correct(i,:);
    b = entrpermin(i,:);
    idx = find(~isnan(a), n, 'last');
    fprintf('%s-', anname(i,idx(1)));
    for j = 1:length(idx)
        fprintf('\t\t%s  %3.2f%%, %1.2f,', sessname(i, idx(j)), round(10000*a(idx(j)))/100, b(idx(j)));
    end
    aa(i,1:length(idx)) = a(idx); 
    bb(i,1:length(idx)) = b(idx);
    fprintf('\n');
end
%% separate groups, evens vs odds
sessn=3; % 
anova1([aa(1:2:end, sessn); aa(2:2:end, sessn)], [0*aa(1:2:end, sessn); 1+0*aa(2:2:end, sessn)])
hold on
plot(aa(2:2:end, sessn)*0 +1.75+ rand(6,1)/2, aa(2:2:end, sessn), 'ro')
plot(aa(1:2:end, sessn)*0 +0.75+ rand(6,1)/2, aa(1:2:end, sessn), 'bo')

anova1([bb(1:2:end, sessn); bb(2:2:end, sessn)], [0*bb(1:2:end, sessn); 1+0*bb(2:2:end, sessn)])
hold on;
plot(bb(2:2:end, sessn)*0 +1.75+ rand(6,1)/2, bb(2:2:end, sessn), 'ro')
plot(bb(1:2:end, sessn)*0 +0.75+ rand(6,1)/2, bb(1:2:end, sessn), 'bo')

s1 = [1:12]; s2 = [6,8,9,10,12]; s1=s1(~ismember(s1,s2));
anova1([aa(s1,sessn); aa(s2,sessn)], [0*aa(s1,sessn); 1+0*aa(s2,sessn)])
hold on
plot(aa(s2,sessn)*0 +1.75+ rand(length(s2),1)/2, aa(s2,sessn), 'ro')
plot(aa(s1,sessn)*0 +0.75+ rand(length(s1),1)/2, aa(s1,sessn), 'bo')

anova1([bb(s1,sessn); bb(s2,sessn)], [0*bb(s1,sessn); 1+0*bb(s2,sessn)])
hold on
plot(bb(s2,sessn)*0 +1.75+ rand(length(s2),1)/2, bb(s2,sessn), 'ro')
plot(bb(s1,sessn)*0 +0.75+ rand(length(s1),1)/2, bb(s1,sessn), 'bo')



figure; 
subplot(2,1,1); hold on; 
plot(aa'); axis([0.5 n+.5 0 1]); legend(anames, 'Location', 'eastoutside'); 
plot([0 n+1], [.5 .5], 'k:'); plot([0 n+1], [.7 .7], 'k:');

subplot(2,1,2); hold on; 
plot(bb'); axis([0.5 n+.5 0 max(bb(:))*1.1]); legend(anames, 'Location', 'eastoutside'); 
plot([0 n+1], [2 2], 'k:'); plot([0 n+1], [1 1], 'k:');

distbins = linspace(0, pi, 9);
% h = binned_statistic1d(zonestruct.choices.choiceDist, zonestruct.choices.isCorrect, distbins, 'mean');
h = binned_statistic1d(all_h(:,1), all_h(:,2), distbins, 'nanmean');
hh = binned_statistic1d(all_h(:,1), all_h(:,2), distbins, 'nanstd');
h2 = binned_statistic1d(all_h2(:,1), all_h2(:,2), distbins, 'nanmean');
hh2 = binned_statistic1d(all_h2(:,1), all_h2(:,2), distbins, 'nanstd');
n1 = sum(all_h(~isnan(all_h(:,1)), 1)>0);
n2 = sum(all_h2(~isnan(all_h2(:,1)), 1)>0);

figure; hold on
plot(distbins(2:end), h, 'b')
plot(distbins(2:end), h2, 'k')
% shadedErrorBar(distbins(2:end), h, hh,  'lineprops', {'b'}); 
% shadedErrorBar(distbins(2:end), h2, hh2,  'lineprops', {'k'}); 
scatter(distbins(2:end), h, 'bd'); 
scatter(distbins(2:end), h2, 'ko'); 
% scatter(zonestruct.choices.choiceDist, zonestruct.choices.isCorrect, 'ko')