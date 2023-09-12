% A_names = {'BM000','BM027','R1','R2','R3','R4','R5','R6','R7'};
A_names = {'BM000','BM027','R1','R2','R3','R4','R6','R7'}';
A_num = length(A_names);
shock_sess = ... % shock sessions
   [5   16;...     % BM000
    10   0;...     % BM027
    17   21;...     % R1
    9   13;...     % R2
    6   10;...    % R3
    9   13;...     % R4
    8   12;...    % R6
    12  16];       % R7
is_female    = [0 1 1 0 1 0 0 1]'==1;
scopo_first = [0 0 1 1 1 1 0 0]'==1;

max_sess = 30;
% f_temp = 'D:\Sample Data\Linear Maze Data\Blair2022_SuppBehavior\BM000\Session3\BM000_Sess3.mat';
fileDir = 'D:\Sample Data\Linear Maze Data\Blair2022_SuppBehavior\';
% SETUP FILES
exp.files       = cell(1,1);
exp.saveName    = cell(1,1);
exp.names       = cell(1,1);
exp.animalidx 	= NaN(1,1);
exp.sessnum     = NaN(1,1);
f_idx = 0;
fprintf('\n~~~~~');
for aLoop = 1:A_num
    fprintf('\nGetting files for animal %s...\t ', A_names{aLoop});
    sessind = 0;
    makeFiles = true;
    while makeFiles == true
        sessind = sessind+1;
        f_temp = sprintf('%s%s\\Session%d\\%s_Sess%d.mat', fileDir, A_names{aLoop}, sessind, A_names{aLoop}, sessind);
        if ~isfile(f_temp)
            temp = dir(sprintf('%s%s\\Session%d*', fileDir, A_names{aLoop}, sessind));
            if ~isempty(temp)
            weirdname = [temp.folder '\' temp.name];
            f_temp = sprintf('%s\\%s_Sess%d.mat', weirdname, A_names{aLoop}, sessind);
            end
        end
        if isfile(f_temp)
            fprintf('%d ', sessind);
            saveName = sprintf('%s%s\\%s_behav_Sess%d.mat', fileDir, A_names{aLoop}, A_names{aLoop}, sessind);
            f_idx = f_idx+1;
            exp.files(f_idx,1) = {f_temp};
            exp.names(f_idx,1) = A_names(aLoop);
            exp.saveName(f_idx,1) = {saveName};
            exp.animalidx(f_idx,1) = aLoop;
            exp.sessnum(f_idx,1) = sessind;
        elseif max_sess<sessind
            makeFiles = false;
        end
    end
    fprintf('Done!');
end
fprintf('\n~~~~~\n');

%% ARRANGING, SEGMENTING, and SAVING BEHAVIOR
Efects.shortrate_pre10 = NaN(max(exp.animalidx), max(exp.sessnum));
Efects.shortlap_prop = NaN(max(exp.animalidx), max(exp.sessnum));
% figure(205); clf
for fileLoop = 1:length(exp.files)
    %%
    behav = [];
    behav.animal = exp.names{fileLoop};
    behav.origFile = exp.files{fileLoop};
    behav.behavFile = exp.saveName{fileLoop};
    behav.sessnum = exp.sessnum(fileLoop);
    behav.animalIdx = exp.animalidx(fileLoop);
    
    temp = load(behav.origFile);
    t = temp.VT_time./1000;
    good_ind = ~isnan(t);%t<=600;%
    x = temp.xpos(good_ind)';
    y = temp.ypos(good_ind)';
    t = t(good_ind)';
    
    xo = x; yo = y;
    xcenter = 170;%175;%350; % center x,y point on the short path
    ycenter = 26;%150; % center x,y point on the short path
    xo = xo-xcenter;
    yo = yo-ycenter;
    x = 250*xo/(285-50);%(295-43);%440;
    y = 125*yo/(150-22);%(149-10);
    % figure; plot3(x,y, t);%1:length(behav.time))
    behav.x = x;
    behav.y = y;
    behav.time = t;
    
    behav.dt = [median(diff(behav.time)); diff([behav.time])];
    behav.spd = sqrt(diff([x(1); x]).^2 + diff([y(1); y]).^2)./behav.dt;
    behav.spd(1) = behav.spd(2);
    
    behav.always_check_return = true;
    plotting = false;
    
    [behav.SR, behav.SL, behav.LR, behav.LL, ~, ~] = linear_segmentation_ScaryMaze_LFOV(behav.x, behav.y, behav.spd, median(behav.dt), behav.always_check_return, plotting);
    pre10_sr = [];
    pre10_sl = [];
    pre10_lr = [];
    pre10_ll = [];
    pre10_sr_spd = [];
    pre10_sl_spd = [];
    pre10_lr_spd = [];
    pre10_ll_spd = [];
    if ~isnan(behav.SR.pass_ind(1))
    pre10_sr = behav.time(behav.SR.pass_ind(:,2))<=600;
    pre10_sr_spd = behav.SR.mean_spd(  behav.time(behav.SR.pass_ind(:,2))<=600 );
    end
    if ~isnan(behav.SL.pass_ind(1))
    pre10_sl = behav.time(behav.SL.pass_ind(:,2))<=600;
    pre10_sl_spd = behav.SL.mean_spd(  behav.time(behav.SL.pass_ind(:,2))<=600 );
    end
    if ~isnan(behav.LR.pass_ind(1))
    pre10_lr = behav.time(behav.LR.pass_ind(:,2))<=600;
    pre10_lr_spd = behav.LR.mean_spd(  behav.time(behav.LR.pass_ind(:,2))<=600 );
    end
    if ~isnan(behav.LL.pass_ind(1))
    pre10_ll = behav.time(behav.LL.pass_ind(:,2))<=600;
    pre10_ll_spd = behav.LL.mean_spd(  behav.time(behav.LL.pass_ind(:,2))<=600 );
    end
    
    behav.numshort = length(pre10_sr) + length(pre10_sl);
    behav.medianshortspeed_pre10 = nanmedian([pre10_sr_spd; pre10_sl_spd])
    behav.shortrate_pre10 = (sum(pre10_sr) + sum(pre10_sl))/10; % 10 mins
    behav.shortrate_post10 = (sum(~pre10_sr) + sum(~pre10_sl))/5; % last 5 min

    behav.numlong = length(pre10_lr) + length(pre10_ll);
    behav.longrate_pre10 = (sum(pre10_lr) + sum(pre10_ll))/10; % 10 mins
    behav.longrate_post10 = (sum(~pre10_lr) + sum(~pre10_ll))/5; % last 5 min
    
    r = behav.shortrate_pre10;
    p = behav.shortrate_pre10/(behav.longrate_pre10+behav.shortrate_pre10);
    
    Efects.shortrate_pre10(behav.animalIdx, behav.sessnum) = r;
    Efects.shortlap_prop(behav.animalIdx, behav.sessnum) = p;
    
    xx = [behav.sessnum-.5, behav.sessnum+.5];
    yy1 = [r, r];
    yy2 = [p p];
%     figure(205);
%     subplot(A_num, 2, 2*(behav.animalIdx-1)+1); hold on
% %     plot3(xx, yy1, yy2); axis([0 max_sess 0 5 0 1])
%     plot([0 max_sess], [2 2], 'r')
%     plot(xx, yy1); axis([0 max_sess 0 5])
%     
%     subplot(A_num, 2, 2*(behav.animalIdx-1)+2); hold on
%     plot([0 max_sess], [2/3 2/3], 'r')
%     plot(xx, yy2);  axis([0 max_sess 0 1])
    % short path rewards / min - immediate effect
    % short path rewards / min - 8 hr retention
    % short path rewards / min
    save(exp.saveName{fileLoop}, 'behav');
end


%% LOAD BEHAVIOR and EVALUATE
% A_names = {'BM000','BM027','R1','R2','R3','R4','R6','R7'};
% shock_sess = ... % shock sessions
%    [5   0;...     % BM000
%     0   0;...     % BM027
%     0   0;...     % R1
%     9   0;...     % R2
%     6   10;...    % R3
%     9   0;...     % R4
%     8   12;...    % R6
%     12  0];       % R7

shock_sess_diff = [-2:2];%-2:2; %-1:1;%
shortrate_scopo = NaN( A_num, length(shock_sess_diff) );
shortspd_scopo = NaN( A_num, length(shock_sess_diff) );
shortprop_scopo = NaN( A_num, length(shock_sess_diff) );
shortrate_saline = NaN( A_num, length(shock_sess_diff) );
shortspd_saline = NaN( A_num, length(shock_sess_diff) );
shortprop_saline = NaN( A_num, length(shock_sess_diff) );

shortrate_scopo_postshock = NaN( A_num, length(shock_sess_diff) );
shortprop_scopo_postshock = NaN( A_num, length(shock_sess_diff) );
shortrate_saline_postshock = NaN( A_num, length(shock_sess_diff) );
shortprop_saline_postshock = NaN( A_num, length(shock_sess_diff) );
% shock_sess_diff = NaN( A_num, length(shock_sess_diff) );

for aLoop = 1:A_num
    if any(shock_sess(aLoop, :)>0)
        for shockLoop = 1:sum(shock_sess(aLoop, :)>0)
            shk_sess = shock_sess(aLoop, shockLoop);
            for sessLoop = 1:length(shock_sess_diff) % shk_sess+shock_sess_diff
                %%
                this_sess = shk_sess + shock_sess_diff(sessLoop);
%                 shock_sess_diff(aLoop, sessLoop) = this_sess - shk_sess;
                session_ind = find(exp.animalidx==aLoop & exp.sessnum==this_sess);
                if ~isempty(session_ind)
                    %%
                fname = exp.saveName{session_ind};
                behav = [];
                load(fname, 'behav');
                r = behav.shortrate_pre10;
                s = behav.medianshortspeed_pre10;
                p = behav.shortrate_pre10/(behav.longrate_pre10+behav.shortrate_pre10);
                t = behav.shortrate_post10;
                g = behav.shortrate_post10/(behav.longrate_post10+behav.shortrate_post10);
                if scopo_first(aLoop)==1 && shockLoop==1
                    shortrate_scopo(aLoop, sessLoop) = r;
                    shortspd_scopo(aLoop, sessLoop) = s;
                    shortprop_scopo(aLoop, sessLoop) = p;                    
                    shortrate_scopo_postshock(aLoop, sessLoop) = t;
                    shortprop_scopo_postshock(aLoop, sessLoop) = g;                    
                elseif scopo_first(aLoop)==0 && shockLoop==2
                    shortrate_scopo(aLoop, sessLoop) = r;
                    shortspd_scopo(aLoop, sessLoop) = s;
                    shortprop_scopo(aLoop, sessLoop) = p;
                    shortrate_scopo_postshock(aLoop, sessLoop) = t;
                    shortprop_scopo_postshock(aLoop, sessLoop) = g;
                else
                    shortrate_saline(aLoop, sessLoop) = r;
                    shortspd_saline(aLoop, sessLoop) = s;
                    shortprop_saline(aLoop, sessLoop) = p;
                    shortrate_saline_postshock(aLoop, sessLoop) = t;
                    shortprop_saline_postshock(aLoop, sessLoop) = g;
                end
% %                 shortrate_scopo(shockLoop, aLoop, sessLoop) = r;
% %                 shortrate_scopo(shockLoop, aLoop, sessLoop) = p;
%                 figure; plot3(behav.x, behav.y, behav.time)
%                 figname = sprintf('%s Sess%d, \nrate=%1.1f, prop=%0.2f', A_names{aLoop}, this_sess, r, p);
%                 title(figname)
%                 figname = sprintf('%s Sess%d, \nrate=%1.1f, prop=%0.2f', A_names{aLoop}, this_sess, r, p);
%                 set(gcf, 'Name', sprintf('%s Sess%d', A_names{aLoop}, this_sess))
                end
            end
            
        end
    end
end
% shock_sess_diff = nanmean(shock_sess_diff,1);


% s2 = squeeze(shortrate(2,any(shock_sess>0,2),:));
presess = shock_sess_diff==-1;
postsess = shock_sess_diff==1;

var_names = {'animal_id', 'is_female', 'scopo_shock_first', 'SAL_pre', 'SAL_post', 'SCOP_pre', 'SCOP_post'};
saline_pre = shortrate_saline(:, presess);
saline_post = shortrate_saline(:, postsess);
scopo_pre = shortrate_scopo(:, presess);
scopo_post = shortrate_scopo(:, postsess);
Short_rate = table(A_names, is_female, scopo_first, saline_pre, saline_post, scopo_pre, scopo_post, 'VariableNames', var_names);

saline_pre = shortprop_saline(:, presess);
saline_post = shortprop_saline(:, postsess);
scopo_pre = shortprop_scopo(:, presess);
scopo_post = shortprop_scopo(:, postsess);
Short_proportion = table(A_names, is_female, scopo_first, saline_pre, saline_post, scopo_pre, scopo_post, 'VariableNames', var_names);
vars2save = {'shock_sess_diff' 'Short_rate' 'Short_proportion'...
    'shortrate_saline' 'shortrate_scopo' 'shortprop_saline' 'shortprop_scopo'...
    'shortrate_saline_postshock' 'shortrate_scopo_postshock' 'shortprop_saline_postshock' 'shortprop_scopo_postshock'};
save('D:\Sample Data\Linear Maze Data\Blair2022_SuppBehavior\SummaryData_Supp_postshock.mat', vars2save{:})
% save('C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\SummaryData_Supp_postshock.mat', vars2save{:})

%% FROM REVIEWS 6_18_2023
SAL_shortrate_preshock = shortrate_saline(:, shock_sess_diff==0);
SAL_shortprop_preshock = shortprop_saline(:, shock_sess_diff==0);
SCOP_shortrate_preshock = shortrate_scopo(:, shock_sess_diff==0);
SCOP_shortprop_preshock = shortprop_scopo(:, shock_sess_diff==0);

SAL_shortrate_postshock = shortrate_saline_postshock(:, shock_sess_diff==0);
SAL_shortprop_postshock = shortprop_saline_postshock(:, shock_sess_diff==0);
SCOP_shortrate_postshock = shortrate_scopo_postshock(:, shock_sess_diff==0);
SCOP_shortprop_postshock = shortprop_scopo_postshock(:, shock_sess_diff==0);
figure; hold on; 
plot([SAL_shortrate_preshock SAL_shortrate_postshock]', 'k')
plot([SCOP_shortrate_preshock SCOP_shortrate_postshock]', 'r')
axis([0 3 -1 6])
%%
figure(207); clf;
pre_trn_post_sess = [2,3,4];
subplot(1,2,1); hold on
bar([1,2,3], nanmean(shortspd_saline(:,pre_trn_post_sess), 1), 'FaceColor', [.5 .5 .5])
plot(shortspd_saline(is_female,pre_trn_post_sess)', 'k:o')
plot(shortspd_saline(~is_female,pre_trn_post_sess)', 'k-o')
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'pre' 'trn' 'post'})
axis([0 4 0 160])
subplot(1,2,2); hold on
bar([1,2,3], nanmean(shortspd_scopo(:,pre_trn_post_sess), 1), 'FaceColor', [1 .5 .6])
plot(shortspd_scopo(~is_female,pre_trn_post_sess)', 'm-o')
plot(shortspd_scopo(is_female,pre_trn_post_sess)', 'm:o')
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'pre' 'trn' 'post'})
axis([0 4 0 160])
vars2save = {'shortspd_saline' 'shortspd_scopo' 'is_female'};
save('D:\Sample Data\Linear Maze Data\Blair2022_SuppBehavior\SummaryData_Supp_runningspeed.mat', vars2save{:})

%%
A_names = Short_rate.animal_id; %%%%%%%%%%%
is_female = Short_rate.is_female; %%%%%%%%%%%
scopo_first = Short_rate.scopo_shock_first; %%%%%%%%%%%
A_num = length(A_names);

all_sess_done = ~any(isnan([Short_rate.SAL_pre Short_rate.SAL_post Short_rate.SCOP_pre Short_rate.SCOP_post]), 2);
sal_done = ~any(isnan([Short_rate.SAL_pre Short_rate.SAL_post]), 2);
scop_done = ~any(isnan([Short_rate.SCOP_pre Short_rate.SCOP_post]), 2);

e_sr_sal = [Short_rate.SAL_pre Short_rate.SAL_post]; % shortrate_saline(:, presess | postsess);
% e_sr_sal = e_sr_sal(~any(isnan(e_sr_sal),2),:);
e_sr_sal = e_sr_sal(sal_done,:);
% [t_sr_sal, p_sr_sal] = ttest2(e_sr_sal(:,1), e_sr_sal(:,2));
% [p_sr_sal] = signrank(e_sr_sal(:,1), e_sr_sal(:,2));
[p_sr_sal] = ranksum(e_sr_sal(:,1), e_sr_sal(:,2));

e_sp_sal = [Short_proportion.SAL_pre Short_proportion.SAL_post]; % e_sp_sal = shortprop_saline(:, presess | postsess);
% e_sp_sal = e_sp_sal(~any(isnan(e_sp_sal),2),:);
e_sp_sal = e_sp_sal(sal_done,:);
% [t_sp_sal, p_sp_sal] = ttest2(e_sp_sal(:,1), e_sp_sal(:,2));
% [p_sp_sal] = signrank(e_sp_sal(:,1), e_sp_sal(:,2));
[p_sp_sal] = ranksum(e_sp_sal(:,1), e_sp_sal(:,2));

e_sr_scop = [Short_rate.SCOP_pre   Short_rate.SCOP_post]; % e_sr_scop = shortrate_scopo(:, presess | postsess); 
% e_sr_scop = e_sr_scop(~any(isnan(e_sr_scop),2),:);
e_sr_scop = e_sr_scop(scop_done,:);
% [t_sr_scop, p_sr_scop] = ttest2(e_sr_scop(:,1), e_sr_scop(:,2));
% [p_sr_scop] = signrank(e_sr_scop(:,1), e_sr_scop(:,2));
[p_sr_scop] = ranksum(e_sr_scop(:,1), e_sr_scop(:,2));

e_sp_scop = [Short_proportion.SCOP_pre   Short_proportion.SCOP_post]; % e_sp_scop = shortprop_scopo(:, presess | postsess); 
% e_sp_scop = e_sp_scop(~any(isnan(e_sp_scop),2),:);
e_sp_scop = e_sp_scop(scop_done,:);
% [t_sp_scop, p_sp_scop] = ttest2(e_sp_scop(:,1), e_sp_scop(:,2));
% [p_sp_scop] = signrank(e_sp_scop(:,1), e_sp_scop(:,2));
[p_sp_scop] = ranksum(e_sp_scop(:,1), e_sp_scop(:,2));


saline_sr_mean = nanmean(shortrate_saline,1);
saline_sp_mean = nanmean(shortprop_saline,1);
scopo_sr_mean  = nanmean(shortrate_scopo,1);
scopo_sp_mean  = nanmean(shortprop_scopo,1);
figure(963); clf; 
d1 = shock_sess_diff(1)-.5;
d2 = shock_sess_diff(end)+.5;
scopocolor = [1 .5 1];
scopomeancolor = [1 0 1];
controlcolor = [.85 .85 .85];
controlmeancolor = [.45 .45 .45];

n_saline = size(e_sr_sal,1); 
n_scopo = size(e_sr_scop,1); 
subplot(221); hold on
rectangle('Position', [shock_sess_diff(presess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
rectangle('Position', [shock_sess_diff(postsess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
title(sprintf('Saline- ranksum(n=%d) p=%0.3f', n_saline, p_sr_sal))
for j = 1:A_num
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    plot(shock_sess_diff, shortrate_saline(j,:), 'Color', controlcolor, 'LineStyle', linetype, 'LineWidth', 1.5); 
end
plot(shock_sess_diff, saline_sr_mean, 'Color', controlmeancolor, 'LineWidth', 2); 
ylabel('short laps / min')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 6])
% set(gca, 'XTickLabel', shock_sess_diff+1)

subplot(222); hold on
rectangle('Position', [shock_sess_diff(presess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
rectangle('Position', [shock_sess_diff(postsess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
title(sprintf('Saline- ranksum(n=%d) p=%0.3f', n_saline, p_sp_sal))
for j = 1:A_num
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    plot(shock_sess_diff, shortprop_saline(j,:), 'Color', controlcolor, 'LineStyle', linetype, 'LineWidth', 1.5); 
end
plot(shock_sess_diff, saline_sp_mean, 'Color', controlmeancolor, 'LineWidth', 2); 
ylabel('short laps / all laps')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 1])
% set(gca, 'XTickLabel', shock_sess_diff+1)
set(gca, 'YTick', [0:.2:1])

subplot(223); hold on
rectangle('Position', [shock_sess_diff(presess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
rectangle('Position', [shock_sess_diff(postsess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
title(sprintf('Scopo- ranksum(n=%d) p=%0.3f', n_scopo, p_sr_scop))
for j = 1:A_num
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    plot(shock_sess_diff, shortrate_scopo(j,:), 'Color', scopocolor, 'LineStyle', linetype, 'LineWidth', 1.5); 
end
plot(shock_sess_diff, scopo_sr_mean, 'Color', scopomeancolor, 'LineWidth', 2); 
ylabel('short laps / min')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 6])
% set(gca, 'XTickLabel', shock_sess_diff+1)

subplot(224); hold on
rectangle('Position', [shock_sess_diff(presess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
rectangle('Position', [shock_sess_diff(postsess)-.2, 0, .4, 6], 'EdgeColor', 'none', 'FaceColor', 'k')
title(sprintf('Scopo- ranksum(n=%d) p=%0.3f', n_scopo, p_sp_scop))
for j = 1:A_num
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    plot(shock_sess_diff, shortprop_scopo(j,:), 'Color', scopocolor, 'LineStyle', linetype, 'LineWidth', 1.5); 
end
plot(shock_sess_diff, scopo_sp_mean, 'Color', scopomeancolor, 'LineWidth', 2); 
ylabel('short laps / all laps')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 1])
% set(gca, 'XTickLabel', shock_sess_diff+1)
set(gca, 'YTick', [0:.2:1])
set(gcf, 'Color', 'w')

fname = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\fig_scopo_shock_paired_supp.fig';
savefig(gcf, fname)
fname = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\fig_scopo_shock_paired_supp.tiff';
temp = getframe(gcf);
imwrite(temp.cdata, fname)
%%
figure(964); clf; 
plotx_prepost = shock_sess_diff(presess|postsess);
d1 = plotx_prepost(1)-.5;
d2 = plotx_prepost(end)+.5;
scopocolor = [1 .2 1];
controlcolor = [.6 .6 .6];
plot_line_diff = .1;
marker_size = 50;
subplot(121); hold on
t1 = sprintf('Saline- ranksum(n=%d) p=%0.3f', n_saline, p_sr_sal);
t2 = sprintf('Scopo- ranksum(n=%d) p=%0.3f', n_scopo, p_sr_scop);
title(sprintf('%s\n%s', t1, t2))

% title(sprintf('Saline- t(n=%d) p=%0.3f', n_saline, p_sr_sal))
for j = 1:size(e_sr_sal,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost-plot_line_diff, e_sr_sal(j,:), marker_size/4, 'o', 'MarkerFaceColor', controlcolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost-plot_line_diff, e_sr_sal(j,:), 'Color', controlcolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sr_sal,1), marker_size, 'o', 'MarkerFaceColor', controlcolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sr_sal,1), 'Color', controlcolor/1.3, 'LineWidth', 2); 
ylabel('short laps / min')
xlabel('sessions to shock')
axis([d1 d2 0 5])
% set(gca, 'XTickLabel', shock_sess_diff+1)

subplot(122); hold on
t1 = sprintf('Saline- ranksum(n=%d) p=%0.3f', n_saline, p_sp_sal);
t2 = sprintf('Scopo- ranksum(n=%d) p=%0.3f', n_scopo, p_sp_scop);
title(sprintf('%s\n%s', t1, t2))
for j = 1:size(e_sp_sal,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost-plot_line_diff, e_sp_sal(j,:), marker_size/4, 'o', 'MarkerFaceColor', controlcolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost-plot_line_diff, e_sp_sal(j,:), 'Color', controlcolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sp_sal,1), marker_size, 'o', 'MarkerFaceColor', controlcolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sp_sal,1), 'Color', controlcolor/1.3, 'LineWidth', 2); 
ylabel('short laps / all laps')
xlabel('sessions to shock')
axis([d1 d2 0 1])
% set(gca, 'XTickLabel', shock_sess_diff+1)
set(gca, 'YTick', [0:.2:1])

subplot(121); hold on
% title(sprintf('Scopo- t(n=%d) p=%0.3f', n_scopo, p_sr_scop))
for j = 1:size(e_sr_scop,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost+plot_line_diff, e_sr_scop(j,:), marker_size/4, 'd', 'MarkerFaceColor', scopocolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost+plot_line_diff, e_sr_scop(j,:), 'Color', scopocolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sr_scop,1), marker_size, 'd', 'MarkerFaceColor', scopocolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sr_scop,1), 'Color', scopocolor/1.3, 'LineWidth', 2); 
ylabel('short laps / min')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 5])
% set(gca, 'XTickLabel', shock_sess_diff+1)

subplot(122); hold on
% title(sprintf('Scopo- t(n=%d) p=%0.3f', n_scopo, p_sp_scop))
for j = 1:size(e_sp_scop,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost+plot_line_diff, e_sp_scop(j,:), marker_size/4, 'd', 'MarkerFaceColor', scopocolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost+plot_line_diff, e_sp_scop(j,:), 'Color', scopocolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sp_scop,1), marker_size, 'd', 'MarkerFaceColor', scopocolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sp_scop,1), 'Color', scopocolor/1.3, 'LineWidth', 2); 
ylabel('short laps / all laps')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 1])
% set(gca, 'XTickLabel', shock_sess_diff+1)
set(gca, 'YTick', [0:.2:1])

set(gcf, 'Color', 'w')

fname = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\fig_pre_post_paired_supp.fig';
savefig(gcf, fname)
fname = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\fig_pre_post_paired_supp.tiff';
temp = getframe(gcf);
imwrite(temp.cdata, fname)

%%%%%%%%%
%% 
e_sr_salscop_pre = [Short_rate.SAL_pre Short_rate.SCOP_pre]; % shortrate_saline(:, presess | postsess);
e_sr_salscop_pre = e_sr_salscop_pre(all_sess_done,:);
[p_sr_salscop_pre] = ranksum(e_sr_salscop_pre(:,1), e_sr_salscop_pre(:,2));

e_sp_salscop_pre = [Short_proportion.SAL_pre Short_proportion.SCOP_pre]; % e_sp_sal = shortprop_saline(:, presess | postsess);
e_sp_salscop_pre = e_sp_salscop_pre(all_sess_done,:);
[p_sp_salscop_pre] = ranksum(e_sp_salscop_pre(:,1), e_sp_salscop_pre(:,2));

e_sr_salscop_post = [Short_rate.SAL_post  Short_rate.SCOP_post]; % e_sr_scop = shortrate_scopo(:, presess | postsess); 
e_sr_salscop_post = e_sr_salscop_post(all_sess_done,:);
[p_sr_salscop_post] = ranksum(e_sr_salscop_post(:,1), e_sr_salscop_post(:,2));

e_sp_salscop_post = [Short_proportion.SAL_post   Short_proportion.SCOP_post]; % e_sp_scop = shortprop_scopo(:, presess | postsess); 
e_sp_salscop_post = e_sp_salscop_post(all_sess_done,:);
[p_sp_salscop_post] = ranksum(e_sp_salscop_post(:,1), e_sp_salscop_post(:,2));



figure(965); clf; 
plotx_prepost = shock_sess_diff(presess|postsess);
d1 = plotx_prepost(1)-.5;
d2 = plotx_prepost(end)+.5;
scopocolor = [1 .2 1];
controlcolor = [.6 .6 .6];
plot_line_diff = .1;
marker_size = 50;
subplot(121); hold on
t1 = sprintf('SAL-SCOP pre (GREY) - ranksum(n=%d) p=%0.3f', sum(all_sess_done), p_sr_salscop_pre);
t2 = sprintf('SAL-SCOP post (MAG) - ranksum(n=%d) p=%0.3f', sum(all_sess_done), p_sr_salscop_post);
title(sprintf('%s\n%s', t1, t2))

% title(sprintf('Saline- t(n=%d) p=%0.3f', n_saline, p_sr_sal))
for j = 1:size(e_sr_salscop_pre,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost-plot_line_diff, e_sr_salscop_pre(j,:), marker_size/4, 'o', 'MarkerFaceColor', controlcolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost-plot_line_diff, e_sr_salscop_pre(j,:), 'Color', controlcolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sr_salscop_pre,1), marker_size, 'o', 'MarkerFaceColor', controlcolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sr_salscop_pre,1), 'Color', controlcolor/1.3, 'LineWidth', 2); 
ylabel('short laps / min')
xlabel('sessions to shock')
axis([d1 d2 0 5])
% set(gca, 'XTickLabel', shock_sess_diff+1)

subplot(122); hold on
t1 = sprintf('SAL-SCOP pre (GREY) - ranksum(n=%d) p=%0.3f', sum(all_sess_done), p_sp_salscop_pre);
t2 = sprintf('SAL-SCOP pre (MAG) - ranksum(n=%d) p=%0.3f', sum(all_sess_done), p_sp_salscop_post);
title(sprintf('%s\n%s', t1, t2))
for j = 1:size(e_sp_salscop_pre,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost-plot_line_diff, e_sp_salscop_pre(j,:), marker_size/4, 'o', 'MarkerFaceColor', controlcolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost-plot_line_diff, e_sp_salscop_pre(j,:), 'Color', controlcolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sp_salscop_pre,1), marker_size, 'o', 'MarkerFaceColor', controlcolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sp_salscop_pre,1), 'Color', controlcolor/1.3, 'LineWidth', 2); 
ylabel('short laps / all laps')
xlabel('sessions to shock')
axis([d1 d2 0 1])
% set(gca, 'XTickLabel', shock_sess_diff+1)
set(gca, 'YTick', [0:.2:1], 'XTick', [-1 1], 'XTickLabel', {'SAL', 'SCOP'})

subplot(121); hold on
% title(sprintf('Scopo- t(n=%d) p=%0.3f', n_scopo, p_sr_scop))
for j = 1:size(e_sr_salscop_post,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost+plot_line_diff, e_sr_salscop_post(j,:), marker_size/4, 'd', 'MarkerFaceColor', scopocolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost+plot_line_diff, e_sr_salscop_post(j,:), 'Color', scopocolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sr_salscop_post,1), marker_size, 'd', 'MarkerFaceColor', scopocolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sr_salscop_post,1), 'Color', scopocolor/1.3, 'LineWidth', 2); 
ylabel('short laps / min')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 5])
set(gca, 'YTick', [0:2:6], 'XTick', [-1 1], 'XTickLabel', {'SAL', 'SCOP'})
% set(gca, 'XTickLabel', shock_sess_diff+1)

subplot(122); hold on
% title(sprintf('Scopo- t(n=%d) p=%0.3f', n_scopo, p_sp_scop))
for j = 1:size(e_sp_salscop_post,1)
    if is_female(j)==1
        linetype = ':';
    else
        linetype = '-';
    end
    scatter(plotx_prepost+plot_line_diff, e_sp_salscop_post(j,:), marker_size/4, 'd', 'MarkerFaceColor', scopocolor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .8); 
    plot(plotx_prepost+plot_line_diff, e_sp_salscop_post(j,:), 'Color', scopocolor, 'LineStyle', linetype, 'LineWidth', 1); 
end
scatter(plotx_prepost, mean(e_sp_salscop_post,1), marker_size, 'd', 'MarkerFaceColor', scopocolor/1.3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); 
plot(plotx_prepost, mean(e_sp_salscop_post,1), 'Color', scopocolor/1.3, 'LineWidth', 2); 
ylabel('short laps / all laps')
xlabel('sessions to shock')
axis square
axis([d1 d2 0 1])
% set(gca, 'XTickLabel', shock_sess_diff+1)
set(gca, 'YTick', [0:.2:1])

set(gcf, 'Color', 'w')

% fname = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\fig_pre_post_paired_supp.fig';
% savefig(gcf, fname)
% fname = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\LinearTrack\figs\fig_pre_post_paired_supp.tiff';
% temp = getframe(gcf);
% imwrite(temp.cdata, fname)

% temp = [Short_rate.SAL_pre Short_rate.SAL_post Short_rate.SCOP_pre Short_rate.SCOP_post]; % shortrate_saline(:, presess | postsess);
% temp = [Short_proportion.SAL_pre Short_proportion.SAL_post Short_proportion.SCOP_pre Short_proportion.SCOP_post]; % shortrate_saline(:, presess | postsess);





