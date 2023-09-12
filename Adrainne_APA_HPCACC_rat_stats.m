% adrianne stats
clear 
ddir = 'D:\Adrianne\spatial properties (early vs late)\';
cd(ddir)
temp = dir('*.mat');

regions = {'ACC', 'HPC'};
vars = {'inforate', 'placefieldnumber', 'shuffle', 'withinsesscorr', 'withinsesscorr_p'};

load([ddir 'ACC_inforate.mat'])
% cols are animal, first 4 are ealy data, last 4 are late
for v = 1:length(vars)
    %%
    temp = NaN(16,2);
    label   = cell(16,2);
    isHPC   = false(16,2);
    isearly = false(16,2);
    for r = 1:2
        fn = sprintf('load(''%s%s_%s.mat'')', ddir, regions{r}, vars{v});
        eval(fn)
        if r==1 % ACC loop
            f = whos('ACC*');
            eval(sprintf('ACC_temp = %s;', f.name))
            for i = 1:4
                temp(i+8*(r-1),1) = nanmean(ACC_temp{1,i});
                temp(i+8*(r-1),2) = nanmean(ACC_temp{2,i});
                label{i+8*(r-1),1} = 'R_C_E';
                label{i+8*(r-1),2} = 'A_C_E';
                isHPC(i+8*(r-1),:) = false;
                isearly(i+8*(r-1),:) = true;
            end
            for i = 5:8
                temp(i+8*(r-1),1) = nanmean(ACC_temp{1,i});
                temp(i+8*(r-1),2) = nanmean(ACC_temp{2,i});
                label{i+8*(r-1),1} = 'R_C_L';
                label{i+8*(r-1),2} = 'A_C_L';
                isHPC(i+8*(r-1),:) = false;
                isearly(i+8*(r-1),:) = false;
            end
        elseif r==2  % HPC loop
            f = whos('HPC*');
            eval(sprintf('HPC_temp = %s;', f.name))
            for i = 1:4
                temp(i+8*(r-1),1) = nanmean(HPC_temp{1,i});
                temp(i+8*(r-1),2) = nanmean(HPC_temp{2,i});
                label{i+8*(r-1),1} = 'R_H_E';
                label{i+8*(r-1),2} = 'A_H_E';
                isHPC(i+8*(r-1),:) = true;
                isearly(i+8*(r-1),:) = true;
            end
            for i = 5:8
                temp(i+8*(r-1),1) = nanmean(HPC_temp{1,i});
                temp(i+8*(r-1),2) = nanmean(HPC_temp{2,i});
                label{i+8*(r-1),1} = 'R_H_L';
                label{i+8*(r-1),2} = 'A_H_L';
                isHPC(i+8*(r-1),:) = true;
                isearly(i+8*(r-1),:) = false;
            end
        end
    end
    fn = sprintf('%s = temp;', vars{v});
    eval(fn)
    fn = sprintf('%s_label = label;', vars{v});
    eval(fn)
    fn = sprintf('%s_isearly = isearly;', vars{v});
    eval(fn)
    
    clearvars ACC* HPC* temp isearly isHPC label
end
%%
close all

disp(['pfield anova'])
plotting = 'on';
anova1(placefieldnumber(:), placefieldnumber_label(:), plotting)
ylim([0 5])
title('N PFIELDS')

disp(['info anova'])
[p,tbl,stats] = anova1(inforate(:), inforate_label(:), plotting)
ylim([0 1])
title('INFO')

disp(['percent pcell anova'])
anova1(shuffle(:), shuffle_label(:), plotting)
ylim([0 1])
title('% PCELL')

disp(['corr anova'])
anova1(withinsesscorr(:), withinsesscorr_label(:), plotting)
ylim([0 1])
title('CORR')

disp(['corr prob anova'])
anova1(withinsesscorr_p(:), withinsesscorr_p_label(:), plotting)
ylim([0 1])
title('CORR P')
%%
e = shuffle_isearly==true;
l = shuffle_isearly==false;
difflabel = shuffle_label(l);

numfielddiff = placefieldnumber(l) - placefieldnumber(e);%
inforatediff = inforate(l) - inforate(e);%
pcelldiff = shuffle(l) - shuffle(e);%
corrdiff = withinsesscorr(l) - withinsesscorr(e);%
corrpdiff = withinsesscorr_p(l) - withinsesscorr_p(e);%


%%



plotting = 'on';
disp(['pfield DIFF anova'])
anova1(numfielddiff, difflabel, plotting)
disp(['info DIFF anova'])
anova1(inforatediff, difflabel, plotting)
disp(['percent DIFF pcell anova'])
anova1(pcelldiff, difflabel, plotting)
disp(['corr DIFF anova'])
anova1(corrdiff, difflabel, plotting)
disp(['corr DIFF prob anova'])
anova1(corrpdiff, difflabel, plotting)


%%
a = contains(difflabel, 'A');
h = contains(difflabel, 'H');

disp(['pfield DIFF ttest2'])
[hy, p] = ttest2(numfielddiff(a), numfielddiff(h))
disp(['info DIFF anova'])
[hy, p] = ttest2(inforatediff(a), inforatediff(h))
disp(['percent DIFF pcell anova'])
[hy, p] = ttest2(pcelldiff(a), pcelldiff(h))
disp(['corr DIFF anova'])
[hy, p] = ttest2(corrdiff(a), corrdiff(h))
disp(['corr DIFF prob anova'])
[hy, p] = ttest2(corrpdiff(a), corrpdiff(h))


%% ADAPTED FROM https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
close all

% data = shuffle; label = shuffle_label;
 data = placefieldnumber; label = shuffle_label;
% data = withinsesscorr; label = shuffle_label;
% data = withinsesscorr_p; label = shuffle_label;
% data = inforate; label = shuffle_label;

% frame_region_exp
F1R1E1=[data(contains(label, 'A_C_E'))]';
F2R1E1=[data(contains(label, 'R_C_E'))]';
F1R2E1=[data(contains(label, 'A_H_E'))]';
F2R2E1=[data(contains(label, 'R_H_E'))]';
F1R1E2=[data(contains(label, 'A_C_L'))]';
F2R1E2=[data(contains(label, 'R_C_L'))]';
F1R2E2=[data(contains(label, 'A_H_L'))]';
F2R2E2=[data(contains(label, 'R_H_L'))]';

Y=[F1R1E1  F2R1E1  F1R2E1  F2R2E1  F1R1E2  F2R1E2  F1R2E2  F2R2E2];
F=[1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2]; % FRAME
R=[1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2]; % REGION
E=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]; % EXPERIENCE
S=[1 2 3 4 1 2 3 4 5 6 7 8 5 6 7 8 1 2 3 4 1 2 3 4 5 6 7 8 5 6 7 8]; % SUBJECT

nesting=[0 0 0 0; ... % This line indicates that factor F is not nested in any other factor.
         0 0 0 0; ... % This line indicates that factor R is not nested in any other factor.
         0 0 0 0; ... % This line indicates that factor E is not nested in any other factor.
         1 1 1 0];    % This line indicates that S (the subject factor) is nested under A
                    % (the 1 in position 1 on the line indicates nesting under the first factor).

[p table stats]=anovan(Y,{F R E S},...
    'model',3,...
    'random', 4,...
    'nested', nesting, ...
    'varnames',{'Frame', 'Region', 'Experience', 'Subj'});























