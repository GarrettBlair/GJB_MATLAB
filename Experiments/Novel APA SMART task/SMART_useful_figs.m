function [outvars] = SMART_useful_figs(zonestruct, prop_split, time_split, all_previous, plotting, mindist)
%%
% all_previous = true;
% prop_split = 5; % proportion of runs
% time_split = 1*60; % seconds
if nargin<6
    valid_distance = zonestruct.choices.choiceDist>=0;% & zonestruct.choices.choiceDist<=(pi-mindist);
else
    valid_distance = zonestruct.choices.choiceDist>=mindist;% & zonestruct.choices.choiceDist<=(pi-mindist);
end

%             valid = zonestruct.choices.choiceDist>=0;


t = zonestruct.choices.isTime;
nt = round(t(end)/time_split);

varsIn = {'zonestruct.choices.isCorrect', 'zonestruct.choices.isLeft'};
varsOut = {'outvars.correct_', 'outvars.left_'};
outvars = [];

for vLoop = 1:length(varsIn)
    v = eval(sprintf('%s', varsIn{vLoop}));
    nV = length(v);
    runav = zeros(nV, 1);
    for ii = 1:length(v)
        runav(ii)       = sum(v(1:ii)==1)/ii;
    end
    runav_xv = t(linspace(1, nV, nV));
    % split percent laps
    if prop_split > 0
        runav_prop = zeros(prop_split, 1);
        for ii = 1:prop_split
            if all_previous == true
                valid       = linspace(0, 1, nV)' <= ii/prop_split;
            else
                v1          = linspace(0, 1, nV)' <= ii/prop_split;
                v2          = linspace(0, 1, nV)' > (ii-1)/prop_split;
                valid       = v1 & v2;
            end
            run             = v == 1;
            runav_prop(ii)  = sum(valid&run&valid_distance)/sum(valid&valid_distance);
        end
        runav_prop_xv = t(round(linspace(nV/prop_split, nV, prop_split)));
    else
        runav_prop    = NaN;
        runav_prop_xv = NaN;
    end
    % split by time
    if time_split > 0
        runav_time = zeros(nt, 1);
        for ii = 1:nt
            if all_previous == true
                valid           = t <= ii*time_split;
            else
                v1          = t <= ii*time_split;
                v2          = t > (ii-1)*time_split;
                valid       = v1 & v2;
            end
            run             = v == 1;
            runav_time(ii)  = sum(valid&run&valid_distance)/sum(valid&valid_distance);
        end
        runav_time_xv = linspace(1, nt, nt)'*time_split;
    else
        runav_time    = NaN;
        runav_time_xv = NaN;
    end
    eval(sprintf('%srunav           = runav;',          varsOut{vLoop}));
    eval(sprintf('%srunav_xv        = runav_xv;',       varsOut{vLoop}));
    eval(sprintf('%srunav_prop      = runav_prop;',     varsOut{vLoop}));
    eval(sprintf('%srunav_prop_xv   = runav_prop_xv;',  varsOut{vLoop}));
    eval(sprintf('%srunav_time      = runav_time;',     varsOut{vLoop}));
    eval(sprintf('%srunav_time_xv   = runav_time_xv;',  varsOut{vLoop}));
    
    if plotting == true
%         if vLoop == 1; figure(); end
        subplot(length(varsIn), 1, vLoop); hold on
        plot(t(linspace(1, nV, nV)), runav, 'k')
        if prop_split>0; plot(t(round(linspace(nV/prop_split, nV, prop_split))), runav_prop, 'b'); end
        if time_split>0; plot(linspace(1, nt, nt)*time_split, runav_time, 'r'); end
    end
    
end
if plotting == true
for vLoop = 1:length(varsIn)
    subplot(length(varsIn), 1, vLoop); hold on
    ylim([0 1.1])
    plot([0 t(end)], [.5 .5], 'k:')
end
end


