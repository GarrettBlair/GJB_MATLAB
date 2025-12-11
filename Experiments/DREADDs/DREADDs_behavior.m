ddir = "E:\RecordingData\GarrettBlair\DREADDs\DAT Files\";
outdir = "E:\RecordingData\GarrettBlair\DREADDs\output\";
cd('E:\RecordingData\GarrettBlair\DREADDs')
anames = {'DRD24508','DRD24509','DRD24520','DRD24521'};
dreadd_region = {'ACC','HPC','HPC','HPC'};
% dreadd_region = {'ACC','HPC','ACC','HPC'};

% ideas
% performance during CCW vs CW2 in saline/c21
% behavior stereotypy diff in exploration under drugs

n = length(anames);
day0 = datetime(2024,8,26);

MAX_TIME_VAL = NaN;
max_tr = 38;
sesstypes = {'TR' 'RET' 'HAB', 'CON'};
inj = ["NONE" "SAL" "C21" "CNO" "DCZ"];

params = [];
APA_rat_imaging_params_current;

sessdate    = string(zeros(n, max_tr +1));
anname      = string(zeros(n, max_tr +1));
region      = string(zeros(n, max_tr +1));
sessname    = string(zeros(n, max_tr +1));
sesstype    = string(zeros(n, max_tr +1));
drugtype    = string(zeros(n, max_tr +1));
drugdose    = zeros(n, max_tr +1);
entrpermin  = NaN(n, max_tr +1);
first_entr  = NaN(n, max_tr +1);
mta         = NaN(n, max_tr +1);
dist        = NaN(n, max_tr +1);
dirpref     = NaN(n, max_tr +1);
trn         = NaN(n, max_tr +1);
daynum      = NaN(n, max_tr +1);

sessnum = ones(4,1)*[0:max_tr];
sessnum_CCW = [11 14 18 21 24 27 29 32];
sessnum_SAL = [8 9 16 17 18 19 34 37]; % Saline session TR#s   [13:15; 17:19];
doseage_SAL = [1 1  1  1  1  1  1  1]; % ml/kg injected
sessnum_C21 = [10 11 12 20 21 22 23 24 25 28 29 30]; % C21 session TR#s   [10:12; 20:22];
doseage_C21 = [.5 .5 .5 .5 .5 .5  1  1  1  3  3  3]; % mg/kg dose
sessnum_CNO = [31 32 33]; % CNO session TR#s
doseage_CNO = [ 3  3  3]; % mg/kg dose
sessnum_DCZ = [35 36 38]; % DCZ session TR#s
doseage_DCZ = [.1 .5 .5]; % mg/kg dose
is_CCW = false(size(sessnum));
is_SAL = false(size(sessnum));
is_C21 = false(size(sessnum));
is_CNO = false(size(sessnum));
is_DCZ = false(size(sessnum));
druglabel = ones(size(sessnum)); % 0=none, 1=saline, 2=C21, 3=CNO, 4=DCZ
% drugdose  = zeros(size(sessnum)); % ammount / concentration injected. All injections were given at 1 ml/kg volume
for i = 1:n
    is_CCW(i,:) = ismember(sessnum(i,:), sessnum_CCW(:));
    for ii = 2:length(inj)
        eval(sprintf('is_%s(i,:) = ismember(sessnum(i,:), sessnum_%s(:));', inj(ii), inj(ii)))
        eval(sprintf('druglabel(i, is_%s(i,:)) = ii;', inj(ii)))
        eval(sprintf('drugdose(i, is_%s(i,:)) = doseage_%s;', inj(ii), inj(ii)))
    end
end
%


for a = 1:n
        fprintf('\n\t____  %s  ____\n', anames{a});
    for s = 0:max_tr
        %
        firststr = sprintf('%s%s', ddir, anames{a});
        sname = [];
        for jj = 1:length(sesstypes)
            r_fn = sprintf('%s%s_%s%d_Room.dat', ddir, anames{a}, sesstypes{jj}, s);
            a_fn = sprintf('%s%s_%s%d_Arena.dat', ddir, anames{a}, sesstypes{jj}, s);
            if isfile(r_fn) && isfile(a_fn)
                sname = [sesstypes{jj} num2str(s)];
                stype = [sesstypes{jj}];
                unique_sname = sprintf('%s_%s%d', anames{a}, sesstypes{jj}, s);
                break
            end
        end
        if ~isempty(sname)
            sessname(a, s+1) = sname;
            sesstype(a, s+1)    = stype;
            disp(unique_sname)
            [room, arena, params_sub] = behavior_DAT_tracking_eval(r_fn, a_fn, params);
            t = room.DAT_fileinfo.datetime;
            sessday = sprintf('%d_%02.0f_%02.0f', t.Year, t.Month, t.Day);
%             sesstime = sprintf('%d_%02.0f_%02.0f', t.Hour, t.Minute, t.Second);
            sessdate(a, s+1) = sessday;
            anname(a, s+1)   = anames(a);
            region(a, s+1)   = dreadd_region(a);
%             sesstime(a, s+1) = sesstime;
            daynum(a, s+1) = floor(days( room.DAT_fileinfo.datetime - day0 ));
            trn(a, s+1) = s;
            dist(a, s+1) = sum(arena.speed.*arena.dt)/100; % total meters
            [theta, rho] = cart2pol(arena.x, arena.y);
            u = unwrap(theta);
            dirpref(a, s+1) = sum(diff(u))/sum(abs(diff(u)));
            drugtype(a, s+1) = inj(druglabel(a, s+1));
            sessmins = (room.timestamps(end)/60000);
            if room.num_entrances == 0
                entrpermin(a, s+1) = room.num_entrances/sessmins;
                max_time = MAX_TIME_VAL; % room.timestamps(end)/1000 + 30;
                first_entr(a, s+1) = max_time;
                mta(a, s+1) = max_time;
            else
                entrpermin(a, s+1) = room.num_entrances/sessmins;
                first_entr(a, s+1) = room.first_entrance;
                et = room.timestamps(room.entrance_start_idx);
                m = max(abs(diff([0; et; room.timestamps(end)])))/1000;
                mta(a, s+1) = m;
            end
        else
            disp('Uknown sesstype')
            sessname(a, s+1)    = '?';
        end
    end
end

% DRD24521 did not do CON38, ended early due to stress
%%
vars = {   'Date',    'Name',   'DREADDs',    'SessName',    'Drug',     'Dose',     'EntrPerMin',  'FirstEntr', 'MaxTimeAvoided', 'DistanceRan', 'DirPreference', 'DayNum', 'TRNum', 'Type', 'isCCW'};
B = table(sessdate(:), anname(:),  region(:), sessname(:), drugtype(:), drugdose(:), entrpermin(:), first_entr(:), mta(:), dist(:), dirpref(:), daynum(:), trn(:), sesstype(:), is_CCW(:), 'VariableNames', vars);
if sum(contains(B.SessName, '?'))>0
    fprintf('\nRemoving %d rows, could not identify session.\n', sum(contains(B.SessName, '?')))
    B = B(contains(B.SessName, '?')==false, :);
end
writetable(B, 'E:\RecordingData\GarrettBlair\DREADDs\DREADDs_summary.csv')


%%
include_HAB = false;
normalize_pre = false;
figure; hold on
subsess = (is_C21 | is_SAL) & contains( sessname, 'TR');

vars = {'entrpermin', 'mta', 'dist', 'dirpref'};
an_color = jet(n*2);
an_color = an_color(n*2:-1:1,:);
% an_color = an_color(floor(n/2):floor(n/2)+n-1, :);
% an_color = [plasma(n); viridis(n)];
an_color = [an_color(2:3,:); an_color(5:6,:)];
an_offset = linspace(-.1, .1, n);
close all

drugcenter = [0 1 2 3 4];
for v = 1:length(vars)
    f1 = figure(100*v); hold on
    f2 = figure(100*v+1); hold on
for i = 1:n
    for j = 1%:2
        %%
        
        sn = sessnum_SAL(j,:)+1;
        
        if include_HAB==true; sn = [1, sn]; end
        ns = length(sn);
%         x_sal = [1:ns] + an_offset(i);
        x_sal = 0*[1:ns] + drugcenter(1) + an_offset(i);
        y_sal = eval(sprintf('%s(i, sn)', vars{v}));
        if normalize_pre == true
            if include_HAB==true
                y_sal = y_sal-y_sal(2);
            else
                y_sal = y_sal-y_sal(1);
            end
        end
        figure(f1)
        plot(x_sal, y_sal, 'k-')
        scatter(x_sal, y_sal, 'o', 'MarkerFaceColor', an_color(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .4)
        
        for jj = 2:length(inj)
%             sn = sessnum_C21(j,:)+1;
            sn = eval(sprintf('sessnum_%s(j,:)+1', inj{jj}));
            if include_HAB==true; sn = [1, sn]; end
            ns = length(sn);
            x = 0*[1:ns] + drugcenter(jj) + an_offset(i);
%             x = [1:ns] + an_offset(i);
            y_c21 = eval(sprintf('%s(i, sn)', vars{v}));
            if normalize_pre == true
                if include_HAB==true
                    y_c21 = y_c21-y_c21(2);
                else
                    y_c21 = y_c21-y_c21(1);
                end
            end
            figure(f1);
%             subplot(1,4,jj); 
            hold on
            plot(x_sal, y_sal, 'k-')
            scatter(x_sal, y_sal, 'o', 'MarkerFaceColor', an_color(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .4)
            plot(x, y_c21, 'r-')
            scatter(x, y_c21, 'o', 'MarkerFaceColor', an_color(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .4)


    %         ydiff = [y_sal(3)-y_sal(2), y_sal(2)-y_sal(1);
    %                 y_c21(3)-y_c21(2), y_c21(2)-y_c21(1)];
    %         ydrug = ydiff(1,:) - ydiff(2,:);
    %         ydrug = [y_sal(3)-y_c21(3), y_sal(2)-y_c21(2)];
            ydrug = [y_sal(2) y_c21(2)];
            figure(f2);
            x=[1,2];
%             subplot(1,4,jj); hold on
            plot(x, ydrug, 'r-')
            scatter(x, ydrug, 'o', 'MarkerFaceColor', an_color(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .4)
        end
    end
end
        figure(f1);
title(vars{v})
xlim([.5 ns+.5])
        figure(f2);
title(vars{v})
xlim([.5 2+.5])
end
% axis = 































