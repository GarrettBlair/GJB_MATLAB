%%
clear
outdir = "E:\SuaKim\SMART task data\behavior output\";
% % anames = {  'SUA34991', 'SUA34992', 'SUA34993', 'SUA34994', 'SUA34995', 'SUA34996',...
% %     'SUA35008', 'SUA35009', 'SUA35010', 'SUA35011', 'SUA35012', 'SUA35013'};
% % load(metricsFileName)
% 
% metersran_per_min = ((distRan./(sessmins*60))*100);
% 

matFile = sprintf('%s%s', outdir, '\SUA_behavior_metrics.mat');
% matFile        = 'SUA_behavior_metrics (1).mat';
% Sessions to plot (by absolute session number in trainNum)
wantedSessNums = [0 1 2 9 10 11 14 15 16 17 18 21 22 23 24 25];
wantedLbls     = {'HAB0','IL1','IL2', 'IL9', 'IL10', 'IL11','IL14','IL15','RET16','IL17','IL18', 'IL21', 'IL22'...
    'CON23', 'CON24', 'CON25'};

wantedSessNums = [0 16 27];
wantedLbls     = {'HAB0' 'RET16' 'RET27'};

saveFigures    = false;           % true to save PNGs
outDir         = 'session_plots';
maxCols        = 3;               % columns per grid
% If you want to include only certain metrics, list them here:
includeOnly    = {};              % e.g., {'entrpermin','arena_displacement_max'}
S = load(matFile);
assert(isfield(S,'trainNum'), 'Variable "trainNum" not found in MAT file.');
trainNum = S.trainNum;                  % animals x sessionSlots
assert(ndims(trainNum)==2, '"trainNum" must be 2-D (animals x sessions).');
[numAnimals, numSessSlots] = size(trainNum);
if isfield(S,'anames')
    an = S.anames;
    if iscell(an), an = string(an);
    elseif ischar(an), an = string(cellstr(an));
    end
    an = an(:);
    if numel(an) == numAnimals
        animalLabels = an;
    else
        animalLabels = "Rat " + string(1:numAnimals);
    end
else
    animalLabels = "Rat " + string(1:numAnimals);
end
% Define which animals should be red
redAnimals = ["SUA34996","SUA35009","SUA35010","SUA35011"];%,"SUA35013"];
% Build a color map: red if label in list, black otherwise
colorMap = zeros(numAnimals,3); % preallocate RGB
for a = 1:numAnimals
    if any(animalLabels(a) == redAnimals)
        colorMap(a,:) = [1 0 0];   % red
    else
        colorMap(a,:) = [0 0 0];   % black
    end
end
% For each wanted session number, for each animal, find the column that matches.
% We'll build a matrix colIdx (animals x numWanted) with NaN where missing.
numWanted = numel(wantedSessNums);
colIdx = nan(numAnimals, numWanted);
for j = 1:numWanted
    sNum = wantedSessNums(j);
    for a = 1:numAnimals
        c = find(trainNum(a,:) == sNum, 1, 'first');
        if ~isempty(c), colIdx(a,j) = c; end
    end
end
% Drop any wanted sessions that literally no animal has:
keepMask = any(~isnan(colIdx), 1);
if ~any(keepMask)
    error('None of the requested session numbers [%s] were found in trainNum.', num2str(wantedSessNums));
end
colIdx      = colIdx(:, keepMask);
wantedSess  = wantedSessNums(keepMask);
wantedLbls  = wantedLbls(keepMask);
numWanted   = numel(wantedSess);
% -------- derive: distRan per minute (distRan ./ sessmins) --------
if isfield(S,'distRan') && isfield(S,'sessmins')
    D = S.distRan;          % distance (e.g., cm)
    T = S.sessmins;         % session time in minutes
    % Ensure T is numeric
    if isstring(T) || ischar(T) || iscellstr(T)
        T = str2double(T);
    end
    % Guard: invalid/zero times -> NaN to avoid divide-by-zero
    T(~isfinite(T) | T<=0) = NaN;
    % Make T broadcast-compatible with D (usually [animals x sessions ...])
    ndD = ndims(D);
    szT = size(T);
    if numel(szT) < ndD
        T = reshape(T, [szT, ones(1, ndD - numel(szT))]);  % add trailing singletons
    end
    % Try to align to [animals x sessions ...]
    [nA, nS] = deal(numAnimals, numSessSlots);
    if ~isequal(size(T), size(D))
        if isequal(size(T), [nA nS])
            % OK
        elseif isequal(size(T), [nS nA])
            T = permute(T, [2 1]);           % swap -> [animals x sessions]
        elseif isequal(size(T), [1 nS]) || isequal(size(T), [nS 1])
            T = reshape(T, [1 nS]);          % session-only -> expand across animals
            T = repmat(T, [nA 1]);
        end
        % Re-add trailing singletons if D has extra dims
        szT = size(T);
        if numel(szT) < ndD
            T = reshape(T, [szT, ones(1, ndD - numel(szT))]);
        end
    end
    % Compute ratio (units: e.g., cm/min if D is cm)
    R = double(D) ./ double(T);
    % Helpful debug print (you can comment out later)
    fprintf('distRan size: %s | sessmins size: %s | ratio finite count: %d\n', ...
        mat2str(size(D)), mat2str(size(T)), nnz(isfinite(R)));
    S.distRan_per_sessmins = R;   % <— new field name is explicit about minutes
end
%includeOnly = {'distRan_per_sessmins'};
allNames = string(fieldnames(S));
exclude = ["trainNum","dayNum","trainNum_sessType","sesstypes","sessmins", ...
           "anames","is_expt","is_fem","dreadd_region","cm","rats_cmap", ...
           "distbins","distRan","c","e","clrs","cam_info_filename", ...
           "sleap_dir","sleap_modelname","sleap_vid_timestamp_dir", ...
           "MAX_TIME_VAL","None", 'swap_performance_all', 'arena_displacement_max'];
if ~isempty(includeOnly)
    candNames = string(includeOnly(:)');
else
    candNames = setdiff(allNames, exclude);
end
plottable = strings(0,1);
sessDimFor = containers.Map('KeyType','char','ValueType','double');
aniDimFor  = containers.Map('KeyType','char','ValueType','double');
for k = 1:numel(candNames)
    fn = char(candNames(k));
    if ~isfield(S, fn), continue; end
    A = S.(fn);
    if ~isnumeric(A), continue; end
    sz = size(A);
    sessDim = find(sz == numSessSlots, 1);   % needs a sessions dimension
    if isempty(sessDim), continue; end
    plottable(end+1) = string(fn); %#ok<AGROW>
    sessDimFor(fn) = sessDim;
    ad = find(sz == numAnimals, 1);
    if isempty(ad), ad = 0; end
    aniDimFor(fn) = ad;
end
assert(~isempty(plottable), 'No numeric variables with a session dimension (%d) found to plot.', numSessSlots);
nPlots  = numel(plottable);
halfIdx = ceil(nPlots/2);
groups  = { plottable(1:halfIdx), plottable(halfIdx+1:end) };
if saveFigures && ~exist(outDir,'dir'), mkdir(outDir); end
for g = 1:2
    vars = groups{g};
    if isempty(vars), continue; end
    nCols = min(maxCols, numel(vars));
    nRows = ceil(numel(vars)/nCols);
    % +1 row for the key tile at bottom
    f = figure(100+g);
    set(f, 'Color','w','Name',sprintf('Metrics Part %d', g));
    clf
    try
        tlo = tiledlayout(f, nRows+1, nCols, 'TileSpacing','compact', 'Padding','compact');
    catch
        tlo=[];
    end
    for k = 1:numel(vars)
        fn = char(vars(k));
        A  = S.(fn);
        sz = size(A);
        sessDim = sessDimFor(fn);
        aniDim  = aniDimFor(fn);
        dims    = 1:ndims(A);
        % Ensure an animal dimension exists; if not, insert a singleton at front
        if aniDim == 0
            A = reshape(A, [1, sz]);          % add leading animals dim
            aniDim  = 1;
            sessDim = sessDim + 1;
            sz = size(A);
            dims = 1:ndims(A);
        end
        % Permute to [animals, sessions, rest...]
        permOrder = [aniDim, sessDim, setdiff(dims, [aniDim, sessDim], 'stable')];
        A = permute(A, permOrder);
        sz = size(A); % now sz(1)=animals, sz(2)=sessions
        % Pull values for requested sessions by per-animal column index
        vals = nan(numAnimals, numWanted);
        for j = 1:numWanted
            for a = 1:numAnimals
                c = colIdx(a,j);
                if ~isnan(c)
                    v = A(a, c, :);
                    if ndims(A) > 2
                        v = mean(v, 3:ndims(A), 'omitnan');  % average extra dims
                    end
                    vals(a,j) = v;
                end
            end
        end
        % If everything is NaN for this metric (for selected sessions), skip it
        if all(isnan(vals), 'all')
            continue;
        end
        % Plot tile
        if isempty(tlo)
            ax = subplot_tight(nRows, nCols, k, [.1 .1]);
        else
            ax = nexttile(tlo);
        end
        hold(ax,'on'); grid(ax,'on');
        for a = 1:numAnimals
            plot(ax, 1:numWanted, vals(a,:), '-', ...
                'Color', colorMap(a,:), 'LineWidth', 1.2, 'MarkerSize', 5);
        end
        hold(ax,'off');
        xticks(ax, 1:numWanted);
        xticklabels(ax, wantedLbls(1:numWanted));
        xlabel(ax, 'Session');
        ylabel(ax, fn, 'Interpreter','none');
        title(ax, fn, 'Interpreter','none');
    end
    % optional save
    if saveFigures
        if ~exist(outDir, 'dir'), mkdir(outDir); end
        outPath = fullfile(outDir, sprintf('metrics_part%d.png', g));
        exportgraphics(f, outPath, 'Resolution', 200);
        fprintf('Saved %s\n', outPath);
    end
end
fprintf('Plotted %d metrics split across %d figure(s).\n', nPlots, 1 + ~isempty(groups{2}));