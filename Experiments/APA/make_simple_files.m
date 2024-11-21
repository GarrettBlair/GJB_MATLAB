clear
% APA_PKCZmouse_imaging_params_current;
fprintf('%s\n', datetime())
%
% animals = {'HPCACC24500', 'HPCACC24502', 'HPCACC24503', 'Acc20832', 'Acc19947', 'Hipp18240'};% animals = {'Acc20832', 'Acc19947'};
animals = {'Acc20832'};% animals = {'Acc20832', 'Acc19947'};
experiment_folder{1} = 'D:/GarrettBlair/APA/';
numAnimals = length(animals);

% DAT_Dir             = sprintf('%sDAT_files/', experiment_folder);

%%
resave_simple   = false;
for animalLoop = 1:numAnimals
    for expLoop = 1:length(experiment_folder)
        animal_name = animals{animalLoop};
        processed_dir = sprintf('%s%s/processed_files/', experiment_folder{expLoop}, animal_name);
        simple_dir = sprintf('%s%s/simple_files/', experiment_folder{expLoop}, animal_name);
        if ~isfolder(simple_dir)
            mkdir(simple_dir)
        end
        fprintf('\n\nREPROCESS files for %s in folder  \n%s...\n', animal_name, processed_dir)
        %             temp = dir([processed_dir '*miniscope*']);
        temp = dir([processed_dir '*@placecells*']);
        temp2 = dir([simple_dir '*@placecells*']);
        nfiles = length(temp);
        for sessionLoop = 1:nfiles
            fname = temp(sessionLoop).name;
            processedFile   = [temp(sessionLoop).folder '/' fname];
            simpleFile      = [simple_dir '/' fname];
            if ~isfile(simpleFile) % || (isfile(simpleFile) && resave_simple==true)
                %%
                fprintf('~~~PROCESSED file rerunning for: %s...\n', simpleFile)
                % first this
%                 simple_file_subfcn(processedFile, simpleFile)
                % then this
%                 update_simple_file_subfcn(simpleFile)
            end
        end
    end
end



function simple_file_subfcn(filename, newname)
%%
% ind1 = strfind(filename, 'processed_files')+length('processed_files') + 1;
% ind2 = strfind(filename, '.mat')+length('.mat') - 1;
% newname = ['D:\GarrettBlair\APA\HPCACC24500\simple_files\' filename(ind1:ind2)];
% load(filename, 'ms')
load(filename, 'ms')

vars2save = {'parentDir' 'filename' 'time_ms' 'dt_ms' 'roomx' 'roomy' 'arenax' 'arenay'...
    'spks' 'craw' 'roomEntranceTimes' 'roomShockTimes' 'arenaEntranceTimes' 'arenaShockTimes'...
    'thisregion' 'sessDate' 'trialname' 'trial_type' 'trial_num'...
    'shock_zone_center' 'shock_zone_size', ...
    'avoid_ind' 'is_entrance' 'shock_zone_dist'...
    'shk_ind' 'entr_ind' 'shock_approach_ind' 'distance_entrance_size'};

parentDir = ms.parentDir;
time_ms = single(ms.timestamps);
dt_ms   = single(ms.dt_corrected/.1000);
roomx   = single(ms.room.x);
roomy   = single(ms.room.y);
arenax  = single(ms.arena.x);
arenay  = single(ms.arena.y);
spks    = single(ms.spks);
craw    = single(normalize_rows(ms.neuron.C+ms.neuron.YrA));


roomEntranceTimes    = single(ms.room.entranceTimes);
roomShockTimes       = single(ms.room.shockTimes);
try
arenaEntranceTimes   = single(ms.arena.entranceTimes);
arenaShockTimes      = single(ms.arena.shockTimes);
catch
arenaEntranceTimes   = [];
arenaShockTimes      = [];
end
% Get some session info
[sessDate, trialname, trial_type, trial_num] = recdata_from_parentdir(ms.parentDir);
if ~isempty(strfind(ms.parentDir, 'Hipp18240/2022_09_08/16_09_56_HAB/'))
    trial_num = 0;
end
if ~isfield(ms, 'cameraName')
    ms.cameraName = 'MiniLFOV';
end
if contains(ms.cameraName, 'ACC_miniscope')
    thisregion='ACC';
elseif contains(ms.cameraName, 'HPC_miniscope')
    thisregion='HPC';
elseif contains(ms.cameraName, 'MiniLFOV') && contains(ms.parentDir, 'Acc')
    thisregion='ACC';
elseif contains(ms.cameraName, 'MiniLFOV') && contains(ms.parentDir, 'Hipp')
    thisregion='HPC';
elseif contains(ms.cameraName, 'MiniLFOV') && contains(ms.parentDir, 'HPC')
    thisregion='HPC';
end

% Calculate shock zone distances
% if strcmp(trial_type, 'TR')
if strcmp(trial_type, 'CON')
    shock_zone_center = 3*pi/2; % typical room shock configuration
    shock_zone_size = pi/6; % size in rad from center to edge
else
    shock_zone_center = pi/2; % typical room shock configuration
    shock_zone_size = pi/6; % size in rad from center to edge
end
distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure

% [shock_approach_ind, is_entrance, ~, ~, shock_zone_dist] = shock_zone_dist_peakfinding(ms.room, ms.timestamps./1000,...
%     shock_zone_center, shock_zone_size, distance_entrance_size, true);
[avoid_ind, is_entrance, ~, ~, shock_zone_dist, shk_ind, entr_ind, shock_approach_ind] = shock_zone_dist_peakfinding(str, time_ms./1000,...
    shock_zone_center, shock_zone_size, distance_entrance_size, false);
% [avoid_ind, is_entrance, ~, ~, shock_zone_dist, shk_ind, entr_ind, shock_approach_ind] = shock_zone_dist_peakfinding(ms.room, time_ms./1000,...
%     shock_zone_center, shock_zone_size, distance_entrance_size, false);
trial_num = single(trial_num);
shock_zone_size = single(shock_zone_size);
shock_zone_center = single(shock_zone_center);
% distance_entrance_size = single(distance_entrance_size);
% shock_zone_dist = single(shock_zone_dist);
% shock_approach_ind = single(shock_approach_ind);
sessDate = char(sessDate);


avoid_ind = single(avoid_ind);
is_entrance = logical(is_entrance);
shock_zone_dist = single(shock_zone_dist);
shk_ind = single(shk_ind);
entr_ind = single(entr_ind);
shock_approach_ind = single(shock_approach_ind);
distance_entrance_size = single(distance_entrance_size);

save(newname, vars2save{:})
end
function update_simple_file_subfcn(filename)
%%

vars2load={'shock_zone_center' 'shock_zone_size' 'time_ms'...
    'roomx' 'roomy' 'roomShockTimes' 'roomEntranceTimes'};
vars2append = {'avoid_ind' 'is_entrance' 'shock_zone_dist'...
    'shk_ind' 'entr_ind' 'shock_approach_ind' 'distance_entrance_size'};

load(filename, vars2load{:});
str = [];
str.x = roomx; 
str.y = roomy; 
str.shockTimes = roomShockTimes; 
str.entranceTimes = roomEntranceTimes;

distance_entrance_size = pi/2; % distance (between 0 to pi) to look at approaches to categorize escape vs failure

[avoid_ind, is_entrance, ~, ~, shock_zone_dist, shk_ind, entr_ind, shock_approach_ind] = shock_zone_dist_peakfinding(str, time_ms./1000,...
    shock_zone_center, shock_zone_size, distance_entrance_size, false);

avoid_ind = single(avoid_ind);
is_entrance = logical(is_entrance);
shock_zone_dist = single(shock_zone_dist);
shk_ind = single(shk_ind);
entr_ind = single(entr_ind);
shock_approach_ind = single(shock_approach_ind);
distance_entrance_size = single(distance_entrance_size);

save(filename, vars2append{:}, '-append')
end