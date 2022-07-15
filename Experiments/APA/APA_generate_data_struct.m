function [ms, behav, params] = APA_generate_data_struct(recording_dir, params)
%% Read in the extracted position from APA_extract_pos_batch.ipynb
% recording_dir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\2022_06_10\18_44_02';
% recording_dir = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/Hipp16942/2022_06_18/16_59_59';
global params_sub
params_sub = params;

fname_position = sprintf('%sexperiment/behav_position_data.csv', recording_dir);
fname_params = sprintf('%sexperiment/behav_ext_params.json', recording_dir);
caimanFilename = sprintf('%sMiniLFOV/caiman_cnmfe_out.mat', recording_dir);
msTSFile = sprintf('%sMiniLFOV/timeStamps.csv', recording_dir);
msCropFile = sprintf('%sMiniLFOV/Crop_params.mat', recording_dir);
msOriFile = sprintf('%sMiniLFOV/headOrientation.csv', recording_dir);
% Read in the coresponding files

crop_params = load(msCropFile);
params_sub.crop_params  = crop_params;
% Generate the miniscope structure
if exist(msTSFile, 'file')==2
    [ms] = make_ms_struct(recording_dir, msTSFile);
else
    error('No timestamp.csv file in directory:\n\t %s', recording_dir);
end
% Get the orientation data from BNO
if exist(msOriFile, 'file')==2
    [ms] = make_ori_struct(ms, msOriFile);
else
    warning('No BNO data:\n\t %s', msOriFile);
    ms.ori = [];
end
% Get the position data if extracted
if exist(fname_position, 'file')==2
    [ms, behav, behav_params] = make_behav_data(ms, fname_position, fname_params);
else
    warning('No extracted behavior data:\n\t %s', fname_position);
    behav = [];
    behav_params = [];
end

ms.caimanFilename = caimanFilename;
params_sub.behav_params = behav_params;
params = params_sub;

end
%%
function [ms] = make_ms_struct(recording_dir, msTSFile)
global params_sub
TS_data = readtable(msTSFile);
ms = [];

recording_dir(strfind(recording_dir, '\')) = '/';

ms.parentDir = recording_dir;
ms.spatialDownsample = params_sub.crop_params.spatialDownSample;
ms.temporalDownsample = params_sub.crop_params.temporalDownSample;
ms.fileName = params_sub.crop_params.tiffStackout;
ms.fileName(strfind(ms.fileName, '/')) = '\';
[ms.height, ms.width] = size(imread(string(ms.fileName),1));
ms.frameNum = TS_data.FrameNumber(1:ms.temporalDownsample:end);
ms.timestamps = TS_data.TimeStamp_ms_(1:ms.temporalDownsample:end);
tiff_numFrames = size(imfinfo(ms.fileName),1);
ms_dt = [median(diff(ms.timestamps)); diff([ms.timestamps])]/1000;
ms.dt = ms_dt;
if params_sub.correct_dt
    % sometimes the camera will disconnect and reconnect, with large jumps
    % in the timestamp file
    bad_dt_thresh = mean(ms_dt)+10*std(ms_dt);
    bad_vals = ms_dt >= bad_dt_thresh;
    if sum(bad_vals)/length(ms_dt)>.01
       warning('Many bad values found in timestamp dt, should check! %d%% bad', ceil(100*sum(bad_vals)/length(ms_dt))) 
    end
    ms_dt(ms_dt >= bad_dt_thresh) = median(ms_dt);
end
ms.dt_corrected = ms_dt;

[matches_bad, ~] = APA_troublesome_sessions({ms.parentDir});
frameMismatch = length(ms.frameNum) ~= tiff_numFrames;
if frameMismatch | any(matches_bad)
    warning('Diff found between imestamp file and tiff file!')
    % Keeping track of known files where this happens
    switch ms.parentDir
        case 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/Hipp16942/2022_06_10/18_25_10/'
            fprintf('\t%s\n', ms.parentDir)
            [ms] = cutoff_session(ms, tiff_numFrames);
        case 'C:/Users/gjb326/Desktop/RecordingData/AlejandroGrau/TestMouse1/2022_07_05/17_09_41/'
            % do nothing, first 3 frames removed            
        otherwise
            error('Unkown issue, update ''APA_troublesome_sessions.mat''')
    end
end
end
%%
function [ms, behav, behav_params] = make_behav_data(ms, fname_position, fname_params)
global params_sub

[behav, behav_params] = read_APA_csv(fname_position, fname_params);

nanind = behav.x==0 & behav.y == 0;
%
x = behav.x;
y = behav.y;
x(nanind) = NaN;
y(nanind) = NaN;
t = behav.timestamps;
if length(unique(t)) ~= length(t)
    ind = find(diff(t)==0)+1;
    t(ind) = t(ind)+1;
end
ts = t./1000; % convert from ms to seconds
if params_sub.nan_interp && any(~nanind)
    nanind = (isnan(x) & isnan(y));
    xn = interp1(ts(~nanind), x(~nanind), ts(nanind), 'linear');
    yn = interp1(ts(~nanind), y(~nanind), ts(nanind), 'linear');
    x(nanind) = xn; y(nanind) = yn;
end
behav.wasnan = nanind;

x = params_sub.arena_radius*(x-behav_params.center(1))./behav_params.radius;
y = params_sub.arena_radius*(y-behav_params.center(2))./behav_params.radius;
% x = movmedian(x, round(params_sub.behav_fps/2));
% y = movmedian(y, round(params_sub.behav_fps/2));
ksize = round(behav_params.behav_fps*params_sub.behav_smoothing_interval);
kern = ones(ksize, 1); kern = kern./sum(kern(:));
x = conv(x, kern, 'same');
y = conv(y, kern, 'same');

ax = behav.arena_x;
ay = behav.arena_y;
nanind = behav.arena_x==0 & behav.arena_y == 0;
ax(nanind) = NaN;
ay(nanind) = NaN;
ax = params_sub.arena_radius*(ax-behav_params.center(1))./behav_params.radius;
ay = params_sub.arena_radius*(ay-behav_params.center(2))./behav_params.radius;

% distsRAT = sqrt((x -  0).^2  + (y - 0).^2);
distsLED = sqrt((ax - 0).^2 + (ay - 0).^2);
nanind = distsLED<=params_sub.arena_radius*1.2 & distsLED<=params_sub.arena_radius*-1.2; % 100 cm ring threshold for excluding arena LED
ax(nanind) = NaN;
ay(nanind) = NaN;
%
if params_sub.nan_interp && any(~nanind)
    nanind = (isnan(ax) & isnan(ay));
    xn = interp1(ts(~nanind), ax(~nanind), ts(nanind), 'spline');
    yn = interp1(ts(~nanind), ay(~nanind), ts(nanind), 'spline');
    ax(nanind) = xn; ay(nanind) = yn;
end
% ax = conv(ax, kern, 'same');
% ay = conv(ay, kern, 'same');
% ax = movmedian(ax, ksize);
% ay = movmedian(ay, ksize);

[theta1,rho1] = cart2pol(x,y);
[theta2,rho2] = cart2pol(ax,ay);
rho2 = movmedian(rho2, ksize);
% ay = movmedian(ay, ksize);
rhodiff = rho2-median(rho2); % small fluctation of maze position from being offcenter
[xx, yy] = pol2cart(theta1-theta2, rho1-rhodiff);
[x, y] = pol2cart(theta1, rho1-rhodiff);
[axx, ayy] = pol2cart(theta2-theta2, rho2-rhodiff);
[ax, ay] = pol2cart(theta2, rho2-rhodiff);

% arena.speed
% arena.x
% arena.y
% arena.dt

behav_dt = [median(diff(behav.timestamps)); diff([behav.timestamps])]/1000;
behav.dt = behav_dt;
if params_sub.correct_dt
    % sometimes the camera will disconnect and reconnect, with large jumps
    % in the timestamp file
    bad_dt_thresh = mean(behav_dt)+10*std(behav_dt);
    bad_vals = behav_dt >= bad_dt_thresh;
    if sum(bad_vals)/length(behav_dt)>.01
       warning('Many bad values found in timestamp dt, should check!') 
    end
    behav_dt(bad_vals) = median(behav_dt);
end
behav.dt_corrected = behav_dt;



behav.spd_roomframe = sqrt(diff([x(1); x]).^2 + diff([y(1); y]).^2)./behav_dt;
behav.spd_arenaframe = sqrt(diff([xx(1); xx]).^2 + diff([yy(1); yy]).^2)./behav_dt;
% adt = [median(diff(behav.timestamps)); diff([behav.timestamps])]/1000;
% axspd = sqrt(diff([ax(1); ax]).^2 + diff([ay(1); ay]).^2)./adt;
% axxspd = sqrt(diff([axx(1); axx]).^2 + diff([ayy(1); ayy]).^2)./adt;
ms.room.x = interp1(t, x, ms.timestamps, 'linear');
ms.room.y = interp1(t, y, ms.timestamps, 'linear');
ms.room.led_x = interp1(t, ax, ms.timestamps, 'linear');
ms.room.led_y = interp1(t, ay, ms.timestamps, 'linear');
ms.arena.x = interp1(t, xx, ms.timestamps, 'linear');
ms.arena.y = interp1(t, yy, ms.timestamps, 'linear');
ms.arena.led_x = interp1(t, axx, ms.timestamps, 'linear');
ms.arena.led_y = interp1(t, ayy, ms.timestamps, 'linear');

dt = ms.dt_corrected;
ms.room.speed = sqrt(diff([ms.room.x(1); ms.room.x]).^2 + diff([ms.room.y(1); ms.room.y]).^2)./dt;
ms.arena.speed = sqrt(diff([ms.arena.x(1); ms.arena.x]).^2 + diff([ms.arena.y(1); ms.arena.y]).^2)./dt;
ms.room.speed_smooth = conv(ms.room.speed,   kern, 'same');
ms.arena.speed_smooth = conv(ms.arena.speed, kern, 'same');
is_moving = ms.arena.speed_smooth>params_sub.min_spd_thresh;

[~, speed_epochs] = get_speed_epochs(ms.arena.speed_smooth, params_sub);

ms.speed_epochs = speed_epochs;
ms.is_moving = is_moving;


if params_sub.plotting
    figure(3); clf;
    subplot(2,2,1)
    hold on;
    scatter3(behav.arena_x, behav.arena_y, ts, 1, 'k.')
    scatter3(behav.x, behav.y, ts, 1, 'b.')
    title('Extracted')
    set(gca, 'View', [-60 60])
    axis square
    
    subplot(2,2,2)
    hold on;
    plot3(ax, ay, ts, 'k-')
    plot3(x, y, ts, 'b-')
    title('Interp & Scaled')
    set(gca, 'View', [-60 60])
    axis square
    
    subplot(2,2,3)
    hold on;
    plot3(x, y, ts, 'b-')
    plot3(xx, yy, ts, 'r-')
    title('Room (b) vs. Arena (r)')
    set(gca, 'View', [-60 60])
    axis square
    
    subplot(4,2,6)
    % hold on;
    plot(ts, behav.spd_roomframe, 'b-')
    title('Room speed (b)')
    axis tight
    subplot(4,2,8)
    plot(ts, behav.spd_arenaframe, 'r-')
    title('Arena speed (r)')
    % set(gca, 'View', [-60 60])
    axis tight
end
end
%%
function [ms] = make_ori_struct(ms, msOriFile)

%%
global params_sub
ORI_Data = readtable(msOriFile);
shared_ts = ismember(ORI_Data.TimeStamp_ms_, ms.timestamps);

ori_ts = ORI_Data.TimeStamp_ms_(shared_ts);
q = [ORI_Data.qw, ORI_Data.qx, ORI_Data.qy, ORI_Data.qx];
q = q(shared_ts, :);
[~, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
ms.ori = ORI_Data(shared_ts, :);
% [ori_str] = LFOV_ori_extract(q, [0:30:360]);

ms.ori.time = ori_ts; 
ms.ori.roll = roll; 
ms.ori.pitch = pitch; 
ms.ori.yaw = yaw;

if params_sub.plotting
    figure;
    hold on;
    plot(roll, 'r.');
    plot(pitch, 'm.');
    plot(yaw, 'y.')
end
end