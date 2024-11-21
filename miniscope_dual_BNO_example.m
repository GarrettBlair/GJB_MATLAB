
% Generate the miniscope structure
miniscopeName = 'HPC_miniscope1';
recording_dir = 'F:\GarrettBlair\APA\HPCACC24504\2024_01_17\14_10_02\';
msTSFile = sprintf('%s%s/timeStamps.csv', recording_dir, miniscopeName);
msOriFile = sprintf('%s%s/headOrientation.csv', recording_dir, miniscopeName);

TS_data = readtable(msTSFile);
if strcmp(TS_data.Properties.VariableNames{1}, 'Var1') % timestamps file was resaved after removing bad frames and col names were not saved
    TS_data.Properties.VariableNames = {'FrameNumber', 'TimeStamp_ms_', 'Buffer'};
end
hpc = [];
hpc.timestamps = TS_data.TimeStamp_ms_;

% [ms] = make_ms_struct(recording_dir, msTSFile, miniscopeName);
    
[hpc] = make_ori_struct(hpc, msOriFile);



miniscopeName = 'ACC_miniscope2';
recording_dir = 'F:\GarrettBlair\APA\HPCACC24504\2024_01_17\14_10_02\';
msTSFile = sprintf('%s%s/timeStamps.csv', recording_dir, miniscopeName);
msOriFile = sprintf('%s%s/headOrientation.csv', recording_dir, miniscopeName);

TS_data = readtable(msTSFile);
if strcmp(TS_data.Properties.VariableNames{1}, 'Var1') % timestamps file was resaved after removing bad frames and col names were not saved
    TS_data.Properties.VariableNames = {'FrameNumber', 'TimeStamp_ms_', 'Buffer'};
end

acc = [];
acc.timestamps = TS_data.TimeStamp_ms_;

% [ms] = make_ms_struct(recording_dir, msTSFile, miniscopeName);
    
[acc] = make_ori_struct(acc, msOriFile);
%%
valid = ms1_time<=110000;
ms1_time = hpc.timestamps(valid);
% ms1_ori = [hpc.ori.roll(valid) hpc.ori.pitch(valid) hpc.ori.yaw(valid)];
% because of the mirror mounting, need to be inverted except yaw
ms1_ori = [-1*hpc.ori.roll(valid) -1*hpc.ori.pitch(valid) hpc.ori.yaw(valid)];
ms2_time = acc.timestamps(valid);
ms2_ori = [acc.ori.roll(valid) acc.ori.pitch(valid) acc.ori.yaw(valid)];

for i = 1:3
    ms2_ori(:,i) = interp1(ms2_time, ms2_ori(:,i), ms1_time, 'linear');
end
ms1_time = ms1_time./1000;
ms2_time = ms1_time;%acc.timestamps(valid);


figure(1);
clf
subplot(2,3,1); 
hold on
plot(ms1_time, ms1_ori(:,1))
plot(ms2_time, ms2_ori(:,1))
axis tight; ylim([-pi-.5 pi+.5])
title('Roll')
legend({'Cam 1', 'Cam 2'})
ylabel('radians')

subplot(2,3,4); 
plot(ms2_time, unwrap(ms2_ori(:,1) - ms1_ori(:,1)))
axis tight; ylim([-pi-.5 pi+.5])
xlabel('Seconds')
ylabel('unwrapped difference')

subplot(2,3,2); 
hold on
plot(ms1_time, ms1_ori(:,2))
plot(ms2_time, ms2_ori(:,2))
axis tight; ylim([-pi-.5 pi+.5])
title('Pitch')
ylabel('radians')

subplot(2,3,5); 
plot(ms2_time, unwrap(ms2_ori(:,2) - ms1_ori(:,2)))
axis tight; ylim([-pi-.5 pi+.5])
xlabel('Seconds')
ylabel('unwrapped difference')

subplot(2,3,3); 
hold on
% plot(ms1_time, unwrap(ms1_ori(:,3)))
% plot(ms2_time, unwrap(ms2_ori(:,3)))
plot(ms1_time, (ms1_ori(:,3)))
plot(ms2_time, (ms2_ori(:,3)))
axis tight; ylim([-pi-.5 pi+.5])
title('Yaw')
ylabel('radians')

subplot(2,3,6); 
plot(ms2_time, unwrap(ms2_ori(:,3) - ms1_ori(:,3)))
axis tight; ylim([-pi-.5 pi+.5])
xlabel('Seconds')
ylabel('unwrapped difference')
%%


function [ms] = make_ms_struct(recording_dir, msTSFile, cameraName)
global params_sub
TS_data = readtable(msTSFile);
if strcmp(TS_data.Properties.VariableNames{1}, 'Var1') % timestamps file was resaved after removing bad frames and col names were not saved
    TS_data.Properties.VariableNames = {'FrameNumber', 'TimeStamp_ms_', 'Buffer'};
end
ms = [];

recording_dir(strfind(recording_dir, '\')) = '/';

ms.parentDir = recording_dir;
ms.spatialDownsample = params_sub.crop_params.spatialDownSample;
ms.temporalDownsample = params_sub.crop_params.temporalDownSample;
% ms.fileName = params_sub.crop_params.tiffStackout;
ms.fileName = [ms.parentDir cameraName '/msCam_MC.tiff']; % params_sub.crop_params.tiffStackout;
if ~isfile(ms.fileName)&& isfile([ms.parentDir cameraName '/msCam.tiff'])
    % I prob changed the directory name since cropping
    f = [ms.parentDir cameraName '/msCam.tiff'];
    params_sub.crop_params.tiffStackout = f;
    ms.fileName = params_sub.crop_params.tiffStackout;
end
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
    if (sum(bad_vals)/length(ms_dt))  >.01
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
    if contains(ms.parentDir, 'Hipp16942/2022_06_10/18_25_10/')
            fprintf('\t%s\n', ms.parentDir)
            [ms] = cutoff_session(ms, tiff_numFrames);
    elseif contains(ms.parentDir, 'PKCZ_imaging/mHPC23454/2023_08_13/17_35_20_RET10/')
            fprintf('\t%s\n', ms.parentDir)
            [ms] = cutoff_session(ms, tiff_numFrames);
    elseif 'C:/Users/gjb326/Desktop/RecordingData/AlejandroGrau/TestMouse1/2022_07_05/17_09_41/'
            % do nothing, first 3 frames removed            
    else
            error('Unkown issue, update ''APA_troublesome_sessions.mat''')
    end
end
end

function [ms] = make_ori_struct(ms, msOriFile)
%%
global params_sub
ORI_Data = readtable(msOriFile);
shared_ts = ismember(ORI_Data.TimeStamp_ms_, ms.timestamps);

ori_ts = ORI_Data.TimeStamp_ms_(shared_ts);
w = ORI_Data.qw(shared_ts);
x = ORI_Data.qx(shared_ts);
y = ORI_Data.qy(shared_ts);
z = ORI_Data.qz(shared_ts);
[roll, pitch, yaw] = quaternion_conversion_LFOV(w, x ,y, z);


if false
q = [ORI_Data.qw, ORI_Data.qx, ORI_Data.qy, ORI_Data.qz];
q = q(shared_ts, :);
[~, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
ms.ori = ORI_Data(shared_ts, :);
end
% [ori_str] = LFOV_ori_extract(q, [0:30:360]);

ms.ori.time = ori_ts; 
ms.ori.roll = roll; 
ms.ori.pitch = pitch; 
ms.ori.yaw = yaw;

end