ddirs = {'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/A_baseline_6mov/camkII_mpfc/2023_09_12/15_49_27/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/A_shocks_12mov/camkII_mpfc/2023_09_13/15_16_35/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/A_retrieval_24mov/camkII_mpfc/2023_09_15/15_14_12/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/B_baseline_6mov/camkII_mpfc/2023_09_12/16_16_30/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/B_shocks_12mov/camkII_mpfc/2023_09_13/15_43_47/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/B_retrieval_24mov/camkII_mpfc/2023_09_15/15_43_03/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/C_baseline_6mov/camkII_mpfc/2023_09_12/16_42_31/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/C_shocks_12mov/camkII_mpfc/2023_09_13/16_11_58/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/C_retrieval_24mov/camkII_mpfc/2023_09_15/16_15_02/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/D_baseline_6mov/camkII_mpfc/2023_09_25/16_11_42/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/D_shocks_12mov/camkII_mpfc/2023_09_26/15_51_39/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/D_retrieval_24mov/camkII_mpfc/2023_09_28/15_34_33/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/E_baseline_6mov/camkII_mpfc/2023_09_25/16_42_58/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/E_shocks_12mov/camkII_mpfc/2023_09_26/16_15_43/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/E_retrieval_24mov/camkII_mpfc/2023_09_28/16_05_37/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/F_baseline_6mov/camkII_mpfc/2023_09_12/17_07_17/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/F_shocks_12mov/camkII_mpfc/2023_09_13/16_38_17/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/F_retrieval_24mov/camkII_mpfc/2023_09_15/16_42_22/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/G_baseline_6mov/camkII_mpfc/2023_09_25/17_08_13/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/G_shocks_12mov/camkII_mpfc/2023_09_26/17_00_14/',...
        'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/miniscope_v4/Mov_Grid/G_retrieval_24mov/camkII_mpfc/2023_09_28/16_33_40/'};
savedir = 'E:/RecordingData/Fernando/';
cmn_name = 'caiman_cnmfe_out';
miniscopeName    = 'My_V4_Miniscope';
behav_name = 'extracted_pos';
webcamName    = 'My_WebCam';

params = [];
params.smin_vals                = -50:5:-10; % smin values used to create the 'deconv_sweep.mat'
% params.nan_interp               = true;
params.remove_all_bad_caiman_segs = false;
params.remove_segs_outside      = true;
params.correct_dt               = true; % correct for large jumps in timestamp file when constructing vmap
params.plotting                 = false;
params.skip_contour_bounding    = false;
params.cameraName               = {miniscopeName};%, 'ACC_miniscope2'}; % will default to 'MiniLFOV'
% params.reuse_contour_crop       = 'Crop_params.mat'; % use the previous ms file contour crop, unless empty
params.reuse_contour_crop       = 'bounding_box.mat'; % use the previous ms file contour crop, unless empty
%%
for i = 1:length(ddirs)
    %%
    recording_dir   = ddirs{i};
    f = ddirs{i}(84:end);
    inds = strfind(f, '/');
    sess_name = f(1:inds(1)-1);
    fname_ms = [savedir 'ms_data/' sess_name '.mat'];
    behav_file = [recording_dir webcamName '/' behav_name '.mat'];
    fname_new = [savedir 'ms_data_simple/' sess_name '.mat'];
    
    if isfile(fname_ms) && ~isfile(fname_new)
        %%
        if isfile(behav_file)
        fprintf('LOAD - %s\n', fname_ms)
        simple_file_subfcn(fname_ms, behav_file, fname_new)
        else
           warning('no behav %s', behav_file) 
        end
    else
        fprintf('SKIP - %s\n', fname_ms)
    end
end
%%
if false
for i = 1:length(ddirs)
    %%
    recording_dir   = ddirs{i};
    f = ddirs{i}(84:end);
    inds = strfind(f, '/');
    sess_name = f(1:inds(1)-1);
    contours_out = [savedir 'contours/' sess_name '.mat'];
    fname_out = [savedir 'ms_data/' sess_name '.mat'];
    
    caimanFilename = sprintf('%s%s/caiman_cnmfe_out.mat', recording_dir, miniscopeName);
    msTSFile = sprintf('%s%s/timeStamps.csv', recording_dir, miniscopeName);
    msCropFile = sprintf('%s%s/Crop_params.mat', recording_dir, miniscopeName);
    msOriFile = sprintf('%s%s/headOrientation.csv', recording_dir, miniscopeName);
    msBadFramesFile = sprintf('%s%s/badFrames.mat', recording_dir, miniscopeName);
    crop_params = load(msCropFile);
    params_sub = params;
    params_sub.crop_params  = crop_params;
    
    if ~isfile(fname_out)
        fprintf('LOAD - %s\n', caimanFilename)
        temp = load(caimanFilename);
        contours = gbContours(temp.fullA, temp.dims, [], .6);
%         ns = size(temp.C,1);
%         cm = jet(ns);
%         cm = cm(randperm(ns),:);
%         im = color_contours_im(c, [], temp.minFrame./255);
%         
%         figure(9); clf;
%         image(im);
%         drawnow
        
        % Read in the coresponding files
        crop_params = load(msCropFile);
        params_sub.crop_params  = crop_params;
        % Generate the miniscope structure
        if isfile(msTSFile)
            [ms] = make_ms_struct(recording_dir, msTSFile, miniscopeName);
        else
            error('No timestamp.csv file in directory:\n\t %s', recording_dir);
        end
        % Get the orientation data from BNO
        goodFrames = true(length(ms.frameNum),1);
        if isfile(msBadFramesFile)
            temp = load(msBadFramesFile);
            for j = 1:length(temp)
                idx1 = eval(sprintf('temp.f%d(1)', j-1));
                idx2 = eval(sprintf('temp.f%d(2)', j-1));
                goodFrames(idx1:idx2) = false;
            end
            ms.msBadFramesFile = msBadFramesFile;
        end
        ms.goodFrames = goodFrames;
        
        if isfile(msOriFile)
            [ms] = make_ori_struct(ms, msOriFile);
        else
            warning('No BNO data:\n\t %s', msOriFile);
            ms.ori = [];
        end
        
        ms.params               = params_sub;
        ms.cameraName           = miniscopeName;
        [ms]                    = extract_caiman_data(ms, ms.params, ms.cameraName);
        cmn = ms.neuron;
        ms = rmfield(ms, 'neuron');
        save(fname_out, 'ms')
        clearvars temp ms cmn
    else
        fprintf('skip %s\n', recording_dir)
    end
    
    if ~isfile(contours_out)
        fprintf('contours %s\n', contours_out)
%         temp = load(fname_out, 'cmn');
%         contours = gbContours(temp.cmn.fullA, temp.cmn.dims, [temp.cmn.idx_components], .6);
        temp = load(fname_out, 'ms');
        contours = gbContours(temp.ms.neuron.fullA, temp.ms.neuron.dims, [], .25, 1);
        save(contours_out, 'contours')
    else
        fprintf('skip %s\n', contours_out)
    end
end
end
%%
if false
cellregdirs = {'E:/RecordingData/Fernando/contours/cellreg_A/',...
    'E:/RecordingData/Fernando/contours/cellreg_B/',...
    'E:/RecordingData/Fernando/contours/cellreg_C/',...
    'E:/RecordingData/Fernando/contours/cellreg_D/',...
    'E:/RecordingData/Fernando/contours/cellreg_E/',...
    'E:/RecordingData/Fernando/contours/cellreg_F/',...
    'E:/RecordingData/Fernando/contours/cellreg_G/'};
n123 = NaN(length(cellregdirs), 3);
o123 = NaN(length(cellregdirs), 3);
p123 = NaN(length(cellregdirs), 3);
for i = 1:length(cellregdirs)
    cmapname = dir([cellregdirs{i} 'cellRegistered*']);
    temp = load([cmapname.folder '/' cmapname.name]);
    cmap = temp.cell_registered_struct.cell_to_index_map;
    nm = sum(cmap>0, 2);
    ns = size(cmap,1);
    n123(i,1) = sum(nm>=1);
    n123(i,2) = sum(nm>=2);
    n123(i,3) = sum(nm>=3);
    o123(i,1) = sum(nm==1);
    o123(i,2) = sum(nm==2);
    o123(i,3) = sum(nm==3);
    p123(i,1) = sum(nm>=1)/ns;
    p123(i,2) = sum(nm>=2)/ns;
    p123(i,3) = sum(nm>=3)/ns;
end
figure(3); clf; 
subplot(1,3,1); hold on
bar(nanmean(n123,1), 'FaceColor', [.7 .7 .7])
plot(n123', 'k')
axis([.25 3.75 -10 650])
ylabel('Number of cells matched')
set(gca, 'XTick', [1,2,3], 'XTickLabel', {'>=1', '>=2', '>=3'}, 'YTick', [100:100:600])
xlabel('# of Sessions')

subplot(1,3,2); hold on
bar(nanmean(o123,1), 'FaceColor', [.9 .7 .7])
plot(o123', 'r')
axis([.25 3.75 -10 450])
ylabel('Number of cells matched')
set(gca, 'XTick', [1,2,3], 'XTickLabel', {'=1', '=2', '=3'}, 'YTick', [100:100:400])
xlabel('# of Sessions')

subplot(1,3,3); hold on
bar(nanmean(p123,1), 'FaceColor', [.7 .7 .9])
plot(p123', 'b')
axis([.25 3.75 -.1 1.1])
ylabel('Percent of cells matched')
set(gca, 'XTick', [1,2,3], 'XTickLabel', {'>=1', '>=2', '>=3'}, 'YTick', [.2:.2:1])
xlabel('# of Sessions')
ttest(o123(:,1), o123(:,3))
end
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
ORI_Data = readtable(msOriFile);
shared_ts = ismember(ORI_Data.TimeStamp_ms_, ms.timestamps);

ori_ts = ORI_Data.TimeStamp_ms_(shared_ts);

q = [ORI_Data.qw, ORI_Data.qx, ORI_Data.qy, ORI_Data.qz];
q = q(shared_ts, :);
[~, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
ms.ori = ORI_Data(shared_ts, :);

ms.ori.time = ori_ts; 
ms.ori.roll = roll; 
ms.ori.pitch = pitch; 
ms.ori.yaw = yaw;

end

function simple_file_subfcn(ms_file, behav_file, newname)
%%
% ind1 = strfind(filename, 'processed_files')+length('processed_files') + 1;
% ind2 = strfind(filename, '.mat')+length('.mat') - 1;
% newname = ['D:\GarrettBlair\APA\HPCACC24500\simple_files\' filename(ind1:ind2)];
load(ms_file, 'ms')
behav = load(behav_file);

vars2save = {'parentDir' 'ms_file' 'behav_file' 'time_ms' 'dt_ms' 'good_ms_frames' 'frameNum'...
    'mousex' 'mousey' 'mouse_ecc' 'ledx' 'ledy' 'headangle' 'barx'...
    'mouse_spd_pxsec' 'led_spd_pxsec' 'bar_spd_pxsec' 'craw' 'spks'};


if contains(behav.beh.parentDir(4:end), ms.parentDir(5:end)) == true
parentDir = ms.parentDir;
time_ms = single(ms.timestamps);
dt_ms   = single(ms.dt);
good_ms_frames   = ms.goodFrames;
frameNum = behav.beh.frameNum+1;

behav_t = behav.beh.timestamps;
behav_dt = behav.beh.dt;
mousex   = single(behav.mouse.x);
mousey   = single(behav.mouse.y);
mouse_ecc   = single(behav.mouse.ecc);
ledx   = single(behav.led.x);
ledy   = single(behav.led.y);

% a = rad2deg(atan2( mousey-ledy, ledx-mousex));
headangle = rad2deg(atan2( ledy-mousey, ledx-mousex));

% for i = 1:4:10000;
%     figure(1); clf
%     hold on
%     plot(mousex(1:i), mousey(1:i), 'k')
%     plot(ledx(i), ledy(i), 'ro')
%     plot([ mousex(i) mousex(i)+cosd(a(i))*10], [ mousey(i) mousey(i)+sind(a(i))*10], 'r-')
%     axis([0 600 0 150])
%     title(a(i))
%     drawnow
%     pause(.1)
% end
barx   = single(behav.barpos(frameNum));

mouse_spd_pxsec = speed_calc_gb(mousex, mousey)./(behav_dt);
led_spd_pxsec = speed_calc_gb(ledx, ledy)./(behav_dt);
bar_spd_pxsec = speed_calc_gb(barx, barx*0)./(behav_dt);

mousex = interp1(behav_t, mousex, time_ms, 'linear');
mousey = interp1(behav_t, mousey, time_ms, 'linear');
mouse_ecc = interp1(behav_t, mouse_ecc, time_ms, 'linear');
ledx = interp1(behav_t, ledx, time_ms, 'linear');
ledy = interp1(behav_t, ledy, time_ms, 'linear');
headangle = interp1(behav_t, headangle, time_ms, 'linear');
barx = interp1(behav_t, barx, time_ms, 'linear');
mouse_spd_pxsec = interp1(behav_t, mouse_spd_pxsec, time_ms, 'linear');
led_spd_pxsec = interp1(behav_t, led_spd_pxsec, time_ms, 'linear');
bar_spd_pxsec = interp1(behav_t, bar_spd_pxsec, time_ms, 'linear');


raw_traces = ms.neuron.C + ms.neuron.YrA;
[dff_filt, ~] = ca_trace_dff_filter(raw_traces, round(1/median(dt_ms)), .1, 1.5, false);
spks = normalize_rows(dff_filt);
spks(isnan(spks)) = 0;
spks    = single(spks);
craw    = single(normalize_rows(raw_traces));

save(newname, vars2save{:})
else
    error('file mismatch!')
end



end