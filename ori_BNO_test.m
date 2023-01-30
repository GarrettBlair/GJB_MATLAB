%%
load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\2022_09_14_H17_52_22_WTR8_@placecells.mat')
x = ms.room.x; y = ms.room.y; t = ms.timestamps./1000;
yaw = ms.ori.yaw;
% ind = t>10 & t<25;
% x = x(ind); y = y(ind); t = t(ind); yaw = yaw(ind);

yaw_offset = -pi/2;
win = 100
for i = 1:300%length(x)-100
    %%
    figure(1); clf
    hold on
    plot3([x(i:i+win)], [y(i:i+win)], [t(i:i+win)], 'k.-')
%     plot([x(i), x(i+1)], [y(i), y(i+1)], 'k.-')
    xd = cos(yaw(i:i+win)+yaw_offset);
    yd = sin(yaw(i:i+win)+yaw_offset);
    for j = 1:length(xd)
    plot3([x(i+j-1) x(i+j-1)+xd(j)], [y(i+j-1) y(i+j-1)+yd(j)], [t(i+j-1) t(i+j-1)], 'r-')
    end
    axis([-40 40 -40 40 0 100])
    drawnow
    pause(.025)
%     plot([x(i), x(i)+xd], [y(i), y(i)+yd], 'r-')
end

%%
recording_dir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\2022_12_14\15_25_31\';
msOriFile = sprintf('%sMiniLFOV/headOrientation.csv', recording_dir);
msTSFile = sprintf('%sMiniLFOV/timeStamps.csv', recording_dir);
%     [mstemp] = make_ms_struct(recording_dir, msTSFile);
%     [mstemp] = make_ori_struct(ms, msOriFile);

ORI_Data = readtable(msOriFile);
TS_data = readtable(msTSFile);
mstemp.timestamps = TS_data.TimeStamp_ms_;%(1:ms.temporalDownsample:end);
% shared_ts = ismember(ORI_Data.TimeStamp_ms_, mstemp.timestamps);

ori_ts = ORI_Data.TimeStamp_ms_;
w = ORI_Data.qw;
x = ORI_Data.qx;
y = ORI_Data.qy;
z = ORI_Data.qz;
q = [ORI_Data.qw, ORI_Data.qx, ORI_Data.qy, ORI_Data.qz];
% q = q(shared_ts, :);
[~, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
[roll2, pitch2, yaw2] = quaternion_conversion_LFOV(w, x ,y, z);
mstemp.ori = ORI_Data;%(shared_ts, :);
% [ori_str] = LFOV_ori_extract(q, [0:30:360]);

mstemp.ori.time = ori_ts; 
mstemp.ori.roll = roll; 
mstemp.ori.pitch = pitch; 
mstemp.ori.yaw = yaw;

figure;
hold on;
% plot(ori_ts./1000, roll, 'r-');
% plot(ori_ts./1000, pitch, 'm-');
% plot(ori_ts./1000, yaw, 'y-')
plot(ori_ts./1000, roll2, 'r-');
plot(ori_ts./1000, pitch2, 'm-');
plot(ori_ts./1000, yaw2, 'y-')

