function [roll2, pitch2, yaw2] = quaternion_conversion_LFOV(w, x ,y, z)

% recording_dir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\2022_12_13\12_47_07\';
% msOriFile = sprintf('%sMiniLFOV/headOrientation.csv', recording_dir);
% msTSFile = sprintf('%sMiniLFOV/timeStamps.csv', recording_dir);
% %     [mstemp] = make_ms_struct(recording_dir, msTSFile);
% %     [mstemp] = make_ori_struct(ms, msOriFile);
% 
% ORI_Data = readtable(msOriFile);
% TS_data = readtable(msTSFile);
% mstemp.timestamps = TS_data.TimeStamp_ms_;%(1:ms.temporalDownsample:end);
% shared_ts = ismember(ORI_Data.TimeStamp_ms_, mstemp.timestamps);
% 
% ori_ts = ORI_Data.TimeStamp_ms_(shared_ts);
% q = [ORI_Data.qw, ORI_Data.qx, ORI_Data.qy, ORI_Data.qx];
% q = q(shared_ts, :);
% % [~, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
%%
q = [];
q.w = w;
q.x = x;
q.y = y;
q.z = z;

% % Conversion according to
% % https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Source_code_2
% roll (x-axis rotation)
sinr_cosp = 2 * (q.w .* q.x + q.y .* q.z);
cosr_cosp = 1 - 2 * (q.x .* q.x + q.y .* q.y);
roll2 = atan2(sinr_cosp, cosr_cosp);

% pitch (y-axis rotation)
sinp = 2 * (q.w .* q.y - q.z .* q.x);
pitch2 = asin(sinp);
% if abs(sinp) >= 1
% %     pitch2 = copysign(pi / 2, sinp); % use 90 degrees if out of range
%     pitch2 = sign(sinp)*pi/2; % use 90 degrees if out of range
% else
%     pitch2 = asin(sinp);
% end
% yaw (z-axis rotation)
siny_cosp = 2 * (q.w .* q.z + q.x .* q.y);
cosy_cosp = 1 - 2 * (q.y .* q.y + q.z .* q.z);
yaw2 = atan2(siny_cosp, cosy_cosp);

% yaw2 = yaw2+pi;
% yaw2(yaw2>pi) = (yaw2(yaw2>pi)-pi).*-1;
% roll2 = roll2+pi;
% roll2(roll2>pi) = (roll2(roll2>pi)-pi).*-1;

% figure;
% hold on;
% plot(roll2, 'r.');
% plot(pitch2, 'm.');
% plot(yaw2, 'y.')