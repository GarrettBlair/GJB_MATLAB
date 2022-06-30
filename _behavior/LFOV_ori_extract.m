function [ori_str] = LFOV_ori_extract(q, degbin)

% ori = ms.ori;
% quat = ori(:,2:5);
% % quat_time = ori(:,1);
% % sharedq = ismember(quat_time, ms.time);
% % sm_quat = quat;
% % q = quat(sharedq,:);
% q = quat;

% [RotationMatrix, roll, pitch, yaw] = quatern2rotMat_Daniel(q);
[~, roll1, pitch1, yaw1] = quatern2rotMat_Daniel(q);
roll = rad2deg(roll1)+180; % offset to correct for position of BNO on camera
yaw = rad2deg(yaw1)+180; % offset to correct for position of BNO on camera
pitch = pitch1; % offset to correct for position of BNO on camera
pitch = pitch; % offset to correct for position of BNO on camera
roll = deg2rad(mod(roll+180, 360)-180);
roll = roll;
yaw = deg2rad(mod(yaw+30, 360)-180);
yaw = yaw;

ori_str.roll = roll;
ori_str.pitch = pitch;
ori_str.yaw = yaw;
if exist('degbin', 'var') && ~isempty(degbin)
[~,yawbin] = histc(yaw, degbin);
[~,rollbin] = histc(roll, degbin);
[~,pitchbin] = histc(pitch, degbin);
ori_str.rollbin = rollbin;
ori_str.pitchbin = pitchbin;
ori_str.yawbin = yawbin;
end




