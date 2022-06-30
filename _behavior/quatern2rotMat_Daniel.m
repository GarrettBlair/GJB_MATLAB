function [R, Roll, Pitch, Yaw] = quatern2rotMat_Daniel(q)

% m00 = 1.0 - 2.0.*qy.*qy - 2.0.*qz.*qz;
% m01 = 2.0.*qx.*qy + 2.0.*qz.*qw;
% m02 = 2.0.*qx.*qz - 2.0.*qy.*qw;
% m10 = 2.0.*qx.*qy - 2.0.*qz.*qw;
% m11 = 1 - 2.0.*qx.*qx - 2.0.*qz.*qz;
% m12 = 2.0.*qy.*qz + 2.0.*qx.*qw;
% m20 = 2.0.*qx.*qz + 2.0.*qy.*qw;
% m21 = 2.0.*qy.*qz - 2.0.*qx.*qw;
% m22 = 1.0 - 2.0.*qx.*qx - 2.0.*qy.*qy;
% 
% RotMat = NaN(3,3,length(qx));
% RotMat(1,:,:) = [m00 m01 m02]';
% RotMat(2,:,:) = [m10 m11 m12]';
% RotMat(3,:,:) = [m20 m21 m22]';
% 
% Roll = atan2(m12, m22);
% c2 = sqrt(m00.*m00 + m01.*m01);
% Pitch = atan2(-m02, c2);
% s1 = sin(Roll);
% c1 = cos(Roll);
% Yaw = atan2(s1.*m20 - c1.*m10, c1.*m11-s1.*m21);

%%
    R = zeros(3,3, size(q,1));
    R(1,1,:) = 1- 2.*q(:,3).^2 - 2.*q(:,4).^2;
    R(1,2,:) = 2.*(q(:,2).*q(:,3)+q(:,1).*q(:,4));
    R(1,3,:) = 2.*(q(:,2).*q(:,4)-q(:,1).*q(:,3));
    R(2,1,:) = 2.*(q(:,2).*q(:,3)-q(:,1).*q(:,4));
    R(2,2,:) = 1-2.*q(:,2).^2-2.*q(:,4).^2;
    R(2,3,:) = 2.*(q(:,3).*q(:,4)+q(:,1).*q(:,2));
    R(3,1,:) = 2.*(q(:,2).*q(:,4)+q(:,1).*q(:,3));
    R(3,2,:) = 2.*(q(:,3).*q(:,4)-q(:,1).*q(:,2));
    R(3,3,:) = 1-2.*q(:,2).^2-2.*q(:,3).^2;
    
    Roll = (atan2(R(2,3,:), R(3,3,:)));
    c2 = sqrt(R(1,1,:).*R(1,1,:) + R(1,2,:).*R(1,2,:));
    Pitch = (atan2(-R(1,3,:), c2));
    s1 = sin(Roll);
    c1 = cos(Roll);
    Yaw = atan2(s1.*R(3,1,:) - c1.*R(2,1,:), c1.*R(2,2,:)-s1.*R(3,2,:));
    Roll=squeeze(Roll);
    Pitch=squeeze(Pitch);
    Yaw=squeeze(Yaw);


%%
% [rows cols] = size(q);
% RotMat = zeros(3,3, rows);
% RotMat(1,1,:) = 1- 2.*q(:,3).^2 - 2.*q(:,4).^2;
% RotMat(1,2,:) = 2.*(q(:,2).*q(:,3)-q(:,1).*q(:,4));
% RotMat(1,3,:) = 2.*(q(:,2).*q(:,4)+q(:,1).*q(:,3));
% RotMat(2,1,:) = 2.*(q(:,2).*q(:,3)+q(:,1).*q(:,4));
% RotMat(2,2,:) = 1-2.*q(:,2).^2-2.*q(:,4).^2;
% RotMat(2,3,:) = 2.*(q(:,3).*q(:,4)-q(:,1).*q(:,2));
% RotMat(3,1,:) = 2.*(q(:,2).*q(:,4)-q(:,1).*q(:,3));
% RotMat(3,2,:) = 2.*(q(:,3).*q(:,4)+q(:,1).*q(:,2));
% RotMat(3,3,:) = 1-2.*q(:,2).^2-2.*q(:,3).^2;
% 
% Roll = squeeze(atan2(RotMat(1,3,:), RotMat(3,3,:)));
% c2 = squeeze(sqrt(RotMat(1,1,:).*RotMat(1,1,:) + RotMat(1,2,:).*RotMat(1,2,:)));
% Pitch = atan2(squeeze(-RotMat(1,3,:)), c2);
% s1 = sin(Roll);
% c1 = cos(Roll);
% Yaw = atan2(s1.*squeeze(RotMat(3,1,:)) - c1.*squeeze(RotMat(2,1,:)), c1.*squeeze(RotMat(2,2,:))-s1.*squeeze(RotMat(3,2,:)));

end




