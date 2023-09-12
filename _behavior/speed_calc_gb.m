function [spd] = speed_calc_gb(x,y)
a(:,1) = x;
b(:,1) = y;

spd = sqrt( diff(a).^2 + diff(b).^2 );
spd = [spd(1); spd];
