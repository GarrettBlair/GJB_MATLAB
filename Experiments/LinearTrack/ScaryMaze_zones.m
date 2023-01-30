function [short, long, left, right] = ScaryMaze_zones(x,y);
%% Draw what points (in cm on maze) are valid

left.bounds = ...
    [-120.01 -30;
    -120.01 15;
    -200 15;
    -200 -30;
    -120.01 -30]; % valid points to include in left goal, Rectanglular shock maze
right.bounds = ...
    [120.01 -30;
    120.01 15;
    200 15;
    200 -30;
    120.01 -30]; % valid points to include in right goal, Rectanglular shock maze
% used in Scarymaze paper
% short.bounds  = ...
%     [-120 -30;
%     -120 5;
%     120 5;
%     120 -30;
%     -120 -30]; % valid points to include in long path, Rectanglular shock maze
% expanded slightly vertical for supp data with different camera a angle
short.bounds  = ...
    [-120 -30;
    -120 15;
    120 15;
    120 -30;
    -120 -30]; % valid points to include in long path, Rectanglular shock maze
% long.bounds  = ...
%     [-200   15.01;
%     -130    15.01;
%     -130    50;
%     130     50
%     130     15.01
%     200   15.01;
%     200   200;
%     -200   200;
%     -200   15.01]; % valid points to include in short path, Rectanglular shock maze
long.bounds  = ...
    [-200   15.01;
    -120    15.01;
    -120    50;
    120     50
    120     15.01
    200   15.01;
    200   200;
    -200   200;
    -200   15.01]; % valid points to include in short path, Rectanglular shock maze


right.inside = inpolygon(x, y, right.bounds(:,1), right.bounds(:,2));
left.inside = inpolygon(x, y, left.bounds(:,1), left.bounds(:,2));
short.inside  = inpolygon(x, y,  short.bounds(:,1),  short.bounds(:,2));
long.inside  = inpolygon(x, y,  long.bounds(:,1),  long.bounds(:,2));

% figure; plot(x,y)
% hold on; plot(left.bounds(:,1), left.bounds(:,2))
% hold on; plot(right.bounds(:,1), right.bounds(:,2))
% hold on; plot(long.bounds(:,1), long.bounds(:,2))
% hold on; plot(short.bounds(:,1), short.bounds(:,2))
% drawnow