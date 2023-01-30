function [SL, SR, LL, LR] = tad_linear_segmentation_rectmaze(x, y, plotting)
%%
t = 1:length(x); % t-t(1);
SL = []; SR = []; LL = []; LR = [];
plotting_sub = false;
short_segs = [100 -100];%100:-100:-100; % x values to sample crossings on
y_bounds = [50 -200]; % [max min] values
%%%%%%%%%%%%%%%% MAKE SL STRUCT
dir_flag = 1;
[SL] = make_pos_struct_temp(x, y, t, short_segs, y_bounds, dir_flag, plotting_sub);

%%%%%%%%%%%%%%%% MAKE SR STRUCT
dir_flag = -1;
[SR] = make_pos_struct_temp(x, y, t, short_segs, y_bounds, dir_flag, plotting_sub);

%%%%%%%%%%%%%%%% MAKE LL STRUCT as a composite of the three segments
y_segs = [50, 100]; 
dir_flag = 1;
x_bounds = [-100 -200]; % [max min] values
[LL.sub_neg] = make_pos_struct_temp(y, x, t, y_segs, x_bounds, dir_flag, plotting_sub);

dir_flag = 1;
y_bounds = [140 100]; % [max min] values
[LL.sub_mid] = make_pos_struct_temp(x, y, t, short_segs, y_bounds, dir_flag, plotting_sub);

dir_flag = -1;
x_bounds = [200 100];
[LL.sub_pos] = make_pos_struct_temp(y, x, t, y_segs, x_bounds, dir_flag, plotting_sub);
LL.index = LL.sub_neg.index | LL.sub_mid.index | LL.sub_pos.index;

LL.x = x; LL.y = y;
LL.x(~LL.index) = NaN; LL.y(~LL.index) = NaN;

%%%%%%%%%%%%%%%% MAKE LR STRUCT as a composite of the three segments
dir_flag = -1;
x_bounds = [-100 -200];
[LR.sub_neg] = make_pos_struct_temp(y, x, t, y_segs, x_bounds, dir_flag, plotting_sub);

dir_flag = -1;
[LR.sub_mid] = make_pos_struct_temp(x, y, t, short_segs, y_bounds, dir_flag, plotting_sub);

dir_flag = 1;
x_bounds = [200 100]; % [max min] values
[LR.sub_pos] = make_pos_struct_temp(y, x, t, y_segs, x_bounds, dir_flag, plotting_sub);

LR.index = LR.sub_neg.index | LR.sub_mid.index | LR.sub_pos.index;
LR.x = x; LR.y = y;
LR.x(~LR.index) = NaN; LR.y(~LR.index) = NaN;

if plotting
    figure(1); clf; hold on
    plot3(x, y, 1:length(x));
    plot3(SL.x, SL.y, 1:length(SL.index), 'b-', 'LineWidth', 3);
    plot3(SR.x, SR.y, 1:length(SR.index), 'r-', 'LineWidth', 3);
    plot3(LL.x, LL.y, 1:length(LL.index), 'c-', 'LineWidth', 3);
    plot3(LR.x, LR.y, 1:length(LR.index), 'm-', 'LineWidth', 3);
end
end
%%
function [POS_STRUCT] = make_pos_struct_temp(x, y, t, x_segs, y_bounds, dir_flag, plotting)
%%
%{
INPUTS: 
    x       - position value to segment along
    y       - position value, used to threshold x locations
    x_segs  - segments to find crossings along, steps through each; eg: x_segs = [100, -100]
    y_bounds- [int] bounding box pair values for valid x data; eg [50 -200], max and min values
    dir_flag- [1 or -1] swap the directionality to look for
    plotting- [logical] plot results
OUTPUT
    structure of position data with indices and positions
%}
dx = diff(x);
y_valid = y(1:end-1) < y_bounds(1) & y(1:end-1) >= y_bounds(2);
% t = 1:length(x);
POS_STRUCT.index= false(size(x)); %
POS_STRUCT.x_segs= x_segs;
POS_STRUCT.y_bounds= y_bounds;
POS_STRUCT.dir_flag= dir_flag;
% FIND BEGINNINGS
for threshLoop = 1:length(x_segs)
    if dir_flag==1
        temp=find(x(1:end-1)>x_segs(threshLoop) & x(2:end)<=x_segs(threshLoop) & y_valid);
    elseif dir_flag==-1
        temp=find(x(1:end-1)<x_segs(threshLoop) & x(2:end)>=x_segs(threshLoop) & y_valid);
    end
    temp=[1; temp; length(x)-1];
    for i=2:(length(temp)-1)
        if dir_flag==1
            d1=temp(i)-find(dx(temp(i):-1:temp(i-1))>0,1,'first');
        elseif dir_flag==-1
            d1=temp(i)-find(dx(temp(i):-1:temp(i-1))<0,1,'first');
        end
        
        if isempty(d1)
            d1=temp(i-1);
        end
        
        if dir_flag==1
            d2=temp(i)+find(dx(temp(i):temp(i+1))>0,1,'first')-1;
        elseif dir_flag==-1
            d2=temp(i)+find(dx(temp(i):temp(i+1))<0,1,'first')-1;
        end
        
        if isempty(d2)
            d2=temp(i+1);
        end
%         spd = sqrt([diff(x(d1:d2)).^2 + diff(y(d1:d2)).^2]);
%         spd_ends = nanmedian([spd(1:5); spd(end-5:end)]);
        
        POS_STRUCT.start_ind(i-1) = d1;
        POS_STRUCT.end_ind(i-1) = d2;
        POS_STRUCT.index(d1:d2)=true;
        POS_STRUCT.index(d1:d2)=true;
    end
end

POS_STRUCT.x= x; POS_STRUCT.x(~POS_STRUCT.index)=NaN;
POS_STRUCT.y= y; POS_STRUCT.y(~POS_STRUCT.index)=NaN;

if plotting
    t(~POS_STRUCT.index)=NaN;
    plot3(POS_STRUCT.x, POS_STRUCT.y, t, 'LineWidth', 2);
end

end
