function [behav] = circular_pos_tracking(behav, radius)
if ~isfield(behav, 'raw_position')
    raw_position = NaN(behav.numFrames, 2);
    for frameNum = 1:behav.numFrames;
        frame = double(msReadFrame(behav,frameNum,false,false,false))/255; % grab frame
        temppp = (frame(:,:,1)) > .7 & frame(:,:,2) > .1 & frame(:,:,3) > .1;
        props = regionprops(temppp,'Area','Centroid');
        if ~isempty(props)
            a = 0;
            for i=1:length(props)
                if (props(i).Area > a)
                    a = props(i).Area;
                    index  = i;
                end
            end
            if index ~=0
                raw_position(frameNum,[1 2]) = props(index).Centroid;
            end
            
        else
            
        end
    end
    
    ux = raw_position(:,1);
    uy = raw_position(:,2);
    behav.raw_position = raw_position;
else
    ux = behav.raw_position(:,1);
    uy = behav.raw_position(:,2);
end

resp = 'n';
while ~strcmp(resp, 'y')
    figure(1);
    clf
    hold on
    plot(behav.raw_position(:,1), behav.raw_position(:,2), 'k.');
    if exist('xv', 'var')
        plot(xv,yv,'r.')
        plot(center(1), center(2), '*m')
    end
    center = ginput(1);
    
    L = linspace(0, 2.*pi, 360);
    xv = cos(L)';
    yv = sin(L)';
    
    xv = (xv*radius + center(1));
    yv = (yv*radius + center(2));
    
    plot(xv,yv,'r.')
    hold off
    resp = input('keep? y/n   ', 's');
end
% [in,on] = inpolygon(ux, uy, xv, yv);
[c,q] = ginput(4);
[in,on] = inpolygon(ux, uy, c, q);
behav.position = NaN(length(ux), 2);
behav.position(in | on,1) = ux(in | on);
behav.position(in | on,2) = uy(in | on);

