function [behav] = LFOV_generate_behav_struct(data_dir, fix_behav_position, nan_interp, rm_foreshortening, animalNum)
% nan_interp = true;
% rm_foreshortening = true;
behav.timestamp_file = sprintf('%s\\BehavCam_0\\timeStamps.csv', data_dir);
behav.tiff_filename = sprintf('%s\\BehavCam_0\\behavCam.tiff', data_dir);

% if strfind(dir, 'barrier')
%     avoid_barrier = true;
% else
%     avoid_barrier = false;
% end
btemp = xlsread(behav.timestamp_file);
tiffinfo = imfinfo(behav.tiff_filename);
numTiffFrames = size(tiffinfo,1);

behav.numFrames = size(btemp, 1);
behav.temporalDownsample = round(behav.numFrames/numTiffFrames);
behav.frameNum = btemp(1:behav.temporalDownsample:end, 1) + 1;
behav.time = btemp(1:behav.temporalDownsample:end, 2);
behav.height = tiffinfo(1).Height;
behav.width = tiffinfo(1).Width;
%%
bigF = zeros(length(behav.frameNum), behav.height*behav.width, 'uint8');
fprintf('Reading in behavior frames...')
nframes = min(numTiffFrames, length(behav.frameNum));
if numTiffFrames ~= length(behav.frameNum)
    fprintf('\n\t Discrepency found in the number of expected tiff images.\n\t\t Expected %d, found %d.\n Extracting the lesser... ', length(behav.frameNum), numTiffFrames)
    behav.frameNum = behav.frameNum(1:nframes);
    behav.time = behav.time(1:nframes);
end
for i = 1:nframes
    bigF(i,:) = reshape(imread(behav.tiff_filename, i), [1, behav.height*behav.width]);
end
fprintf(' Done!\n')
%% Extract body coordinates
frameSub = bigF(1:5:end,:);
led_thresh = quantile(frameSub(:), .9999)*.9;%.99985);
% body_thresh = quantile(frameSub(:), .00005);
se = fspecial('disk', 3);
% figure(9);
%%
min_area = 20; % Minimum area to conisder as target
max_area = 700; % Maximum area to conisder as target
max_accel = 30; % Maximum distance acceleration to consider
max_gap = 3; % Maximum frame gap to evaluate the distance metric
x = NaN(nframes, 1);
y = NaN(nframes, 1);
jump_dist = NaN(nframes, 1);
led_val = NaN(nframes, 1);
%a = quantile(bigF(:), .05);
fprintf('Extracting position using LED...')

% figure 
for i = 1:nframes
    f = reshape(squeeze(bigF(i,:)), [behav.height, behav.width]);
    ft = f;
    rat_roi = f>led_thresh;% Looking for bright spots works best
    smooshed_rat = filter2(se, rat_roi)>0; % smooth out slightly
    ft(smooshed_rat==0) = 0;
    im = regionprops(smooshed_rat, ft, {'Area', 'Centroid', 'MaxIntensity'});
    a = min_area;
    val = led_thresh;
    for j = 1:size(im, 1) % search through the regions to find
        % max area within min_area and max_area and the brightest spot
        if im(j).Area > a && im(j).Area <= max_area && im(j).MaxIntensity > val
            xn = im(j).Centroid(1);
            yn = im(j).Centroid(2);
            % check the distance between this centroid and the last known
            % position
            if i==1
                xprev = xn;
                yprev = yn;
                jprev = sqrt((xn-xprev)^2 + (yn-yprev)^2);
            else
                last_good = find(~isnan(x(1:i-1)), 1, 'last');
                xprev = x(last_good);
                yprev = y(last_good);
                jprev = jump_dist(last_good);
                % if no last know, or is greater than max_gap samples away,
                % take this spot
                if isempty([xprev yprev]) || (i-last_good)>=max_gap
                    xprev = xn;
                    yprev = yn;
                    jprev = sqrt((xn-xprev)^2 + (yn-yprev)^2);
                end
            end
            % dist between current and last known
            jn = sqrt((xn-xprev)^2 + (yn-yprev)^2);
            if  abs(jn-jprev) <= max_accel % if change acceleration is below threshold, accept. 
                %This helps to eliminate large jumps while not punishing
                %running periods
                x(i) = xn;
                y(i) = yn;
                a = im(j).Area;
                val = im(j).MaxIntensity;
                jump_dist(i) = jn;
                led_val(i) = val;
            end
        end
    end
%     if i > (nframes*flag)
%         fprintf('%d%%..', round(flag*100))
%         flag = flag+.1;
%     end
% cla
% imagesc(f, [0 255]); axis image
% hold on
% plot(x(i), y(i), 'ro')
% title(i)
% drawnow
end
%%
fprintf(' Done!\n')
behav.x = x;
behav.y = y;
behav.spd = sqrt(diff([x(1); x]).^2 + diff([y(1); y]).^2);

behav.raw_x = behav.x;
behav.raw_y = behav.y;
if fix_behav_position % correct for camera position and angle
    x = behav.raw_x;
    y = behav.raw_y;
    t = behav.time;
    
    if nan_interp
        nanind = (isnan(x) & isnan(y));
        xn = interp1(find(~nanind), x(~nanind), find(nanind), 'linear');
        yn = interp1(find(~nanind), y(~nanind), find(nanind), 'linear');
        x(nanind) = xn; y(nanind) = yn;
        behav.wasnan = nanind;
    end
    behav.wasnan = nanind;
    % 4 degree rotation about x to fix cam
    % This is specific to the Hipp5-18 animals on the scary maze from early
    % 2020
    if animalNum <30
        c = rotx(4);
    else % camera changed slightly between animals
        c = rotx(5);
    end
    c = c(2:end, 2:end);
    xn = [x y]*c;
    x = xn(:,1); y = xn(:,2);
    

    
    if rm_foreshortening  % correct for camera distortion
        %% Scaling position to remove foreshortening
        xo = x; yo = y;
        if animalNum < 30 % camera changed slightly between animals
            xcenter = 338; % center x,y point on the short path
            ycenter = 190; % center x,y point on the short path
            ymin = 190-ycenter;
            yextreme = 346-ycenter;
            xmin = 120-xcenter;
            xextreme = 220-xcenter;
        else
            xcenter = 350; % center x,y point on the short path
            ycenter = 150; % center x,y point on the short path
            ymin = 190-ycenter;
            yextreme = 346-ycenter;
            xmin = 120-xcenter;
            xextreme = 220-xcenter;
        end
        xo = xo-xcenter;
        yo = yo-ycenter;
        xs = (xextreme - xmin)/(xextreme-xcenter);
        ys = (yo-ymin)/(yextreme-ymin);
        xo = xo +  xo*xs.*ys;
        
        % Shorten y and x dims
        x = 250*xo/440;
        y = 125*yo/(149-10);
    end
    figure; plot3(x,y, 1:length(behav.time))
    set(gcf, 'Name', data_dir)
    drawnow
    %%
    behav.x = x;
    behav.y = y;
    behav.t = t;
    
    behav_dt = [median(diff(behav.time)); diff([behav.time])]/1000;
    behav.spd = sqrt(diff([x(1); x]).^2 + diff([y(1); y]).^2)./behav_dt;


end