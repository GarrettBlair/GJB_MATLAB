

topdir = 'G:\.shortcut-targets-by-id\1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY\miniscope_v4\Mov_Grid\D_baseline_6mov\camkII_mpfc\2023_09_25\16_11_42\My_WebCam\';
topdir = 'G:\.shortcut-targets-by-id\1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY\miniscope_v4\Mov_Grid\D_shocks_12mov\camkII_mpfc\2023_09_26\15_51_39\My_WebCam\';
vids = [0:21];
msTSFile = [topdir 'timeStamps.csv'];
[beh] = make_ms_struct(topdir, msTSFile);
px = []; py = [];
%%
numv = length(vids);
vid = VideoReader( [topdir num2str(vids(1)) '.avi']);
Yr = uint8(zeros(vid.height, vid.width, numv*vid.NumFrames));
Yg = uint8(zeros(vid.height, vid.width, numv*vid.NumFrames));
Yb = uint8(zeros(vid.height, vid.width, numv*vid.NumFrames));
step = 1;
for i = 1:length(vids)
    vname = [topdir num2str(vids(i)) '.avi'];
    vid = VideoReader(vname); %#ok<TNMLP>
    disp([num2str(vids(i))])
    frames2read = min([1000, vid.NumFrames]);

    % while v.hasFrame
    for ind = 1:step:frames2read % 1:step:726
        fnum = (i-1)*1000 + ind;
        im = vid.read(ind);
        Yr(:,:,fnum) = im(:,:,1);
        Yg(:,:,fnum) = im(:,:,2);
        Yb(:,:,fnum) = im(:,:,3);
        lastf = fnum;
    end
end
Yr = Yr(:,:,1:lastf);
Yg = Yg(:,:,1:lastf);
Yb = Yb(:,:,1:lastf);
% Yrbg = std(single(Yr), [], 3) + mean(Yr, 3);
Yrbg = movmean(Yr, 30, 3);
Ygbg = movmean(Yg, 30, 3);
Ybbg = movmean(Yb, 30, 3);
for fnum = 5000:4:lastf; figure(1); im = cat(3, Yr(:,:,fnum), Yg(:,:,fnum), Yb(:,:,fnum)); bg = cat(3, Yrbg(:,:,fnum), Ygbg(:,:,fnum), Ybbg(:,:,fnum)); image([im; single(im) - bg]); drawnow; end
% Yrbg = min(Yr,[], 3);
%%
led_sz =3; % radius
mouse_sz = 8; % radius
grid_sz = [80 40]; % h x w

    


LED_se = strel('disk', led_sz);
hsmall = gausswin(led_sz); hsmall  = hsmall*hsmall'; hsmall = hsmall./sum(hsmall(:));
mouse_se = strel('disk', mouse_sz);
grid_se = strel('rectangle', [grid_sz]);
line = strel('rectangle', [3,1]);


led_x = NaN(numv*1000,1);
led_y = NaN(numv*1000,1);
led_val = NaN(numv*1000,1);
led_dist = NaN(numv*1000,1);
led_size = NaN(numv*1000,1);
mouse_x = NaN(numv*1000,1);
mouse_y = NaN(numv*1000,1);
mouse_ecc = NaN(numv*1000,1);
mouse_size = NaN(numv*1000,1);
headangle = NaN(numv*1000,1);
barpos = NaN(numv*1000,1);


%
minbararea = 60;
led_sizemin = 2;
led_sizemax = 16;
led_valmin = 200;
minmousearea = 100;
searchsize = 50;
%
plotting=true;
tic
step = 1;
plotstep = 1:10:numv*1000;
leddistthresh = 2;
strict_dist_min = 6;
lastfnum =[]; px = [];

mouse_gauss = gausswin(2*searchsize + 1, 1.5); mouse_gauss = mouse_gauss*mouse_gauss'; mouse_gauss = mouse_gauss./(sum(mouse_gauss(:)));
% mouse_gauss = mouse_gauss.^2;
px=[]
% px = [32   611];
% py =[40   141];
% lastx=26; lasty=52; lastfnum=-2;
% vo = VideoWriter("avi_project\mouse_track_demo");
% vo.FrameRate = 30;
% vo.open();

for fnum = 1:step:lastf
        if isempty(px)
            figure;
            image(im)
            [px, py] = ginput(2);
            px = round([min(px) max(px)]);
            py = round([min(py) max(py)]);
            
        % end
        % if isempty(lastfnum)
        %     figure;
        %     image(im)
            [lastx, lasty] = ginput(1);
            lastx = round(lastx)-px(1);
            lasty = round(lasty)-py(1);
            lastfnum = -3;

        end
        im = cat(3, Yr(:,:,fnum), Yg(:,:,fnum), Yb(:,:,fnum));
        im = im(py(1):py(2), px(1):px(2), :);
        bg = cat(3, Yrbg(:,:,fnum), Ygbg(:,:,fnum), Ybbg(:,:,fnum));
        bg = bg(py(1):py(2), px(1):px(2), :);
        imr = double(im(:,:,1));
        img = double(im(:,:,2));
        imb = double(im(:,:,3));
        % imr = conv2(imr, hsmall, 'same');
        validLED = false(size(imr));
        colorLED = zeros(size(imr));
        imbar = mean(double(im), 3)>230;

        h = median(imr, 2);
        w = median(imr, 1);
        m = h*w; m = max(imr(:))*m./(max(m(:)));
        % bg = conv2(im(:,:,1), ones(15, 15)./225, 'same');
        % mim = mean(double(im), 3);
        immouse = conv2(imr, 1*mouse_se.Neighborhood, 'same');
        rdiffmouse = m - immouse;
        rdiff_mouse = rdiffmouse >= quantile(rdiffmouse(:), .99);
        % rdiff_mouse = conv2(rdiff_mouse, ones(3,3)./9, 'same')>0;
        h = median(imb, 2);
        w = median(imb, 1);
        m = h*w; m = max(imb(:))*m./(max(m(:)));
        rdiff = m - imb;
        % rdiff_LED = rdiff.*imr;
        % rdiff_LED = rdiff_LED <= quantile(rdiff_LED(:), .001);
        % 
        % % rdiff_LED = rdiff <= quantile(rdiff(:), .01);
        % % rdiffLED = conv2(rdiff, mean(ledkern,3), 'same');
        % 
        % 
        % 
        % imled = imclose(rdiff_LED, LED_se);
        % imled = imerode(imled, line);
        immouse = imclose(rdiff_mouse, mouse_se);
        imbar = imclose(imbar, grid_se);
        a = regionprops(imbar, {'Area', 'Eccentricity', 'Centroid'});
        if ~isempty(a)
            bigind = find([a.Area] == max([a.Area]), 1);
            if a(bigind).Area >= minbararea
                barpos(fnum) = a(bigind).Centroid(1);
            end
        end


        a = regionprops(immouse, {'Area', 'Eccentricity', 'Centroid'});
        if ~isempty(a)

            bigind = find([a.Area] == max([a.Area]), 1);
            mousefound = false;
            if a(bigind).Area >= minmousearea
                mousefound = true;
                mouse_x(fnum) = a(bigind).Centroid(1);
                mouse_y(fnum) = a(bigind).Centroid(2);
                mouse_ecc(fnum) = a(bigind).Eccentricity;
                mouse_size(fnum) = a(bigind).Area;

                if (fnum - lastfnum) > 30
                    lastx = []; lasty =  [];
                    lastfnum = [];
                end

                ww = round(mouse_x(fnum))-searchsize:round(mouse_x(fnum))+searchsize;
                hh = round(mouse_y(fnum))-searchsize:round(mouse_y(fnum))+searchsize;
                subG = mouse_gauss;
                if ~isempty(lastx) && abs(lastx - mouse_x(fnum))<searchsize && abs(lasty - mouse_y(fnum))<searchsize
                    % ww = round(lastx-searchsize/2:lastx+searchsize/2);
                    % hh =  round(lasty-searchsize/2:lasty+searchsize/2);
                    xx = mean([lastx mouse_x(fnum)]);
                    yy = mean([lasty mouse_y(fnum)]);
                    ww = round(xx-searchsize/2:xx+searchsize/2);
                    hh =  round(yy-searchsize/2:yy+searchsize/2);
                    subG = imresize(subG, .5);
                end
                www = find(ww>1 & ww <= size(imb,2));
                hhh = find(hh>1 & hh <= size(imb,1));
                ww = ww(ww>1 & ww <= size(imb,2));
                hh = hh(hh>1 & hh <= size(imb,1));

                rdiff_LED = rdiff.*imr;
                rdiff_LED = rdiff_LED(hh,ww).*subG(hhh,www);
                % figure(13); clf; imagesc(rdiff_LED);
                        % rdiff_LED = rdiff_LED <= min(rdiff_LED(:))*.90;

                rdiff_LED = rdiff_LED <= quantile(rdiff_LED(:), .005);
                imled = imclose(rdiff_LED, LED_se) == 1;
                validLED(hh,ww) = imled;
                % validLED(hh,ww) = imled(hh,ww)==1;
                % colorLED(validLED==1) = imb(validLED==1);

                a = regionprops(validLED==1, {'Area', 'Centroid', 'PixelIdxList'});
                if ~isempty(a)
                    % brightest
                    % vals = NaN(length(a),1);
                    % for aa = 1:length(a)
                    %     avVal = max(max(imr( a(aa).PixelIdxList )));
                    %     vals(aa) = avVal;
                    % end
                    % dists = abs(255-vals);

                    % closest
                    c = cat(1, a.Centroid);
                    if ~isempty(lastx) && (fnum - lastfnum)<=30 %% fnum>step+1 && 
                    % ax = led_x(fnum-step) - c(:,1);
                    % ay = led_y(fnum-step) - c(:,2);
                    ax = lastx - c(:,1);
                    ay = lasty - c(:,2);
                    prevdists = sqrt( (ax.^2) + (ay.^2) );
                    prevdists(prevdists>leddistthresh*(fnum - lastfnum) & prevdists > strict_dist_min) = NaN;
                    else
                    prevdists = ones(length(a),1);
                    lastx=[]; lasty=[]; 
                    end
                    ax = mouse_x(fnum) - c(:,1);
                    ay = mouse_y(fnum) - c(:,2);
                    
                    dists = sqrt( (ax.^2) + (ay.^2) );
                    dists_heuristic = dists.*prevdists;
                    if any(~isnan(dists_heuristic))
                    [~, ind_ord] = sort(dists_heuristic, 'ascend');
                    aa=1; ledfound=false;
                    while ledfound==false && (aa<= length(a))% for aa = 1:length(a)
                        idx = ind_ord(aa);
                        thissize = a(idx).Area;
                        % avVal = mean(mean(imr( a(idx).Image==1 )));
                        avVal = mean(mean(imr( a(idx).PixelIdxList )));
                        if thissize >= led_sizemin && thissize <= led_sizemax && avVal >= led_valmin % imr(a(idx).Image)
                            led_x(fnum) = a(idx).Centroid(1);
                            led_y(fnum) = a(idx).Centroid(2);
                            lastx = led_x(fnum); 
                            lasty = led_y(fnum); lastfnum = fnum;
                            led_val(fnum) = avVal;
                            led_dist(fnum) = dists(idx);
                            led_size(fnum) = a(idx).Area;
                            ledfound=true;

                            x0 = mouse_x(fnum);
                            y0 = mouse_y(fnum);
                            x1 = led_x(fnum);
                            y1 = led_y(fnum);
                            x2 = mouse_x(fnum)+ sqrt( (x0-x1).^2 + (y0-y1).^2 );
                            y2 = mouse_y(fnum);

                            x10 = x1-x0;
                            % y10 = y1-y0;
                            y10 = y0-y1; % because the pxal convention swaps the yaxis
                            x20 = x2-x0;
                            y20 = y2-y0;
                            % x10 = x0-x1;
                            % y10 = y0-y1;
                            x21 = x2-x1;
                            y21 = y2-y1;
                            % ang = rad2deg( atan(sin(y10)/cos(x10)) ); 
                            ang =  rad2deg( atan2(y10, x10) ); 

                            % disp([x10 y10])

                            % if x10>0 && y10>0
                            %     ang = ang+90;
                            % elseif x10>0 && y10<0
                            %     ang = ang+180;
                            % elseif x10<0 && y10<0
                            %     ang = ang+270;
                            % end
                            
                            % ang = rad2deg( atan2(abs(x10*y20-x20*y10),x10*y10+x20*y20) );


                            headangle(fnum) = ang;

                        end
                        aa=aa+1;
                    end
                    aa = aa-1;
                    end
                end
            end
        end
        % if fnum < 3000 && fnum > 1000
        %     plotting = true;
        % else
        %     plotting = false;
        % end
        if plotting==true && any(ismember(fnum, plotstep))
            %%
            
            figure(1);clf; 
            set(gcf, 'Position', [25   528   759   230])
            subplot_tight(2,3,1:2)
            image(im); axis image; title(['Frame  ' num2str(fnum)]); axis image off
            subplot_tight(2,3,4:5)
            imL = im; imL = imL+.6*(255-imL);
            image(imL); axis image off; %title(num2str(fnum))
            hold on
            plot(mouse_x(~isnan(mouse_x)), mouse_y(~isnan(mouse_y)), '-', 'Color', [.2 .2 .2])
            plot(led_x(~isnan(led_x)), led_y(~isnan(led_y)), 'r-')
            plot([barpos(fnum) barpos(fnum)], [0 size(im,1)], 'm--')
            % plot(barpos(1:fnum), 150, 'm-')
            % subplot(2,2,2)
            % imagesc(colorLED); axis image; title('LED')
            % hold on
            % scatter(led_x(fnum), led_y(fnum), 100, 'm')
            % subplot(2,2,3)
            % imagesc(immouse); axis image; title('mouse')
            % hold on
            % scatter(mouse_x(fnum), mouse_y(fnum), 100, 'mo')
            subplot_tight(1,3, 3, [.05 .125])
            % imagesc(imbar); axis image; title('bar')
            % hold on
            % scatter(barpos(fnum), size(im,1)/2, 100, 'mo')
            % drawnow;
            isub = im;
                        imL = im; imL = imL+.2*(255-imL);
            image(imL); hold on; axis image
            % plot([0 x10]-searchsize+size(isub,2), [0 -y10]-searchsize+size(isub,1), 'ro-');
            scatter([mouse_x(fnum), led_x(fnum)], [mouse_y(fnum), led_y(fnum)], 'ko', 'MarkerFaceColor', 'k')
            plot([mouse_x(fnum), led_x(fnum)], [mouse_y(fnum), led_y(fnum)], 'r-', 'LineWidth', 2)
            % plot([0 x20]-searchsize+size(isub,2), [0 y20]-searchsize+size(isub,1), 'ko-');
            ha = headangle(fnum);
            if ~isnan(ha)
            title(sprintf('Head Ang %4.fdeg', ha))
            else
            title(sprintf('Head Ang      %sdeg', ' '))
            end
            axis image off % ([0 searchsize*2 0 searchsize*2])
                ww = round(mouse_x(fnum))-searchsize:round(mouse_x(fnum))+searchsize;
                hh = round(mouse_y(fnum))-searchsize:round(mouse_y(fnum))+searchsize;
                ww = ww(ww>1 & ww <= size(imb,2));
                hh = hh(hh>1 & hh <= size(imb,1));
            axis([ww(1), ww(end) hh(1) hh(end)])
            drawnow
            % temp = getframe(gcf);
            % frameout = temp.cdata;
            % vo.writeVideo(frameout)
        end
end


toc

%%
ha = headangle;
ha(isnan(ha))= interp1(find(~isnan(ha)), ha(~isnan(ha)), find(isnan(ha)), 'nearest');
da = abs(diff(unwrap(deg2rad(ha))));
da = [da(1); da];






%% topdir = 'G:\.shortcut-targets-by-id\1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY\miniscope_v4\Mov_Grid\';
clear
topdir = 'G:\.shortcut-targets-by-id\1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY\miniscope_v4\Mov_Grid\ms_data\';
names = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
types = {'_baseline_6mov' '_shocks_12mov' '_retrieval_24mov'};
numA= length(names);
numsessper = length(types);

nfiles = numA*numsessper;
%
if false
for i  = 1:numA
    for j = 1:numsessper
        fname = [topdir names{i} types{j} '.mat'];
        disp([names{i} types{j} '.mat'])
        clearvars ms
        load(fname)
        raw = ms.neuron.YrA + ms.neuron.C;
        dt = 1/median(abs(diff(ms.timestamps./1000)));
        [spks, ~] = ca_trace_dff_filter(raw, dt, .1, 1.5, true);
        drawnow;
        spks(~ms.goodFrames) = 0;
        spks = normalize_rows(spks);   
        ms.spks = spks;
        save(fname, 'ms', '-append')
    end
end
end

%%
% spks_1secbin = bin_spks_time(ms.spks>0, 1, ms.timestamps./1000, false);
% spks_1secbin(isnan(spks_1secbin)) = 0;
% [c,p] = corr(spks_1secbin','Type','Pearson');
% segs_corr = [];
% segs_corr.c = c;
% segs_corr.p = p;
% segs_corr.tau = '1sec';
% segs_corr.binspks = spks_1secbin;

function [ms] = make_ms_struct(recording_dir, msTSFile)
TS_data = readtable(msTSFile);
ms = [];

recording_dir(strfind(recording_dir, '\')) = '/';

ms.parentDir = recording_dir;
ms.frameNum = TS_data.FrameNumber;
ms.timestamps = TS_data.TimeStamp_ms_;
ms_dt = [median(diff(ms.timestamps)); diff([ms.timestamps])]/1000;
ms.dt = ms_dt;
end