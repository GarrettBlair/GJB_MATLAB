if false

topdir = 'G:\.shortcut-targets-by-id\1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY\miniscope_v4\Mov_Grid\D_baseline_6mov\camkII_mpfc\2023_09_25\16_11_42\My_WebCam\';
topdir = 'G:\.shortcut-targets-by-id\1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY\miniscope_v4\Mov_Grid\D_shocks_12mov\camkII_mpfc\2023_09_26\15_51_39\My_WebCam\';
vids = [0:21];
msTSFile = [topdir 'timeStamps.csv'];
[beh] = make_ms_struct(topdir, msTSFile);
px = []; py = [];
figure(1);clf;
set(gcf, 'Position', [25   528   759   230])

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
    for ind = 1:step:frames2read
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
Yrbg = movmean(Yr, 5, 3);
Ygbg = movmean(Yg, 5, 3); % 30
Ybbg = movmean(Yb, 5, 3);
% for fnum = 5000:4:lastf; figure(1); im = cat(3, Yr(:,:,fnum), Yg(:,:,fnum), Yb(:,:,fnum)); bg = cat(3, Yrbg(:,:,fnum), Ygbg(:,:,fnum), Ybbg(:,:,fnum)); image([im; single(im) - bg]); drawnow; end
%%
% MOUSE PARAMS
mouse_sz = 8; % radius

% MOVING GRID PARAMS
minbararea = 60;
bar_thresh = 230;
grid_sz = [80 40]; % h x w

% MS LED PARAMS

led_sizemin = 5;
led_sizemax = 30;
led_valmin = 200;
minmousearea = 100;
searchsize = 50;

led_sz =3; % radius

    


LED_se = strel('disk', led_sz);
smallgauss = gausswin(led_sz*2 - 1, 1); smallgauss  = smallgauss*smallgauss'; smallgauss = smallgauss./sum(smallgauss(:));
mouse_se = strel('disk', mouse_sz);
grid_se = strel('rectangle', grid_sz);
line = strel('rectangle', [3,1]);


mouse_x = NaN(numv*1000,1);
mouse_y = NaN(numv*1000,1);
mouse_ecc = NaN(numv*1000,1);
mouse_size = NaN(numv*1000,1);
barpos = NaN(numv*1000,1);


%
plotting=false;
tic
step = 1;
plotstep = 1:10:numv*1000;
plotstepreport = 1:500:numv*1000;
leddistthresh = 2;
strict_dist_min = 6;
led_dff_thresh = 30;
lastfnum =1;

mouse_gauss = gausswin(2*searchsize + 1, 1.5); mouse_gauss = mouse_gauss*mouse_gauss'; mouse_gauss = mouse_gauss./(sum(mouse_gauss(:)));
% mouse_gauss = mouse_gauss.^2;
px=[]
% px = [32   611];
% py =[40   141];
% lastx=26; lasty=52; lastfnum=-2;
px = [27   608];
py =[29 135];
% vo = VideoWriter("avi_project\mouse_track_demo");
% vo.FrameRate = 30;
% vo.open();
%%
debug = true;
for fnum = 1:step:lastf
    %%
        im = double(cat(3, Yr(:,:,fnum), Yg(:,:,fnum), Yb(:,:,fnum)));
        if isempty(px)
            figure;
            image(uint8(im))
            [px, py] = ginput(2);
            px = round([min(px) max(px)]);
            py = round([min(py) max(py)]);
        end
        bg = double(cat(3, Yrbg(:,:,fnum), Ygbg(:,:,fnum), Ybbg(:,:,fnum)));
        led_dff = mean(im - bg, 3);

        imr = im(:,:,1);
        % img = im(:,:,2);
        % imb = im(:,:,3);
        % bgr = bg(:,:,1);
        % bgg = bg(:,:,2);
        % bgb = bg(:,:,3);
        % % imr = conv2(imr, hsmall, 'same');
        imbar = mean(double(im), 3)>bar_thresh;
        imbar = imopen(imbar, LED_se);
        imbar = imclose(imbar, grid_se);
        imbar = imbar(py(1):py(2), px(1):px(2));
        a = regionprops(imbar, {'Area', 'Eccentricity', 'Centroid'});
        if ~isempty(a)
            bigind = find([a.Area] == max([a.Area]), 1);
            if a(bigind).Area >= minbararea
                barpos(fnum) = a(bigind).Centroid(1);
            end
        end

        immouse = imclose(imr, 1*LED_se.Neighborhood);
        immouse = imopen(immouse, 1*mouse_se.Neighborhood);
        immouse = conv2(immouse, smallgauss, 'same');
        rdiffmouse = -immouse; % imr - immouse;
        rdiff_mouse = rdiffmouse >= quantile(rdiffmouse(:), .995) & rdiffmouse<quantile(imr(:), .5);

        rdiffmouse = rdiffmouse(py(1):py(2), px(1):px(2));
        rdiff_mouse = rdiff_mouse(py(1):py(2), px(1):px(2));
        a = regionprops(rdiff_mouse, {'Area', 'Eccentricity', 'Centroid'});
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
            end
        end
        if any(ismember(fnum, plotstepreport))
            disp(fnum)
        end
        if debug && any(ismember(fnum, plotstep))
            %%
            figure(99); clf
            subplot_tight(2,1,1)
            imagesc(rdiffmouse, [-255 0]) ; hold on
            plot(mouse_x(fnum), mouse_y(fnum), 'ko')
            axis image off
            subplot_tight(2,1,2)
            cla
            image(cat(3, imbar, rdiff_mouse, rdiff_mouse));
            hold on
            scatter(mouse_x(fnum), mouse_y(fnum), 'bx')
            plot([barpos(fnum) barpos(fnum)], [0 size(imbar,1)], 'm--', 'LineWidth', 3)
            axis off image
            title(fnum)
            drawnow
        end
end

end

%% interoplate bar an mouse pos
x = mouse_x(1:lastf);
y = mouse_y(1:lastf);

naninds = isnan(x)& isnan(y);
x(naninds)= interp1(find(~naninds), x(~naninds), find(naninds), 'linear');
y(naninds)= interp1(find(~naninds), y(~naninds), find(naninds), 'linear');
beh_sooth_sz = 7; halfsz = ceil(beh_sooth_sz/2);
smallgauss1d = gausswin(beh_sooth_sz); smallgauss1d = smallgauss1d./sum(smallgauss1d(:)); %%%%%%%%%%
xs = conv(x, smallgauss1d, 'same');
ys = conv(y, smallgauss1d, 'same');
xs(1:halfsz) = x(1:halfsz); xs(end-halfsz:end) = x(end-halfsz:end); 
ys(1:halfsz) = y(1:halfsz); ys(end-halfsz:end) = y(end-halfsz:end); 
speed = speed_calc_gb(xs,ys); 
%%
plotting=true;
tic
step = 1;
plotstep = 1:1000:numv*1000;
leddistthresh = 10;
strict_dist_min = 6;
led_dff_thresh = 1e3;%;
lastfnum =1;

minbararea = 60;
led_sizemin = 4;
led_sizemax = 75;
led_valmin = 200;
minmousearea = 100;
searchsize = 50;
% minLEDval = 150;
minLEDval = 200;

led_x = NaN(numv*1000,1);
led_y = NaN(numv*1000,1);
led_meanval = NaN(numv*1000,1);
led_maxval = NaN(numv*1000,1);
led_dist = NaN(numv*1000,1);
led_size = NaN(numv*1000,1);
headangle = NaN(numv*1000,1);
smallgauss = gausswin(5); smallgauss = smallgauss*smallgauss'; smallgauss = smallgauss./sum(smallgauss(:));

lastx=[]
% px = [27   608];
% py =[29 135];
% lastx=29; lasty=25; lastfnum=-2;
debug = false;
futurewindow = 15;
posdistthresh= 100;
hsmall = strel('disk', 25);
pathconv = hsmall.Neighborhood/sum(hsmall.Neighborhood(:));
hsmall = strel('disk', 9);
ledconv = hsmall.Neighborhood./sum(hsmall.Neighborhood(:));
cmap = lbmap(256, 'RedBlue');
for fnum =  1:step:lastf % 17365:step:17537
    %%
    im = double(cat(3, Yr(:,:,fnum), Yg(:,:,fnum), Yb(:,:,fnum)));
    %
    bg = double(cat(3, Yrbg(:,:,fnum), Ygbg(:,:,fnum), Ybbg(:,:,fnum)));
    % led_dff = mean(im, 3);
    % led_dff = mean(im - bg, 3);
    imr = im(:,:,1);
    % img = im(:,:,2);
    imb = im(:,:,3);
    % bgr = bg(:,:,1);
    % bgg = bg(:,:,2);
    % bgb = bg(:,:,3);
    meanim = mean(double(im), 3);
    imbar = meanim >bar_thresh;
    imbar = imopen(imbar, LED_se);
    imbar = imclose(imbar, grid_se);
    imbar = conv2(imbar, ledconv, 'same')>0;
    % imbar = imbar(py(1):py(2), px(1):px(2));

    if isempty(lastx)
        figure(1); clf;
        im2 = uint8(im(py(1):py(2), px(1):px(2),:));
        image(im2)
        [lastx, lasty] = ginput(1);
        lastx = round(lastx);
        lasty = round(lasty);
        lastfnum = fnum-3;
        % ledVal = mean(mean(mean(im2(lasty-2:lasty+2, lastx-2:lastx+2, :))));
    end

    futureinds = fnum:fnum+futurewindow;
    futureinds = futureinds(futureinds<=lastf & futureinds>0);
    % spdinds = fnum-2:fnum+2;
    % spdinds = spdinds(spdinds<=lastf & spdinds>0);
    im2 = im(py(1):py(2), px(1):px(2),1);

    % Use the next x indices to create a likeliehood map of where the LED might be
    posdist = sqrt( (xs(fnum) - xs(fnum:end)).^2 + (ys(fnum) - ys(fnum:end)).^2 );
    % intspd = mean(speed(spdinds));
    % frame_distthresh = leddistthresh*(fnum - lastfnum) + strict_dist_min + floor(intspd);

    farinds = fnum:fnum+find(posdist>posdistthresh, 1, 'first');
    farinds = farinds(farinds<=lastf);
    farmap = zeros(size(im2));
    futuremap = zeros(size(im2));
    pathmap = zeros(size(im2));
    for j = length(farinds):-1:1
        farmap(round(ys(farinds(j))), round(xs(farinds(j)))) =  length(farinds) - j + 1;
    end
    for j = length(futureinds):-1:1
        futuremap(round(ys(futureinds(j))), round(xs(futureinds(j)))) = length(futureinds) - j + 1;
    end
    futuremap = conv2(futuremap, pathconv, 'same');
    farmap = conv2(farmap, pathconv, 'same');
    futuremap = futuremap./max(futuremap(:));
    farmap = farmap./max(farmap(:));
    validmap = conv2(futuremap>0 & farmap>0, pathconv, 'same')>0;
    pathmap(validmap) = farmap(validmap) + futuremap(validmap);
    pathmap = conv2(pathmap, pathconv, 'same');
    ledmap = zeros(size(im2));
    if ~isempty(lastx)
        ledmap(round(lasty), round(lastx)) = 10;
    end
    ledmap = conv2(ledmap, ledconv, 'same')>0;
    ledmap = conv2(ledmap, ledconv, 'same');
    ledmap = 1*ledmap./max(ledmap(:));


    pathmap = pathmap+ledmap + ledmap;
    validmap = pathmap>0;


    %% imagesc(led_dff)
    showLED = false;

    imavbg = mean(im-bg, 3);
    imav = mean(im, 3);
    imLED = imav.*imb;
    imLED = conv2(imLED, ones(3,3)./9, 'same'); % small uniform smoothing
    imLED = imLED - conv2(imLED, ones(11,3)./(11*3), 'same'); % remove vertical bar noise
    d = imLED > 5000 & ~imbar;
    validLED=d(py(1):py(2), px(1):px(2)).*validmap;
    validLED = validLED>0;
    imLED = imLED(py(1):py(2), px(1):px(2));
    if showLED==true && debug
        if isbg
        imav2 = normalize_matrix(imavbg(py(1):py(2), px(1):px(2)));
        else
        imav2 = normalize_matrix(imav(py(1):py(2), px(1):px(2)));
        end
        imLED2 = normalize_matrix(imLED);
        figure(111); clf; set(gcf, 'Position', [163   412   560   189])
        imagesc([imav2; imLED2; validLED]);
        title(['bg is ' num2str(1*isbg)])
    end
        imav=imav(py(1):py(2), px(1):px(2)).*validmap;
        imb=imb(py(1):py(2), px(1):px(2)).*validmap;
        imr=imr(py(1):py(2), px(1):px(2)).*validmap;
    a = regionprops(validLED>0, {'Area', 'Centroid', 'PixelIdxList'});
    %%
    if ~isempty(a)
        % closest
        c = cat(1, a.Centroid);
        sz = cat(1, a.Area);
        if ~isempty(lastx) && (fnum - lastfnum)<=10 %% fnum>step+1 &&
            ax = lastx - c(:,1);
            ay = lasty - c(:,2);
            prevdists = sqrt( (ax.^2) + (ay.^2) );
            prevdists(prevdists>posdistthresh) = NaN;
        else
            prevdists = ones(length(a),1);
            % lastx=[]; lasty=[];
        end
        
        ax = xs(fnum) - c(:,1);
        ay = ys(fnum) - c(:,2);

        dists = sqrt( (ax.^2) + (ay.^2) );
        % small angular difference
        x10 = c(:,1) - xs(fnum);
        y10 = ys(fnum) - c(:,2);
        ang =  rad2deg( atan2(y10, x10) );
        lastang = headangle(find(~isnan(headangle), 1, 'last'));
        if ~isempty(lastang)
            angdiff=abs(ang-lastang);
        else
            angdiff = dists*0 + 1;
        end

        prevdists(dists>posdistthresh/2) = NaN; % was 50
        dists_heuristic = prevdists.*angdiff;
        % figure(112); clf
        % imshow(imav2, [0 1]);
        % hold on
        % scatter(xs(fnum), ys(fnum), 'co');
        % scatter(c(:,1), c(:,2), sz, 'mo')
        % drawnow
        if any(~isnan(dists_heuristic))
            [dists_heuristic, ind_ord] = sort(dists_heuristic, 'ascend');
            ind_ord = ind_ord(~isnan(dists_heuristic));
            aa=1; ledfound=false;
            while ledfound==false && (aa<= length(ind_ord))% for aa = 1:length(a)
                idx = ind_ord(aa);
                thissize = a(idx).Area;
                maxVal = max(max(imb( a(idx).PixelIdxList )));
                avVal = mean(mean(imav( a(idx).PixelIdxList )));
                % scatter(c(idx,1), c(idx,2), sz, 'mx')
                if thissize >= led_sizemin && thissize <= led_sizemax && maxVal>=minLEDval  % && avVal >= led_valmin % imr(a(idx).Image)
                % scatter(c(idx,1), c(idx,2), sz, 'kx')
                    led_x(fnum) = a(idx).Centroid(1);
                    led_y(fnum) = a(idx).Centroid(2);
                    lastx = led_x(fnum);
                    lasty = led_y(fnum); lastfnum = fnum;
                    led_meanval(fnum) = avVal;
                    led_maxval(fnum) = maxVal;
                    led_dist(fnum) = dists(idx);
                    led_size(fnum) = a(idx).Area;
                    ledfound=true;
                    headangle(fnum) = ang(idx);
                end
                aa=aa+1;
            end
            aa = aa-1;
        end
    end
    if debug && any(ismember(fnum, plotstep))
        %%
        figure(99); clf
        set(gcf, 'Position', [136    80   560   319])
        subplot_tight(3,1,1)
        imagesc(imLED2, [0 1]); axis image; colormap(cmap)
        subplot_tight(3,1,2)
        imagesc(validLED + pathmap, [0 4]); axis image
        hold on
        scatter(xs(fnum), ys(fnum), 100, 'xk');            scatter(xs(fnum), ys(fnum), 100, 'ko')
        if ~isempty(a)
            scatter(c(:,1), c(:,2),2*sz, 'kx')
        end
        % if ~isempty(lastx)
        %     pos = [lastx - frame_distthresh lasty - frame_distthresh frame_distthresh*2 frame_distthresh*2];
        % end
        % rectangle('Position', pos, 'EdgeColor', 'w', 'FaceColor', 'none')
        subplot_tight(3,1,3)
        cla
        imbar2 = imbar(py(1):py(2), px(1):px(2));

        image(cat(3, imbar2+ validLED, validLED, imbar2*0));
        hold on
        scatter(xs(fnum), ys(fnum), 'wx')
        scatter(xs(fnum), ys(fnum), 'wo')
        plot(xs(futureinds), ys(futureinds), 'w-')
        plot(xs(farinds), ys(farinds), 'm-')
        plot([barpos(fnum) barpos(fnum)], [0 size(imbar2,1)], 'm--', 'LineWidth', 3)
        axis off image
        drawnow
    end

    if plotting==true && any(ismember(fnum, plotstep))
        %%
        im2 = uint8(im(py(1):py(2), px(1):px(2),:));
        figure(1);clf;
        set(gcf, 'Position', [40   645   666   203])
        subplot_tight(2,3,1:2)
        image(im2); axis image; title(['Frame  ' num2str(fnum)]); axis image off
        subplot_tight(2,3,4:5)
        imL = im2; imL = imL+.8*(255-imL);
        image(imL); axis image off; %title(num2str(fnum))
        hold on
        previnds = fnum-30:fnum; previnds = previnds(previnds>0);
        plot(xs(1:fnum), ys(1:fnum), '-', 'Color', [.8 .7 .9])
        plot(xs(previnds), ys(previnds), '-', 'Color', [.2 .2 .2])
        plot(led_x(previnds), led_y(previnds), 'r-')
        plot([barpos(fnum) barpos(fnum)], [0 size(im2,1)], 'm--')

        subplot_tight(1,3, 3, [.05 .125])
        isub = im2;
        imL = im2; imL = imL+.2*(255-imL);
        image(imL); hold on; axis image
        % plot([0 x10]-searchsize+size(isub,2), [0 -y10]-searchsize+size(isub,1), 'ro-');
        scatter([xs(fnum), led_x(fnum)], [ys(fnum), led_y(fnum)], 'ko', 'MarkerFaceColor', 'k')
        plot([xs(fnum), led_x(fnum)], [ys(fnum), led_y(fnum)], 'r-', 'LineWidth', 2)
        % plot([0 x20]-searchsize+size(isub,2), [0 y20]-searchsize+size(isub,1), 'ko-');
        ha = headangle(fnum);
        if ~isnan(ha)
            title(sprintf('Head Ang %4.fdeg', ha))
        else
            title(sprintf('Head Ang      %sdeg', ' '))
        end
        axis image off % ([0 searchsize*2 0 searchsize*2])
        lastgood = fnum;%find(~isnan(xs),1,'last');
        ww = round(xs(lastgood))-searchsize:round(xs(lastgood))+searchsize;
        hh = round(ys(lastgood))-searchsize:round(ys(lastgood))+searchsize;
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
lx = led_x(1:lastf);
ly = led_y(1:lastf);
lh = headangle(1:lastf);

naninds = isnan(lx)& isnan(ly);
lx(naninds)= interp1(find(~naninds), lx(~naninds), find(naninds), 'linear');
ly(naninds)= interp1(find(~naninds), ly(~naninds), find(naninds), 'linear');
naninds = isnan(lh);
lh(naninds)= interp1(find(~naninds), lh(~naninds), find(naninds), 'linear');

smallgauss = gausswin(7); smallgauss = smallgauss./sum(smallgauss(:));
lxs = conv(lx, smallgauss, 'same');
lys = conv(ly, smallgauss, 'same');
lhs = conv(lh, smallgauss, 'same');
lxs(1:4) = lx(1:4); lxs(end-4:end) = lx(end-4:end); 
lys(1:4) = ly(1:4); lys(end-4:end) = ly(end-4:end); 
lhs(1:4) = lh(1:4); lhs(end-4:end) = lh(end-4:end); 
% speed = speed_calc_gb(lxs,lys); 
%%
if false
vo = VideoWriter('test.avi');
vo.open()
searchsize = 75;
redfade = [1 0 0] + (ones(31,1).*linspace(1,0,31)')*[0 1 1];
figure(1);clf;
set(gcf, 'Position', [40   645   666   203], 'Color', 'k')
startf = 7100;
for fnum = startf:4:lastf
        %%
        im = double(cat(3, Yr(:,:,fnum), Yg(:,:,fnum), Yb(:,:,fnum)));
        im2 = uint8(im(py(1):py(2), px(1):px(2),:));
        figure(1);clf;
        set(gcf, 'Position', [40   645   666   203])
        subplot_tight(2,3,1:2, [.1 0])
        image(im2); axis image; 
        title(sprintf('Frame- %d   Time- %5.0f sec', fnum-startf, (beh.timestamps(fnum)-beh.timestamps(startf))./1000 ), 'Color', 'w'); % ['Frame  ' num2str(fnum)]); 
        axis image off
        subplot_tight(2,3,4:5, [0.1 0])
        imL = im2; imL = imL+.8*(255-imL);
        image(imL); axis image off; %title(num2str(fnum))
        hold on
        previnds = fnum-30:fnum; 
        clr = redfade(previnds>0,:);
        previnds = previnds(previnds>0);
        plot(xs(startf:fnum), ys(startf:fnum), '-', 'Color', [.8 .7 .9])
        plot(xs(previnds), ys(previnds), '-', 'Color', [.2 .2 .2])
        scatter(lxs(previnds), lys(previnds), 10, clr, '.')
        plot([barpos(fnum) barpos(fnum)], [0 size(im2,1)], 'm--')

        subplot_tight(1,3, 3, [.05 .125])
        % isub = im2;
        imL = uint8(im); imL = imL+.2*(255-imL);
        image(imL); hold on; axis image
        % plot([0 x10]-searchsize+size(isub,2), [0 -y10]-searchsize+size(isub,1), 'ro-');
        tx = xs(fnum) + px(1); 
        tlx = lxs(fnum) + px(1); 
        ty = ys(fnum) + py(1);
        tly = lys(fnum) + py(1);
        plot([tx, tlx], [ty tly], '-', 'LineWidth', 2, 'Color', [.7 .7 .7])
        % scatter([tx, tlx], [ty tly], 'ko', 'MarkerFaceColor', 'k')

        % plot(xs(1:fnum)+px(1), ys(1:fnum)+py(1), '-', 'Color', [.8 .7 .9])
        plot(xs(previnds)+px(1), ys(previnds)+py(1), '-', 'Color', [.0 .0 .0], 'LineWidth', 2)
        scatter(lxs(previnds)+px(1), lys(previnds)+py(1), 30, clr, '.')
        
        % plot([0 x20]-searchsize+size(isub,2), [0 y20]-searchsize+size(isub,1), 'ko-');
        ha = lhs(fnum);
        if ~isnan(ha)
            title(sprintf('Head Ang %4.f deg', ha), 'Color', 'w')
        else
            title(sprintf('Head Ang      %s deg', ' '), 'Color', 'w')
        end
        axis image off % ([0 searchsize*2 0 searchsize*2])
        lastgood = fnum;%find(~isnan(xs),1,'last');
        ww = round(mean([tx,tlx]))-searchsize:round(mean([tx,tlx]))+searchsize;
        hh = round(mean([ty,tly]))-searchsize:round(mean([ty,tly]))+searchsize;
        wbad = ww>1 & ww <= size(imL,2);
        hbad = hh>1 & hh <= size(imL,1);
        ww = ww(wbad);
        hh = hh(hbad);
        axis([ww(1), ww(end) hh(1) hh(end)])
        % temp = getframe(gca);
        drawnow
        temp = getframe(gcf);
        frameout = temp.cdata;
        vo.writeVideo(frameout)
end
vo.close();
end
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