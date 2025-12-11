%%
clear
% fname = "C:\Users\gjb326\Desktop\TRACKER DOCS\Tracker Videos\NAPA_24510_IL7_TrackerVideo.avi";
% dname = "C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL7.mat";
fname = "C:\Users\gjb326\Desktop\TRACKER DOCS\Tracker Videos\NAPA_24510_IL4_TrackerVideo.avi";
dname = "C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL4.mat";
d=load(dname);
v = VideoReader(fname);
nf = ceil(v.Duration*v.FrameRate);

% nf=2000;
% nf = ceil(timetoget*v.FrameRate);
ds = 10;
Y = uint8(false(v.Height, v.Width, nf));
vt = zeros(nf,1);

idx = 0;
badcount = 0;
tic
% v.open
while v.hasFrame % && v.CurrentTime<timetoget
    im = rgb2gray(v.readFrame);
    if any(im(:)>0)
        idx = idx+1;
        Y(:,:,idx) = im;
        vt(idx) = v.CurrentTime;
    else
        badcount = badcount+1;
    end

end
badcount;
Y = Y(:,:,1:idx);
vt = vt(1:idx);
toc
%%
d.params=d.params_sub;
px = [119.5383  113.5642  109.8946  111.7294  117.6926  158.0587  155.7651  157.1413  158.9761  162.6457 119.5383  ];
py = [297.5677  276.8502  247.7317  218.6132  186.5082  203.6806  223.0930  241.7587  256.6912  274.6103 297.5677  ];

nframes=size(Y,3);
xfs = 1:nframes;
rx = d.room.x;
ry = d.room.y;
ax = d.arena.x;
ay = d.arena.y;
behavtime = d.room.timestamps/1000;
dt = round(1/median(abs(diff(vt))));

% x = d.room.x(1:4:nf);
% y = d.room.y(1:4:nf);
% (x-params.arena_center(1))./params.pixpercm
ad1 = angdiff(d.arena.pol_theta, d.room.pol_theta);
ad = unwrap(ad1);
dad = abs(diff(ad)); 
ad(dad>.01)=NaN;
nanind = find(isnan(ad));
ad(nanind) = interp1(find(~isnan(ad)), ad(~isnan(ad)), nanind, 'linear');
s=300;
ad(300/2:end-s/2) = conv(ad, ones(s,1)/s,'valid');
ad = mod(ad+pi, -2*pi)+pi;

%
rx = 1.35*rx*d.params.pixpercm/d.room.DAT_fileinfo.ArenaDiameter_m + d.room.DAT_fileinfo.ArenaCenterXY(1) + 190;
ry = -1.35*ry*d.params.pixpercm/d.room.DAT_fileinfo.ArenaDiameter_m + d.room.DAT_fileinfo.ArenaCenterXY(2) + 110;
ax = 1.35*ax*d.params.pixpercm/d.arena.DAT_fileinfo.ArenaDiameter_m + d.arena.DAT_fileinfo.ArenaCenterXY(1) + 190;
ay = -1.35*ay*d.params.pixpercm/d.arena.DAT_fileinfo.ArenaDiameter_m + d.arena.DAT_fileinfo.ArenaCenterXY(2) + 110;
showcropx = [81, 560];
showcropy = [1, 480];

rx = rx-showcropx(1); ry = ry-showcropy(1);
ax = ax-showcropx(1); ay = ay-showcropy(1);
% ztop=[225, 22; 200, 100];
% zbot=[225, 350; 200, 100];
% % close all
% figure(10); clf; set(gcf, 'Color', 'k', 'Position', round([-1600 80 561*1.75 517*1.75]/2)*2); clf;% imshow(max(Y,[],3))
% set(gcf,'renderer','painters')
% colormap plasma
% tail = round(dt*3);
% c_tail = magma(tail);
% c = bone(length(rx));
% c = reshape(c, [size(c,1),1,3]);
% c(end,:)=NaN;
% c(1,:)=NaN;
% c_tail = reshape(c_tail, [size(c_tail,1),1,3]);
% c_tail(end,:)=NaN;
% c_tail(1,:)=NaN;
%%
for i= 1:1:nframes
    %%
    idx = i; % xfs(i);
    clf
    im = double(squeeze(Y(:,:,idx)));
    vidt = vt(idx);
    behav_idx = find(min(abs(behavtime-vidt)) == (abs(behavtime-vidt)), 1);
    if ~any(im(:)>0)
        im = lastim;
    else
        lastim = im;
    end
    im_room = im(showcropy(1):showcropy(2), showcropx(1):showcropx(2));
    imr = imrotate(im, -rad2deg(ad(idx)), 'nearest', 'crop');
    im_arena = imr(showcropy(1):showcropy(2), showcropx(1):showcropx(2));
    if i==1
        YA = uint8(false(size(im_arena,1), size(im_arena,2), nframes));
    end
    YA(:,:,i) = im_arena;
 
%     subplot_tight(1,2,1, [.01, .01])
%     imshow(im_room, [0 250]); 
%     title('\color[rgb]{1 .4 .4}Room Reference Frame', 'interpreter', 'tex')
%     subplot_tight(1,2,2, [.01, .01])
%     imshow(im_arena, [0 250]); 
%     drawnow
end

%%
bg = double(std(single(YA),[],3));
bg2 = double(mean(YA,3));
bg2(bg2>75) = 75; % IR lights in the center make detecting the rat very bad, so clip the values a bit
bg(bg2>20) = 20;
% Create a grid of coordinates
[xt, yt] = meshgrid(1:size(im_arena,2), 1:size(im_arena, 1));

% Center of the image
center_x = size(im_arena, 2) / 2;
center_y = size(im_arena, 1) / 2;

% Shift coordinates to the center
xt = xt - center_x;
yt = yt - center_y;

% Convert to polar coordinates
rho = sqrt(xt.^2 + yt.^2);
theta = atan2d(yt, xt); % Angle in degrees

% Normalize theta to [0, 360)
% theta(theta < 0) = theta(theta < 0) + 360;

% Create the binary mask
maxextent = 210;
mask = ((rho > 165) & (rho < maxextent)) | (abs(xt)<20 & rho<maxextent);

% Display the mask
clf
imshow(mask.*bg2, []);




%
behavtime = d.room.timestamps/1000;

figure(1); clf; colormap viridis
ratsize = 21;
kern = gausswin(ratsize); kern = kern.*(kern');
kern = kern./sum(kern(:));
xi = NaN(nframes,1); yi = NaN(nframes,1);
%%
lastx = xi(1); lasty = yi(1); lastind=1;

fs = 1:10:nframes;

plotting = true;
for ii= 1:length(fs)
    %%
    i = fs(ii);
    vidt = vt(i);
    behav_idx = find(min(abs(behavtime-vidt)) == (abs(behavtime-vidt)), 1);

    clf
    previ = i-10:i; previ = previ(previ>0);
    im = double(squeeze(YA(:,:,i)));
    %     im2 = abs(im-mean(double(squeeze(YA(:,:,previ))), 3));
    %     im2_z = conv2(im2, kern, 'same').*mask;
    %     im2_mask = im2_z >= .5*max(im2_z(:));
    %     subplot_tight(1,3,2, [.01, .01])
    %     imshow(im2_mask, [-5 5]); 
    
    im_z = (im - bg2)./bg;
    im_z = conv2(im_z, kern, 'same');
    im_mask = im_z.*mask;
    im_mask = im_mask <= .8*min(im_mask(:));
    if i-lastind > 100
        lastx = NaN;
        lasty = NaN;
    end
    
    im_mask = conv2(im_mask, kern, 'same')>.25;
    max_dist = 100;
    min_size = 30;
    [xi(i), yi(i), largestComponentMask] = find_largest_BW(im_mask, lastx, lasty, max_dist, min_size);
    if any(largestComponentMask(:)) >0
        lastx = xi(i);
        lasty = yi(i);
        lastind=i;
    end
    if plotting==true
        subplot_tight(1,2,1, [.01, .01]);
        imshow(im+im_mask.*250, [0 250]);  hold on
        plot(ax(1:behav_idx), ay(1:behav_idx), 'k')
        plot(xi(fs(1:ii)), yi(fs(1:ii)), 'r-')
        set(gca, 'YDir', 'reverse')
        
        subplot_tight(1,2,2, [.01, .01]); hold on
        imagesc(im_mask+10*largestComponentMask, [-2 2]); axis image off
        %     colorbar
        plot(xi(fs(1:ii)), yi(fs(1:ii)), 'r-')
        set(gca, 'YDir', 'reverse')
        title(i)
        drawnow
        pause(.2)
    end
end
%%
figure; hold on; 
plot3(ax, ay, 1:length(ax))
plot3(xi, yi, 1:length(xi))

%%


save_vid = false;
if save_vid
    vv = VideoWriter('C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\test2.avi');
    vv.Quality=80;
    vv.FrameRate = 30;
    vv.open();
end
%
function [xi, yi, largestComponentMask] = find_largest_BW(im_mask, x, y, max_dist, min_size)

connectedComponents = bwlabel(im_mask);
% Measure properties of the connected components
stats = regionprops(logical(connectedComponents), 'Area', 'Centroid');
largestComponentMask = im_mask*0;
xi = NaN;
yi = NaN;
% Find the largest connected component
if ~isempty(stats)
    % Extract areas
    stats = stats([stats.Area]>min_size);
    if ~isempty(stats)
        areas = [stats.Area];
        if isnan(x)
            % Find the largest area
            [~, Idx] = max(areas);
        else
            % Find the closest area
            dxy = NaN(length(stats),1);
            for i = 1:length(stats)
                dxy(i) = sqrt( (x-stats(i).Centroid(1))^2 + (y-stats(i).Centroid(2))^2 );
            end
            dxy(dxy>max_dist) = NaN;
            if any(~isnan(dxy))
                Idx = find(dxy == nanmin(dxy));
            else
                Idx = NaN;
            end
        end
        % Get the centroid of the largest connected component
        if ~isnan(Idx)
            largestCentroid = stats(Idx).Centroid;
            xi = largestCentroid(1);
            yi = largestCentroid(2);
            % Display results
            %         fprintf('Centroid of the largest connected component: (%.2f, %.2f)\n', ...
            %             largestCentroid(1), largestCentroid(2));
            % Optional: Visualize the largest connected component
            largestComponentMask = connectedComponents == Idx;
        end
    end
end

end