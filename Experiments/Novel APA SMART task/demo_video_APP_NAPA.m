%%
clear
fname = "C:\Users\gjb326\Desktop\TRACKER DOCS\Tracker Videos\NAPA_24510_IL5_TrackerVideo.avi";
dname = "C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL5.mat";
% fname = "C:\Users\gjb326\Desktop\TRACKER DOCS\Tracker Videos\NAPA_24511_IL9_TrackerVideo.avi";
% dname = "C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24511_IL9.mat";
d=load(dname);
v = VideoReader(fname);
% nf = ceil(v.Duration*v.FrameRate);

% nf=2000;
timetoget =60*30;
nf = round(timetoget*v.FrameRate);
ds = 10;
Y = uint8(false(v.Height, v.Width, nf));
vt = zeros(nf,1);

idx = 0;
badcount = 0;
tic

[Y, vt, badcount] = load_tracker_avi(fname, []);
toc

% % v.open
% while v.hasFrame && v.CurrentTime<timetoget
%     im = rgb2gray(v.readFrame);
%     if any(im(:)>0)
%         idx = idx+1;
%         Y(:,:,idx) = im;
%         vt(idx) = v.CurrentTime;
%     else
%         badcount = badcount+1;
%     end
% 
% end
% badcount
% Y = Y(:,:,1:idx);

% bg=
%%
px = [93.9813   91.3854   89.1604   86.9354   85.8229   85.4521   85.4521   86.9354   89.5313   92.1271   94.3521   86.9354   46.8854   43.9188   42.4354   40.9521   39.4688   37.6146   36.5021   35.7604   34.2771   34.4646   34.7604   36.1313   37.6146   38.7271   40.9521   42.4354   43.5479   44.2896   46.1438   46.8854];
py = [293.2146  282.8312  274.6729  263.1771  249.0854  239.8146  232.3979  217.5646  202.7312  196.0562  187.5271  183.4479  169.7271  178.9979  186.0437  191.9771  199.0229  207.5521  217.1937  224.6104  233.1396  243.8937  253.4229  265.0313  274.3021  280.2354  290.6187  295.8104  300.2604  303.2271  308.4188  309.9021];
% px = [119.5383  113.5642  109.8946  111.7294  117.6926  158.0587  155.7651  157.1413  158.9761  162.6457 119.5383  ];
% py = [297.5677  276.8502  247.7317  218.6132  186.5082  203.6806  223.0930  241.7587  256.6912  274.6103 297.5677  ];


% correct_runav = 0*d.zonestruct.choices.isCorrect;
% for ii = 1:length(d.zonestruct.choices.isCorrect)
%     correct_runav(ii) = sum(d.zonestruct.choices.isCorrect(1:ii)==1)/ii;
% end
% left_runav = 0*d.zonestruct.choices.isLeft;
% for ii = 1:length(d.zonestruct.choices.isLeft)
%     left_runav(ii) = sum(d.zonestruct.choices.isLeft(1:ii)==1)/ii;
% end
all_previous = true;
[outvars] = SMART_useful_figs(d.zonestruct, 0, 60, all_previous, false);
correct_runav = outvars.correct_runav_time;
left_runav = outvars.left_runav_time;
runav_time = outvars.correct_runav_time_xv;

nframes=size(Y,3);
xfs = 1:nframes;
rx = d.room.x;
ry = d.room.y;
ax = d.arena.x;
ay = d.arena.y;
behavtime = d.room.timestamps/1000;
dt = round(1/median(abs(diff(vt))));

% ad1 = angdiff(d.arena.pol_theta, d.room.pol_theta);
% ad = unwrap(ad1);
% dad = abs(diff(ad)); 
% ad(dad>.01)=NaN;
% nanind = find(isnan(ad));
% ad(nanind) = interp1(find(~isnan(ad)), ad(~isnan(ad)), nanind, 'linear');

s=300;
ad = d.arena.angular_offset;
ad(s/2:end-s/2) = conv(ad, ones(s,1)/s,'valid');
ad = mod(ad+pi, -2*pi)+pi;
%
scale = 1.2;
rx = scale*rx*d.params_sub.pixpercm/d.room.DAT_fileinfo.ArenaDiameter_m + d.room.DAT_fileinfo.ArenaCenterXY(1) + 190;
ry = -scale*ry*d.params_sub.pixpercm/d.room.DAT_fileinfo.ArenaDiameter_m + d.room.DAT_fileinfo.ArenaCenterXY(2) + 110;
ax = scale*ax*d.params_sub.pixpercm/d.arena.DAT_fileinfo.ArenaDiameter_m + d.arena.DAT_fileinfo.ArenaCenterXY(1) + 190;
ay = -scale*ay*d.params_sub.pixpercm/d.arena.DAT_fileinfo.ArenaDiameter_m + d.arena.DAT_fileinfo.ArenaCenterXY(2) + 110;
showcropx = [81, 560];
showcropy = [1, 480];

% rx = rx-showcropx(1); ry = ry-showcropy(1);
% ax = ax-showcropx(1); ay = ay-showcropy(1);
rx = rx-showcropx(1); ry = ry-showcropy(1);
ax = ax-showcropx(1); ay = ay-showcropy(1);

ztop=[225, 22; 200, 100];
zbot=[225, 350; 200, 100];
% close all
figure(10); clf; set(gcf, 'Color', 'k', 'Position', round([400 80 561*1.75 517*1.75]/2)*2); clf;% imshow(max(Y,[],3))
set(gcf,'renderer','painters')
colormap plasma
tail = round(dt*3);
c_tail = magma(tail);
c = bone(length(rx));
c = reshape(c, [size(c,1),1,3]);
c(end,:)=NaN;
c(1,:)=NaN;
c_tail = reshape(c_tail, [size(c_tail,1),1,3]);
c_tail(end,:)=NaN;
c_tail(1,:)=NaN;




save_vid = true;
if save_vid
    vv = VideoWriter('C:\Users\gjb326\Desktop\TRACKER DOCS\figures\vids\test3.avi');
    vv.Quality=80;
    vv.FrameRate = 30;
    vv.open();
end
%

frame2show = [dt*2:4:dt*90, dt*90+1:20:dt*60*9.5, dt*60*9.5+1:4:dt*60*11, dt*60*11+1:20:dt*60*26.5];
% frame2show = [2:4:dt*60];
lightgreen = [.4 1 .4];
lightblue  = [.4 .4 1];
lightred   = [1 .4 .4];
for i= 5:1:length(frame2show) % 40000+tail:50:nframes
    %%
    idx = frame2show(i); % xfs(i);
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

    subplot_tight(3,2,[1,3], [.01, .01])
%     imagesc(im_room, [0 250]); 
    imshow(im_room, [0 250]); 
    title('\color[rgb]{1 .4 .4}Room Reference Frame', 'interpreter', 'tex')
    hold on
    reinforce_window = any(abs(behav_idx - d.room.entrance_start_idx)<=15);
    
    marks = d.room.entrance_start_idx(d.room.entrance_start_idx<=behav_idx);
    if reinforce_window
%         patch(px-showcropx(1), py-showcropy(1), 'g', 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0)
        patch(px, py, 'g', 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0)
        coloridx = 2;
    else
%         patch(px-showcropx(1), py-showcropy(1), 'r', 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', .5)
        patch(px, py, 'r', 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', .5)
        coloridx = 1;
    end
    
    imr = imrotate(im, -rad2deg(ad(idx)), 'nearest', 'crop');
    im_arena = imr(showcropy(1):showcropy(2), showcropx(1):showcropx(2));
    im_arena(im_arena==0) = .2*255;
    %     if length(c)== length(rx)
    inds = 1:behav_idx-1;
    patch(rx(inds), ry(inds),c(inds,:,:).*0,'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .05)
    inds_tail = behav_idx-tail:behav_idx-1;
    inds_tail = inds_tail(inds_tail>0 & inds_tail<nframes);
    patch(rx(inds_tail), ry(inds_tail), c_tail(1:length(inds_tail),:,:),'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .5)
    scatter(rx(behav_idx), ry(behav_idx), 30, 'wo', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor','w', 'MarkerFaceAlpha', .9)
    scatter(rx(marks), ry(marks), 20, 'wo', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor', lightgreen, 'MarkerFaceAlpha', .5)
    text(size(im_room,2)*.7, size(im_room,2)*.05, sprintf('Time:  %4d', round(behavtime(behav_idx))), 'Color', 'w', 'FontSize', 20)
    
    axis image off
    
%     [px2_t, px2_r] = cart2pol(px-showcropx(1)-size(im_arena,2)/2, py-showcropy(1)-size(im_arena,1)/2);
    [px2_t, px2_r] = cart2pol(px-size(im_arena,2)/2, py-size(im_arena,1)/2);
    px2_t = mod(px2_t + ad(idx), 2*pi);
    [px2, py2] = pol2cart(px2_t, px2_r);
    px2=px2 + size(im_arena,2)/2;
    py2=py2 + size(im_arena,1)/2;
    subplot_tight(3,2,[2,4], [.01, .01])
%     imagesc(im_arena, [0 250]); axis image off
    imshow(im_arena, [0 250]); axis image off
    title('\color[rgb]{.4 .4 1}Arena Reference Frame')
    hold on

    if reinforce_window
        patch(px2, py2, 'g', 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0)
    else
        patch(px2, py2, 'r', 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', .5)
    end
    
    patch(ax(inds), ay(inds),c(inds,:,:).*0,'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .05)
    patch(ax(inds_tail), ay(inds_tail), c_tail(1:length(inds_tail),:,:), 'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .5)
    scatter(ax(behav_idx), ay(behav_idx), 30, 'o', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor','w', 'MarkerFaceAlpha', .9) 
    scatter(ax(marks), ay(marks), 30, 'o', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor', lightgreen, 'MarkerFaceAlpha', .5)
%     text(20, 20, sprintf('%2.2f', ad(idx)/pi), 'Color', 'w')
    yy1=.9; yy2 = .42;
    if abs(ad(idx)/pi) < .5 % left is correct
        
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.691, .619],  'Y', [yy1, yy1], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.691, .619],  'Y', [yy2, yy2], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [yy1, yy1], 'Color', lightgreen, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [yy2, yy2], 'Color', lightgreen, 'LineWidth', 3, 'HeadWidth', 10);
    else
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.799, .871],  'Y', [yy1, yy1], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.799, .871],  'Y', [yy2, yy2], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [yy1, yy1], 'Color', lightgreen, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [yy2, yy2], 'Color', lightgreen, 'LineWidth', 3, 'HeadWidth', 10);
    end
    
    subplot_tight(3,4,[10 11], [.07 .01]);    hold on

    bt = behavtime(behav_idx)/60;
    plot([0 30], [.5 .5], 'k:', 'LineWidth', 2)
    plot([bt bt], [0 2], 'k-', 'LineWidth', 2)
    if bt>=10
        plot([10 10], [0 2], 'r', 'LineWidth', 2)
        text(10.5, .3, 'Rotation on', 'Color', 'r')
    end

    pL = left_runav; % d.zonestruct.metrics.sequence_prob_Left;
    pC = correct_runav; % d.zonestruct.metrics.sequence_prob_Correct;
    pX = runav_time/60; % d.zonestruct.metrics.sequence_prob_times / 60;
    plotidx = 1:ceil(bt/1);
%     plotidx = 1:ceil(bt/5);
    plotidx = plotidx(1:min(max(plotidx), length(pX)));
    plotidx = pX<=bt;
    plot(pX(plotidx), pL(plotidx), '-', 'Color',  lightblue, 'LineWidth', 4.5)
    plot(pX(plotidx), pC(plotidx), '-', 'Color', lightred, 'LineWidth', 3)
    L1 = scatter(pX(find(plotidx==1,1,'last')), pL(find(plotidx==1,1,'last')), 150, 's', 'MarkerFaceColor', lightblue, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
    L2 = scatter(pX(find(plotidx==1,1,'last')), pC(find(plotidx==1,1,'last')), 100, 'o', 'MarkerFaceColor', lightred, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
%     rectangle('Position', [bt 0 50 1], 'FaceColor', 'w', 'EdgeColor', 'none')
%     axis([-1, 27, .4, 1.1])
    axis([-1, 27, .2, .8])
    ylabel('Probability')
    xlabel('Time (minutes)')
    legend([L1, L2], 'Left turn', 'Correct turn', "Location", 'northwest')
    title('\color[rgb]{1 .4 .4}Room \color[rgb]{1 1 1}vs \color[rgb]{.4 .4 1}Arena \color[rgb]{1 1 1} choice performance', 'FontSize', 12)
    set(gca, 'YColor', 'w', 'XColor', 'w', 'YTick', [.25 .5 .75 1], 'XTick', [0:5:25], 'FontName', 'Arial', 'FontWeight', 'bold')
    drawnow()
    %%
    if save_vid==true
%         temp = getframe(gcf);
        temp = getframe(gcf);
        writeVideo(vv, temp)
    end
end
if save_vid==true
    vv.close()
end
%% Version with two choice ROIs inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10); clf;% imshow(max(Y,[],3))
tail = 120;
c = gray(tail);
% c = gray(length(xfs));
c = reshape(c, [size(c,1),1,3]);
c(end,:)=NaN;
c(1,:)=NaN;
cmap_im = viridis(256);

% c(end,:)=NaN;
% c = reshape(c, [tail,1,3]);
for i = 20000:ds:22000 % 20000+1+tail:ds:nframes
    %%
    idx = xfs(i);
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

    subplot_tight(5,2,[3:2:7], [.01, .01])
    imagesc(im_room, [0 250]); 
    hold on
    reinforce_window = any(abs(behav_idx - d.room.entrance_start_idx)<=15);
    
    if reinforce_window
        patch(px-showcropx(1), py-showcropy(1), 'g', 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0)
        coloridx = 2;
    else
        patch(px-showcropx(1), py-showcropy(1), 'r', 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', .5)
        coloridx = 1;
    end
    
    imr = imrotate(im, -rad2deg(ad(idx)), 'bilinear', 'crop');
    im_arena = imr(showcropy(1):showcropy(2), showcropx(1):showcropx(2));
    imtop = imr(ztop(1,2):ztop(1,2)+ztop(2,2), ztop(1,1):ztop(1,1)+ztop(2,1));
    imbot = imr(zbot(1,2):zbot(1,2)+zbot(2,2), zbot(1,1):zbot(1,1)+zbot(2,1));

    if length(tail)== length(rx)
        inds = 1:idx-1;
        patch(rx(inds), ry(inds),c(inds),'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .15)
    else
        patch(rx(inds), ry(inds),c,'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .5)
    end
    scatter(rx(behav_idx), ry(behav_idx), 30, 'wo', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor','w', 'MarkerFaceAlpha', .9)
    
    axis image off
    
    [px2_t, px2_r] = cart2pol(px-showcropx(1)-size(im_arena,2)/2, py-showcropy(1)-size(im_arena,1)/2);
    px2_t = mod(px2_t + ad(idx), 2*pi);
    [px2, py2] = pol2cart(px2_t, px2_r);
    px2=px2 + size(im_arena,2)/2;
    py2=py2 + size(im_arena,1)/2;
    subplot_tight(5,2,[4:2:8], [.01, .01])
    imagesc(im_arena, [0 250]); axis image off
    hold on

    if reinforce_window
        patch(px2, py2, 'g', 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0)
        coloridx = 2;
    else
        patch(px2, py2, 'r', 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', .5)
        coloridx = 1;
    end
    if length(tail)== length(rx)
        inds = 1:idx-1;
        patch(ax(inds), ay(inds),c(inds),'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .15)
    else
        patch(ax(inds), ay(inds),c,'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .5)
    end
    scatter(ax(behav_idx), ay(behav_idx), 30, 'wo', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor','w', 'MarkerFaceAlpha', .9)
    inds = behav_idx-tail:behav_idx-1;
    patch(ax(inds), ay(inds),c,'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .5)
    %     inds = 1:idx-1;
    %     patch(x(inds), y(inds),c(inds),'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .15)
    scatter(ax(behav_idx), ay(behav_idx), 30, 'wo', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor','w', 'MarkerFaceAlpha', .9)
    
    if mean(px2_t) > pi/2 && mean(px2_t)<3*pi/2 % left is correct
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.692, .622],  'Y', [.97, .97], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.692, .622],  'Y', [.03, .03], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [.97, .97], 'Color', lightblue, 'LineWidth', 3, 'HeadWidth', 10);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [.03, .03], 'Color', lightblue, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.692, .622],  'Y', [.77, .77], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.692, .622],  'Y', [.23, .23], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [.77, .77], 'Color', lightblue, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [.23, .23], 'Color', lightblue, 'LineWidth', 3, 'HeadWidth', 10);
    else
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.798, .872],  'Y', [.97, .97], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.798, .872],  'Y', [.03, .03], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [.97, .97], 'Color', lightred, 'LineWidth', 3, 'HeadWidth', 10);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [.03, .03], 'Color', lightred, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.798, .872],  'Y', [.77, .77], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.798, .872],  'Y', [.23, .23], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [.77, .77], 'Color', lightred, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [.23, .23], 'Color', lightred, 'LineWidth', 3, 'HeadWidth', 10);
    end
    subplot_tight(5,2,2, [.01, .01])
    imagesc(imtop, [0 250]); axis image off
    subplot_tight(5,2,10, [.01, .01])
    imagesc(imbot, [0 250]); axis image off 
    colormap bone

    drawnow()
end
% hold on
% plot(x,y,'r')

%%















% fname = "C:\Users\gjb326\Desktop\TRACKER DOCS\Tracker Videos\NAPA_24510_IL5_TrackerVideo.avi";
% dname = "C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24510_IL5.mat";
fname = "C:\Users\gjb326\Desktop\TRACKER DOCS\Tracker Videos\NAPA_24511_IL10_TrackerVideo.avi";
dname = "C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\@activeplacepref_NAPA_24511_IL10.mat";
d=load(dname);
v = VideoReader(fname);
% nf = ceil(v.Duration*v.FrameRate);

% nf=2000;
timetoget =60*30;
nf = round(timetoget*v.FrameRate);
ds = 10;
Y = uint8(false(v.Height, v.Width, nf));
vt = zeros(nf,1);

idx = 0;
badcount = 0;
tic
% v.open
while v.hasFrame && v.CurrentTime<timetoget
    im = rgb2gray(v.readFrame);
    if any(im(:)>0)
        idx = idx+1;
        Y(:,:,idx) = im;
        vt(idx) = v.CurrentTime;
    else
        badcount = badcount+1;
    end

end
badcount
Y = Y(:,:,1:idx);

toc
%%
%%%%% making the reinformenet  zone mask
[xt, yt] = meshgrid(1:size(im_arena,2), 1:size(im_arena, 1));
center_x = size(im_arena, 2) / 2;
center_y = size(im_arena, 1) / 2;
xt = xt - center_x;
yt = yt - center_y;
rho = sqrt(xt.^2 + yt.^2);
theta = atan2d(yt, xt); % Angle in degrees

maxextent = 205;
theta = mod(theta+360, 360);
mask = ((rho > 155) & (rho < maxextent)) & ((theta>180-20)&(theta<180+20));
% [py, px] = ind2sub(size(mask),find(mask>0));
% Display the mask
figure
imshow(mask.*im_arena, []);
[px2, py2, q] = ginput(32);

%%
d.params = d.params_sub;
px = [93.9813   91.3854   89.1604   86.9354   85.8229   85.4521   85.4521   86.9354   89.5313   92.1271   94.3521   86.9354   46.8854   43.9188   42.4354   40.9521   39.4688   37.6146   36.5021   35.7604   34.2771   34.4646   34.7604   36.1313   37.6146   38.7271   40.9521   42.4354   43.5479   44.2896   46.1438   46.8854];
py = [293.2146  282.8312  274.6729  263.1771  249.0854  239.8146  232.3979  217.5646  202.7312  196.0562  187.5271  183.4479  169.7271  178.9979  186.0437  191.9771  199.0229  207.5521  217.1937  224.6104  233.1396  243.8937  253.4229  265.0313  274.3021  280.2354  290.6187  295.8104  300.2604  303.2271  308.4188  309.9021];
% px = interp1(1:length(px), px, linspace(1, length(px), 100), 'pchip'); % 100:1:length(px)*100
% py = interp1(1:length(py), py, linspace(1, length(py), 100), 'pchip'); % 100:1:length(px)*100
% px = [119.5383  113.5642  109.8946  111.7294  117.6926  158.0587  155.7651  157.1413  158.9761  162.6457 119.5383  ];
% py = [297.5677  276.8502  247.7317  218.6132  186.5082  203.6806  223.0930  241.7587  256.6912  274.6103 297.5677  ];

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
px = px-10; %py = py-showcropy(1);

%
lr_swap = d.zonestruct.choices.isLeft(1:end-1)==1 & d.zonestruct.choices.isLeft(2:end)==-1;
rl_swap = d.zonestruct.choices.isLeft(2:end)==1   & d.zonestruct.choices.isLeft(1:end-1)==-1;
dir_swap = [0; lr_swap | rl_swap];
dir_swap_idx = find(dir_swap==1);

lose1_win2  = d.zonestruct.choices.isCorrect(dir_swap_idx-1)==-1 & d.zonestruct.choices.isCorrect(dir_swap_idx)==1;
win1_win2   = d.zonestruct.choices.isCorrect(dir_swap_idx-1)==1  & d.zonestruct.choices.isCorrect(dir_swap_idx)==1;
isIDX = d.zonestruct.choices.isIDX(dir_swap_idx);
choice_inds = sort(isIDX(win1_win2));
ts_roi = behavtime(choice_inds);

% f_inds = []; time_win = 10; % seconds pre-post
% for ii = 1:length(ts_roi)
%     t = ts_roi(ii)-time_win*2;
%     vid_idx1 = find(min(abs(vt-t)) == (abs(vt-t)), 1);
%     t = ts_roi(ii)+time_win;
%     vid_idx2 = find(min(abs(vt-t)) == (abs(vt-t)), 1);
%     f_inds = cat(2, f_inds, vid_idx1:vid_idx2);
% end

f_inds = [];
turn_inds = sort([d.zonestruct.choices.L_start_end; d.zonestruct.choices.R_start_end]);
correct_inds = sort(d.zonestruct.choices.Correct_start_end);
error_inds = sort(d.zonestruct.choices.Error_start_end);
for ii = 1:length(choice_inds)
    closest = find( min(abs(turn_inds(:,1)- choice_inds(ii))) == abs(turn_inds(:,1)- choice_inds(ii)));
%     prev = find( min(abs(turn_inds(1:closest-1) - choice_inds(ii))) == abs(turn_inds(1:closest-1) - choice_inds(ii)));
%     next = closest + find( min(abs(turn_inds(closest+1:end) - choice_inds(ii))) == abs(turn_inds(closest+1:end) - choice_inds(ii)));
    thisind = turn_inds(closest,1);
    lasterror = find( error_inds(:,2) < thisind, 1, 'last');
    nexterror = find( error_inds(:,1) > thisind, 1, 'first');
%     if ~any(isempty([closest, prev, next]))
%         f_inds = cat(2, f_inds, prev-time_win:next-time_win);
%     end
    if ~any(isempty([closest, lasterror, nexterror]))
%         length(f_inds)
        f_inds = cat(2, f_inds, error_inds(lasterror,2)-time_win:error_inds(nexterror,1)-time_win);
        length(f_inds);
    end
end

f_inds = unique(f_inds);
%
% close all
figure(10); clf; set(gcf, 'Color', 'k', 'Position', round([400   463   982   521]/2)*2); clf;% imshow(max(Y,[],3))
set(gcf,'renderer','painters')
colormap plasma
tail = round(dt*time_win);
c_tail = magma(tail+1);
c = bone(length(rx));
c = reshape(c, [size(c,1),1,3]);
c(end,:)=NaN;
c(1,:)=NaN;
c_tail = reshape(c_tail, [size(c_tail,1),1,3]);
c_tail(end,:)=NaN;
c_tail(1,:)=NaN;


save_vid = true;
if save_vid
    vv = VideoWriter('C:\Users\gjb326\Desktop\TRACKER DOCS\behavior output\test_winpredict.avi');
    vv.Quality = 80;
    vv.FrameRate = 30;
    vv.open();
end
%

frame2show = f_inds; % [2:4:dt*90, dt*90+1:20:dt*60*9.5, dt*60*9.5+1:4:dt*60*11, dt*60*11+1:20:dt*60*25];
% frame2show = [2:4:dt*60];
lightgreen = [.4 1 .4];
lightblue  = [.4 .4 1];
lightred   = [1 .4 .4];
for i = 1:4:length(frame2show) % 40000+tail:50:nframes
    %%
    clf
    idx = frame2show(i); % xfs(i);
    im = double(squeeze(Y(:,:,idx)));
    vidt = vt(idx);
    behav_idx = find(min(abs(behavtime-vidt)) == (abs(behavtime-vidt)), 1);
    bt = behavtime(behav_idx)/60;

    if ~any(im(:)>0)
        im = lastim;
    else
        lastim = im;
    end
    im_room = im(showcropy(1):showcropy(2), showcropx(1):showcropx(2));

    subplot_tight(1,2,[1], [.01, .01])
%     imagesc(im_room, [0 250]); 
    imshow(im_room, [0 250]); 
    title('\color[rgb]{1 .4 .4}Room Reference Frame', 'interpreter', 'tex')
    hold on
    reinforce_window = any(abs(behav_idx - d.room.entrance_start_idx)<=15);
    
    marks = d.room.entrance_start_idx(d.room.entrance_start_idx<=behav_idx);
    if reinforce_window
%         patch(px-showcropx(1), py-showcropy(1), 'g', 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0)
        patch(px, py, 'g', 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0)
        coloridx = 2;
    else
%         patch(px-showcropx(1), py-showcropy(1), 'r', 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', .5)
        patch(px, py, 'r', 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', .25)
        coloridx = 1;
    end
    
    imr = imrotate(im, -rad2deg(ad(idx)), 'bicubic', 'crop');
    im_arena = imr(showcropy(1):showcropy(2), showcropx(1):showcropx(2));
%     im_arena(im_arena<=10) = 0*255;
    %     if length(c)== length(rx)
    inds = 1:behav_idx-1;
%     patch(rx(inds), ry(inds),c(inds,:,:).*0,'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .05)
    inds_tail = behav_idx-tail:behav_idx;
    inds_tail = inds_tail(inds_tail>0 & inds_tail<nframes);
%     patch(rx(inds_tail), ry(inds_tail), c_tail(1:length(inds_tail),:,:),'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .5)
    patch(rx(inds_tail), ry(inds_tail), c_tail(1:length(inds_tail),:,:),'EdgeColor','interp', 'LineWidth', 2, 'EdgeAlpha', .5)
    scatter(rx(behav_idx), ry(behav_idx), 30, 'wo', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor','w', 'MarkerFaceAlpha', .9)
    scatter(rx(marks), ry(marks), 20, 'wo', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor', lightgreen, 'MarkerFaceAlpha', .5)
    text(size(im_room,2)*.7, size(im_room,2)*.05, sprintf('Time:  %4d', round(behavtime(behav_idx))), 'Color', 'w', 'FontSize', 20)
    
    axis image off
    
%     [px2_t, px2_r] = cart2pol(px-showcropx(1)-size(im_arena,2)/2, py-showcropy(1)-size(im_arena,1)/2);
    [px2_t, px2_r] = cart2pol(px-size(im_arena,2)/2, py-size(im_arena,1)/2);
    px2_t = mod(px2_t + ad(idx), 2*pi);
    [px2, py2] = pol2cart(px2_t, px2_r);
    px2=px2 + size(im_arena,2)/2;
    py2=py2 + size(im_arena,1)/2;
    subplot_tight(1,2,[2], [.01, .01])
%     imagesc(im_arena, [0 250]); axis image off
    imshow(im_arena, [0 250]); axis image off
    title('\color[rgb]{.4 .4 1}Arena Reference Frame')
    hold on

    if reinforce_window
        patch(px2, py2, 'g', 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0)
    else
        patch(px2, py2, 'r', 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', .25)
    end
    
%     patch(ax(inds), ay(inds),c(inds,:,:).*0,'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .05)
    patch(ax(inds_tail), ay(inds_tail), c_tail(1:length(inds_tail),:,:), 'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .5)
    scatter(ax(behav_idx), ay(behav_idx), 30, 'o', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor','w', 'MarkerFaceAlpha', .9) 
    scatter(ax(marks), ay(marks), 30, 'o', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor', lightgreen, 'MarkerFaceAlpha', .5)
%     text(20, 20, sprintf('%2.2f', ad(idx)/pi), 'Color', 'w')
    yy1=.92; yy2 = .08;
    if abs(ad(idx)/pi) < .5 % left is correct
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.691, .619],  'Y', [yy1, yy1], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.691, .619],  'Y', [yy2, yy2], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [yy1, yy1], 'Color', lightgreen, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [yy2, yy2], 'Color', lightgreen, 'LineWidth', 3, 'HeadWidth', 10);
    else
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.799, .871],  'Y', [yy1, yy1], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.799, .871],  'Y', [yy2, yy2], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [yy1, yy1], 'Color', lightgreen, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [yy2, yy2], 'Color', lightgreen, 'LineWidth', 3, 'HeadWidth', 10);
    end
    
%     subplot_tight(3,4,[10 11], [.07 .01])
%     pL = d.zonestruct.metrics.sequence_prob_Left;
%     pC = d.zonestruct.metrics.sequence_prob_Correct;
%     pX = d.zonestruct.metrics.sequence_prob_times / 60;
%     plotidx = 1:ceil(bt/5);
%     plotidx = plotidx(1:min(max(plotidx), length(pX)));
%     hold on
%     plot(pX(plotidx), pL(plotidx), '-', 'Color',  lightblue, 'LineWidth', 3)
%     plot(pX(plotidx), pC(plotidx), '--', 'Color', lightred, 'LineWidth', 3)
%     L1 = scatter(pX(plotidx), pL(plotidx), 150, 's', 'MarkerFaceColor', lightblue, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
%     L2 = scatter(pX(plotidx), pC(plotidx), 100, 'o', 'MarkerFaceColor', lightred, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
%     plot([0 30], [.5 .5], 'k:', 'LineWidth', 2)
%     plot([bt bt], [0 2], 'k-', 'LineWidth', 2)
%     bt = 5*ceil(bt/5);
%     axis([-1, 27, .4, 1.1])
%     ylabel('Probability')
%     xlabel('Epoch (minutes)')
%     legend([L1, L2], 'Correct turn', 'Left turn', "Location", 'northwest')
%     title('\color[rgb]{1 .4 .4}Allocentric \color[rgb]{1 1 1}vs \color[rgb]{.4 .4 1}Egocentric \color[rgb]{1 1 1}performance', 'FontSize', 12)
%     set(gca, 'YColor', 'w', 'XColor', 'w', 'YTick', [.5 .75 1], 'XTick', [0:5:25])
    %%
    drawnow()
    if save_vid==true
%         temp = getframe(gcf);
        temp = getframe(gcf);
        writeVideo(vv, temp)
    end
end
if save_vid==true; 
    vv.close(); 
end
%% Version with two choice ROIs inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10); clf;% imshow(max(Y,[],3))
tail = 120;
c = gray(tail);
% c = gray(length(xfs));
c = reshape(c, [size(c,1),1,3]);
c(end,:)=NaN;
c(1,:)=NaN;
cmap_im = viridis(256);

% c(end,:)=NaN;
% c = reshape(c, [tail,1,3]);
for i = 20000:ds:22000 % 20000+1+tail:ds:nframes
    %%
    idx = xfs(i);
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

    subplot_tight(5,2,[3:2:7], [.01, .01])
    imagesc(im_room, [0 250]); 
    hold on
    reinforce_window = any(abs(behav_idx - d.room.entrance_start_idx)<=15);
    
    if reinforce_window
        patch(px-showcropx(1), py-showcropy(1), 'g', 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0)
        coloridx = 2;
    else
        patch(px-showcropx(1), py-showcropy(1), 'r', 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', .5)
        coloridx = 1;
    end
    
    imr = imrotate(im, -rad2deg(ad(idx)), 'bilinear', 'crop');
    im_arena = imr(showcropy(1):showcropy(2), showcropx(1):showcropx(2));
    imtop = imr(ztop(1,2):ztop(1,2)+ztop(2,2), ztop(1,1):ztop(1,1)+ztop(2,1));
    imbot = imr(zbot(1,2):zbot(1,2)+zbot(2,2), zbot(1,1):zbot(1,1)+zbot(2,1));

    if length(tail)== length(rx)
        inds = 1:idx-1;
        patch(rx(inds), ry(inds),c(inds),'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .15)
    else
        patch(rx(inds), ry(inds),c,'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .5)
    end
    scatter(rx(behav_idx), ry(behav_idx), 30, 'wo', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor','w', 'MarkerFaceAlpha', .9)
    
    axis image off
    
    [px2_t, px2_r] = cart2pol(px-showcropx(1)-size(im_arena,2)/2, py-showcropy(1)-size(im_arena,1)/2);
    px2_t = mod(px2_t + ad(idx), 2*pi);
    [px2, py2] = pol2cart(px2_t, px2_r);
    px2=px2 + size(im_arena,2)/2;
    py2=py2 + size(im_arena,1)/2;
    subplot_tight(5,2,[4:2:8], [.01, .01])
    imagesc(im_arena, [0 250]); axis image off
    hold on

    if reinforce_window
        patch(px2, py2, 'g', 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0)
        coloridx = 2;
    else
        patch(px2, py2, 'r', 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', .5)
        coloridx = 1;
    end
    if length(tail)== length(rx)
        inds = 1:idx-1;
        patch(ax(inds), ay(inds),c(inds),'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .15)
    else
        patch(ax(inds), ay(inds),c,'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .5)
    end
    scatter(ax(behav_idx), ay(behav_idx), 30, 'wo', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor','w', 'MarkerFaceAlpha', .9)
    inds = behav_idx-tail:behav_idx-1;
    patch(ax(inds), ay(inds),c,'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .5)
    %     inds = 1:idx-1;
    %     patch(x(inds), y(inds),c(inds),'EdgeColor','interp', 'FaceAlpha', 0, 'LineWidth', 2, 'EdgeAlpha', .15)
    scatter(ax(behav_idx), ay(behav_idx), 30, 'wo', 'MarkerEdgeColor','none', ...
        'MarkerFaceColor','w', 'MarkerFaceAlpha', .9)
    
    if mean(px2_t) > pi/2 && mean(px2_t)<3*pi/2 % left is correct
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.692, .622],  'Y', [.97, .97], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.692, .622],  'Y', [.03, .03], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [.97, .97], 'Color', lightblue, 'LineWidth', 3, 'HeadWidth', 10);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [.03, .03], 'Color', lightblue, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.692, .622],  'Y', [.77, .77], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.692, .622],  'Y', [.23, .23], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [.77, .77], 'Color', lightblue, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.69, .62], 'Y', [.23, .23], 'Color', lightblue, 'LineWidth', 3, 'HeadWidth', 10);
    else
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.798, .872],  'Y', [.97, .97], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.798, .872],  'Y', [.03, .03], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [.97, .97], 'Color', lightred, 'LineWidth', 3, 'HeadWidth', 10);
%         ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [.03, .03], 'Color', lightred, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.798, .872],  'Y', [.77, .77], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.798, .872],  'Y', [.23, .23], 'Color', [0 0 0], 'LineWidth', 4, 'HeadWidth', 11);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [.77, .77], 'Color', lightred, 'LineWidth', 3, 'HeadWidth', 10);
        ar = annotation("arrow", 'Units', 'normalized', 'X', [.8, .87],  'Y', [.23, .23], 'Color', lightred, 'LineWidth', 3, 'HeadWidth', 10);
    end
    subplot_tight(5,2,2, [.01, .01])
    imagesc(imtop, [0 250]); axis image off
    subplot_tight(5,2,10, [.01, .01])
    imagesc(imbot, [0 250]); axis image off 
    colormap bone

    drawnow()
end




