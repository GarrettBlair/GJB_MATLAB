%%
topdir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\EWL_testing\Hipp16942\2022_07_01\12_54_51\MiniLFOV';
tiffname = sprintf('%s\\msCam.tiff', topdir);
bgname = sprintf('%s\\bg_filt41.tiff', topdir);
notefname = sprintf('%s', 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\EWL_testing\Hipp16942\2022_07_01\12_54_51\notes.csv');
msTSFile = sprintf('%s/timeStamps.csv', topdir);
tiffout = sprintf('%s\\test.tiff', topdir);
ds = 2;

Y = imread_big(tiffname);
bg = imread(bgname);
bg = imresize(bg, .5);
% bg = mean(Y,3);
% bg = uint8(bg);
n = readtable(notefname);
ms = readtable(msTSFile);
% f1 = ;

ms_ts = ms.TimeStamp_ms_(1:ds:end);
ewl_ts = n.Var1;
ewl = n.Var4;


ewl_ms = interp1(ewl_ts, ewl, ms_ts, 'last', 'extrap');
ewl_ms(ewl_ms>127) = 127;
ewl_ms(ewl_ms<-127) = -127;
ewl_ms = floor(ewl_ms);
ewl_ms(ewl_ms==-1) = 0;
%% Making video of changing EWL
v = Fast_Tiff_Write(tiffout);
%
figure(9); clf
set(gcf, 'Position', [276   422   716   374], 'Color', 'w')
for i = 1:2:size(Y,3)
    %%
    im = squeeze(Y(:,:,i));
    im = imresize(im, .5);
    imdff = 4*(im-(bg));
    figure(9); clf
    subplot_tight(1,2,1, [.15, .075]); cla
    hold on;
    scatter(ms_ts(i)./1000, ewl_ms(i), 300, '.r');
    plot(ms_ts./1000, ewl_ms, 'k');
    xlabel('Seconds')
    ylabel('EWL value')
    axis square tight
    set(gca, 'XTick', [15:15:75], 'YTick', -127:127:127, 'YLim', [-150 150])
    
    subplot_tight(1,2,2, [.15, .075]); cla
    
    imshow(im, [0 255])
    text(size(im,2)*.65, size(im,1)*.05 - 40, sprintf('EWL: %3.0f', ewl_ms(i)), 'Color', 'r', 'FontSize', 14)
    text(size(im,2)*.05, size(im,1)*1.075, sprintf('LED:  3'), 'Color', 'r', 'FontSize', 10)
    text(size(im,2)*.05, size(im,1)*1.025, sprintf('Gain: Medium'), 'Color', 'r', 'FontSize', 10)
    drawnow
    temp = getframe(gcf);
    im_out = temp.cdata;
    im_out = permute(im_out, [2 1 3]);
    v.WriteIMG(im_out)
end
%
v.close()

%% Second, longer recording
% uses a 500x500 square, topleft edge at y=176 x=446
nframes = 100;
ewlNEG0 = 1240;%1;
ewlNEG127 = 770;
ewlNEG71  = 1650;
ewlNEG33 = 2200;
ewlPOS127 = 407;
ewlNEG0_2 = 3600;%1;
ewlvals = [0 -33 -71 -127 127 0];
ewlNEG0   = ewlNEG0:ewlNEG0+nframes-1;
ewlNEG127 = ewlNEG127:ewlNEG127+nframes-1;
ewlNEG71  = ewlNEG71:ewlNEG71+nframes-1;
ewlNEG33  = ewlNEG33:ewlNEG33+nframes-1;
ewlPOS127   = ewlPOS127:ewlPOS127+nframes-1;
ewlNEG0_2   = ewlNEG0_2:ewlNEG0_2+nframes-1;
roi_start_y = 176;
roi_start_x = 446;

topdir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\EWL_testing\Hipp16942\2022_07_01\12_59_52\MiniLFOV';
tiffname = sprintf('%s\\msCam.tiff', topdir);

notefname = sprintf('%s', 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\EWL_testing\Hipp16942\2022_07_01\12_59_52\notes.csv');
msTSFile = sprintf('%s/timeStamps.csv', topdir);
ds = 2;
n = readtable(notefname);
ms = readtable(msTSFile);
ms_ts = ms.TimeStamp_ms_(1:ds:end);
ewl_ts = n.Var1;
ewl = n.Var4;
ewl_ms = interp1(ewl_ts, ewl, ms_ts, 'nearest', 'extrap');
ewl_ms(ewl_ms>127) = 127;
ewl_ms(ewl_ms<-127) = -127;
figure;
plot(ewl_ms, 'k');


Y = imread_big(tiffname);
Y1 = Y(roi_start_y:roi_start_y+499, roi_start_x:roi_start_x+499, ewlNEG0);
Y2 = Y(roi_start_y:roi_start_y+499, roi_start_x:roi_start_x+499, ewlNEG33);
Y3 = Y(roi_start_y:roi_start_y+499, roi_start_x:roi_start_x+499, ewlNEG71);
Y4 = Y(roi_start_y:roi_start_y+499, roi_start_x:roi_start_x+499, ewlNEG127);
Y5 = Y(roi_start_y:roi_start_y+499, roi_start_x:roi_start_x+499, ewlPOS127);
Y6 = Y(roi_start_y:roi_start_y+499, roi_start_x:roi_start_x+499, ewlNEG0_2);

% Y1bg = imresize(min(Y1,[],3), .5);
% Y2bg = imresize(min(Y2,[],3), .5);
% Y3bg = imresize(min(Y3,[],3), .5);
% Y4bg = imresize(min(Y4,[],3), .5);
% Y1bg = uint8(imresize(mean(Y1,3), .5));
% Y2bg = uint8(imresize(mean(Y2,3), .5));
% Y3bg = uint8(imresize(mean(Y3,3), .5));
% Y4bg = uint8(imresize(mean(Y4,3), .5));
%% Making video of multiple planes
tiffout = sprintf('%s\\test2.tiff', topdir);
v = Fast_Tiff_Write(tiffout);
%
figure(10); clf
q = 3;
iz = 250;
set(gcf, 'Color', 'w', 'Position', [300, 400, iz*3 + q*2, iz*2 + q])
hborder = ones(iz,q)*255;
wborder = ones(q,2*q+ 3*iz)*255;
xs = 20:iz+q:20+2*(iz+q);
ys = 20:iz+q:20+1*(iz+q);
for i = 1:nframes
    %%
    for j = 1:6
    eval(sprintf('im%d = squeeze(Y%d(:,:,i));', j, j))
    eval(sprintf('im%d = imresize(im%d, .5);', j, j))
    end
%     im = 15.*[im1-Y1bg, hborder, im2-Y2bg; wborder; im3-Y3bg, hborder, im4-Y4bg]+30;
%     im = [im1, hborder, im2; wborder; im3, hborder, im4];
    im = [im1, hborder, im2, hborder, im3; wborder; im4, hborder, im5, hborder, im6];
%     imshow(im, [50 200]);
    imshow(im, [0 255]);
    text(xs(1), ys(1), sprintf('EWL: %3.0f', ewlvals(1)), 'Color', 'w', 'FontSize', 14)
    text(xs(2), ys(1), sprintf('EWL: %3.0f', ewlvals(2)), 'Color', 'w', 'FontSize', 14)
    text(xs(3), ys(1), sprintf('EWL: %3.0f', ewlvals(3)), 'Color', 'w', 'FontSize', 14)
    text(xs(1), ys(2), sprintf('EWL: %3.0f', ewlvals(4)), 'Color', 'w', 'FontSize', 14)
    text(xs(2), ys(2), sprintf('EWL: %3.0f', ewlvals(5)), 'Color', 'w', 'FontSize', 14)
    text(xs(3), ys(2), sprintf('EWL: %3.0f', ewlvals(6)), 'Color', 'w', 'FontSize', 14)
    drawnow
    temp = getframe(gcf);
    im_out = temp.cdata;
    im_out = permute(im_out, [2 1 3]);
    v.WriteIMG(im_out)
end
%
v.close()
