% d1 = 'F:\GarrettBlair\APA\HPCACC24504\2024_01_24\16_07_41_TR14\HPC_miniscope1\msCam_MC.tiff';
% d2 = 'F:\GarrettBlair\APA\HPCACC24504\2024_01_24\16_07_41_TR14\ACC_miniscope2\msCam_MC.tiff';
% % f1 = 'F:\GarrettBlair\APA\HPCACC24504\2024_01_24\16_07_41_TR14\experiment\2024_01_24_H16_07_41_TR14_@placecells_HPC_miniscope1.mat';
% % f2 = 'F:\GarrettBlair\APA\HPCACC24504\2024_01_24\16_07_41_TR14\experiment\2024_01_24_H16_07_41_TR14_@placecells_ACC_miniscope2.mat';
% f1 = 'D:\GarrettBlair\APA\HPCACC24504\processed_files\2024_01_24_H16_07_41_TR14_@placecells_HPC_miniscope1.mat';
% f2 = 'D:\GarrettBlair\APA\HPCACC24504\processed_files\2024_01_24_H16_07_41_TR14_@placecells_ACC_miniscope2.mat';
% ms1 = load(f1);
% ms2 = load(f2);


% d1 = 'D:\GarrettBlair\APA\HPCACC24502_TR18_msCam_MC.tiff';
% f1 = 'D:\GarrettBlair\APA\HPCACC24502\simple_files\2023_06_26_H11_39_11_TR18_@placecells_HPC_miniscope1.mat';
d1 = 'F:\GarrettBlair\APA\HPCACC24504\2024_02_19\14_51_00_TDF23\HPC_miniscope1\msCam_MC.tiff';
f1 = 'F:\GarrettBlair\APA\HPCACC24504\simple_files\2024_02_19_H14_51_00_TDF23_@placecells_HPC_miniscope1.mat';
d2 = 'D:\GarrettBlair\APA\Acc20832_NEW5_msCam_MC.tiff';
f2 = 'D:\GarrettBlair\APA\Acc20832\simple_files\2023_01_26_H16_10_07_NEW5_@placecells.mat';
ms1 = load(f1);
ms2 = load(f2);

% c1 = load("D:\GarrettBlair\APA\HPCACC24502\matching_contours\2023_06_26_H11_39_11_TR18_@contours_HPC_miniscope1.mat");
c1 = load("F:\GarrettBlair\APA\HPCACC24504\matching_contours\2024_02_19_H14_51_00_TDF23_@contours_HPC_miniscope1.mat");
c2 = load("D:\GarrettBlair\APA\Acc20832\matching_contours\2023_01_26_H16_10_07_NEW5_contours.mat");
c2.contours = c2.contours(:, 1:2:end, 1:2:end);


% Y1 = zeros(ms1.ms.height, ms1.ms.width, length(ms1.ms.frameNum));
% Y2 = zeros(ms1.ms.height, ms1.ms.width, length(ms1.ms.frameNum));
Y1 = imread_big(d1);
Y2 = imread_big(d2);
Y2 = Y2(1:2:end, 1:2:end, :);
%
raw_traces = ms1.craw; % ms1.ms.neuron.C + ms1.ms.neuron.YrA;
[dff_filt, ~] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, false);
spks1 = normalize_rows(dff_filt);
% spks1(isnan(spks1)) = 0;
% 
raw_traces = ms2.craw; %  = ms2.ms.neuron.C + ms2.ms.neuron.YrA;
[dff_filt, ~] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, false);
spks2 = normalize_rows(dff_filt);
% spks2(isnan(spks2)) = 0;
% t1 = ms1.ms.timestamps./1000;
% t2 = ms2.ms.timestamps./1000;

% spks1 = ms1.spks;
% spks2 = ms2.spks;
t1 = ms1.time_ms./1000;
t2 = ms2.time_ms./1000;

t1 = t1-t1(1);
t2 = t2-t2(1);

integral_time = .5;
[spks1, gt1] = average_spks_time(spks1, integral_time, t1', false, 'sum');
[spks2, gt2] = average_spks_time(spks2, integral_time, t2', false, 'sum');

ns = min(size(spks2,2), size(spks1,2));
Y1sub = zeros(size(Y1,1), size(Y1,2), ns);
Y2sub = zeros(size(Y2,1), size(Y2,2), ns);
for i =1:ns
    Y1sub(:,:,i) = squeeze(mean(Y1(:,:,i==gt1), 3));
    Y2sub(:,:,i) = squeeze(mean(Y2(:,:,i==gt2), 3));
end
% Y1n=Y1sub;
% Y2n=Y2sub;
Y1sub = Y1sub - quantile(Y1sub(:), .05);
Y1sub(Y1sub<0) = 0;
Y1sub = Y1sub ./ quantile(Y1sub(:), .9999);
Y1sub(Y1sub>1) = 1;

Y2sub = Y2sub - quantile(Y2sub(:), .05);
Y2sub(Y2sub<0) = 0;
Y2sub = Y2sub ./ quantile(Y2sub(:), .9999);
Y2sub(Y2sub>1) = 1;

bg1 = nanmean(Y1sub, 3);
bg2 = nanmean(Y2sub, 3);

clearvars Y2 Y2n Y1 Y1n
%%
figure(1); clf;
set(gcf, 'Color', 'w')
colormap magma
accscale = .2;
hpcscale = .1;
% hpccrop = [30, 80]; hpccrop = cat(2, hpccrop, hpccrop+150);
% acccrop = [25, 25]; acccrop = cat(2, acccrop, acccrop+150);
% hpccrop = [1, 1]; hpccrop = cat(2, hpccrop, [size(Y1sub,1)-1, size(Y1sub,2)-1]);
% acccrop = [1, 1]; acccrop = cat(2, acccrop, [size(Y2sub,1)-1, size(Y2sub,2)-1]);
hpccrop = [65, 85]; hpccrop = cat(2, hpccrop, hpccrop+100);
acccrop = [28, 25]; acccrop = cat(2, acccrop, acccrop+100);
makevid=false;
if makevid==true
v = VideoWriter('C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\_visualization\output\imaging_dff.avi');
v.Quality = 100;
v.FrameRate = 5*(1/integral_time);
v.open();
end
for i = 1:1:900%ns
    %%
    figure(1); clf;
    subplot_tight(2,2,1,[.05 .05])
    imagesc(Y1sub(:,:,i), [.4 1]); axis image off
    axis([hpccrop(2) hpccrop(4) hpccrop(1) hpccrop(3)])
    t_string = sprintf('Time = %i sec', round(i*integral_time));
    text( hpccrop(4)-20, hpccrop(3)+10, t_string)
    
    subplot_tight(2,2,2,[.05 .05])
%     imagesc(Y1sub(:,:,i)-double(ms1.ms.neuron.meanFrame), [-hpcscale hpcscale])
    imagesc(Y1sub(:,:,i)-bg1, [-hpcscale hpcscale]); axis image off
    axis([hpccrop(2) hpccrop(4) hpccrop(1) hpccrop(3)])
    
    subplot_tight(2,2,3,[.05 .05])
    imagesc(Y2sub(:,:,i), [.2 1]); axis image off
    axis([acccrop(2) acccrop(4) acccrop(1) acccrop(3)])
    
    subplot_tight(2,2,4,[.05 .05])
%     imagesc(Y2sub(:,:,i)-double(ms2.ms.neuron.meanFrame), [-accscale accscale])
    imagesc(Y2sub(:,:,i)-bg2, [-accscale accscale]); axis image off
    axis([acccrop(2) acccrop(4) acccrop(1) acccrop(3)])
    drawnow
    if makevid==true
        temp = getframe(gcf);
        v.writeVideo(temp.cdata)
    end
    pause(.01)

end
if makevid==true
v.close()
end


%%
figure(2); clf;
set(gcf, 'Color', 'w')
colormap magma
[ns1, h1, w1] = size(c1.contours);
[ns2, h2, w2] = size(c2.contours);
conts1 = zeros(h1*w1, ns1);
conts2 = zeros(h2*w2, ns2);
for i = 1:ns1
    conts1(:,i) = reshape(squeeze(c1.contours(i,:,:)), [h1*w1, 1]);
end
for i = 1:ns2
    conts2(:,i) = reshape(squeeze(c2.contours(i,:,:)), [h2*w2, 1]);
end

% cont1 = normalize_cols(ms1.ms.neuron.fullA);
% cont2 = normalize_cols(ms2.ms.neuron.fullA);
cont1 = normalize_cols(conts1);
cont2 = normalize_cols(conts2);
% cont1(cont1<.6)=0;
% cont2(cont2<.6)=0;

ksize=10;
slowkern = [zeros(1,ksize), linspace(1,0,ksize)]; slowkern = slowkern/sum(slowkern);
slowspks1 = conv2(1, slowkern, spks1, 'same');
slowspks2 = conv2(1, slowkern, spks2, 'same');

% slowspks1 = spks1>0; % normalize_rows(ms1.craw);
% slowspks2 = spks2>0; % normalize_rows(ms2.craw);
for i =100:1:200%ns
    %%
    a1 = cont1*slowspks1(:,i);
    a2 = cont2*slowspks2(:,i);
%     a1 = reshape(a1, [ms1.ms.neuron.dims(1), ms1.ms.neuron.dims(2)]);
%     a2 = reshape(a2, [ms2.ms.neuron.dims(1), ms2.ms.neuron.dims(2)]);
    a1 = reshape(a1, [h1, w1]);
    a2 = reshape(a2, [h2, w2]);

    figure(2); clf;
    subplot_tight(2,3,1,[.05 .05])
    imagesc(Y1sub(:,:,i));
    axis([hpccrop(2) hpccrop(4) hpccrop(1) hpccrop(3)])
    set(gca, 'YTick',[], 'XTick', [])
    
    t_string = sprintf('Time = %i sec', round(i*integral_time));
    text( hpccrop(4)-20, hpccrop(3)+10, t_string)
    
    subplot_tight(2,3,2,[.05 .05])
%     df = Y1sub(:,:,i)-mean(Y1sub(:,:,i-10:i-1), 3);
    df = Y1sub(:,:,i)-bg1;
    imagesc(df, [-hpcscale hpcscale])
    axis([hpccrop(2) hpccrop(4) hpccrop(1) hpccrop(3)])
    set(gca, 'YTick',[], 'XTick', [])
    
    df(df<-hpcscale)=-hpcscale;
    df(df> hpcscale)= hpcscale;
    df = normalize_matrix(df)*.5 + .25;
    subplot_tight(2,3,3,[.05 .05])
    image(cat(3, df, df+a1, df+a1))
    axis([hpccrop(2) hpccrop(4) hpccrop(1) hpccrop(3)])
    set(gca, 'YTick',[], 'XTick', [])
    
    
    subplot_tight(2,3,4,[.05 .05])
    imagesc(Y2sub(:,:,i))
    axis([acccrop(2) acccrop(4) acccrop(1) acccrop(3)])
    set(gca, 'YTick',[], 'XTick', [])
    
    subplot_tight(2,3,5,[.05 .05])
%     df = Y2sub(:,:,i)-mean(Y2sub(:,:,i-10:i-1), 3);
    df = Y2sub(:,:,i)-bg2;
    imagesc(df, [-accscale accscale])
    axis([acccrop(2) acccrop(4) acccrop(1) acccrop(3)])
    set(gca, 'YTick',[], 'XTick', [])
    
    df(df<-accscale)=-accscale;
    df(df> accscale)= accscale;
    df = normalize_matrix(df)*.5 + .25;
    subplot_tight(2,3,6,[.05 .05])
    image(cat(3, df+a2, df, df+a2))
    axis([acccrop(2) acccrop(4) acccrop(1) acccrop(3)])
    set(gca, 'YTick',[], 'XTick', [])
    
    drawnow
    pause(.01)

end
%%
figure(3); clf;
set(gcf, 'Color', 'w')
colormap magma
nsegs = 100;
backwin=60;

dfax1 = subplot_tight(2,4,1,[.05 .05]); cla
dfax2 = subplot_tight(2,4,5,[.05 .05]); cla
contax1 = subplot_tight(2,4,2,[.05 .05]); cla
contax2 = subplot_tight(2,4,6,[.05 .05]); cla

ax1 = subplot_tight(2,4,3:4,[.05 .05]);
s = normalize_rows(spks1(1:nsegs,:));
stacked_traces(s, 10, {'Color', [0, .1, .3], 'LineWidth', .1}, [.9 1 1])

ax2 = subplot_tight(2,4,7:8,[.05 .05]);
s = normalize_rows(spks2(1:nsegs,:));
stacked_traces(s, 10, {'Color', [.3, 0, .1], 'LineWidth', .1}, [1 .9 1])
%     xlim([100-backwin, 100+backwin*5])
makevid=false;
if makevid==true
v = VideoWriter('C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\_visualization\output\dff_contours_traces.avi');
v.Quality = 100;
v.FrameRate = 5*(1/integral_time);
v.open();
end
cropspd = 10; % higher slower zoom out
for i =60:2:600 % backwindow*2:1:ns-backwindow*5
    %%
    b1 = cont1*(.1*ones(ns1,1));
    b2 = cont2*(.1*ones(ns2,1));
    a1 = cont1*(slowspks1(:,i));
    a2 = cont2*(slowspks2(:,i));
%     a1 = reshape(a1, [ms1.ms.neuron.dims(1), ms1.ms.neuron.dims(2)]);
%     a2 = reshape(a2, [ms2.ms.neuron.dims(1), ms2.ms.neuron.dims(2)]);
    b1 = reshape(b1, [h1, w1]);
    b2 = reshape(b2, [h2, w2]);
    a1 = reshape(a1, [h1, w1]);
    a2 = reshape(a2, [h2, w2]);

    inds = i-backwin:i+backwin*5;
    inds = inds( inds>0 & inds<=ns);
%     figure(3); clf;
    
    df = Y1sub(:,:,i)-bg1;
    subplot(dfax1);
    cla
    imagesc(df, [-hpcscale hpcscale])
    axis image off
    crop = [max(hpccrop(2)-i/cropspd, 1), min(hpccrop(4)+i/cropspd, min(w1,h1)),...
            max(hpccrop(1)-i/cropspd, 1), min(hpccrop(3)+i/cropspd,  min(w1,h1))];
%     axis([hpccrop(2) hpccrop(4) hpccrop(1) hpccrop(3)])
    axis(crop)
    subplot(contax1); cla
    df(df<-hpcscale)=-hpcscale;
    df(df> hpcscale)= hpcscale;
    df = normalize_matrix(df)*.5 + .25;
    df = df/2;
    image(cat(3, df+b1, df+a1+b1, df+a1+b1))
    axis image
    axis(crop)
    
    set(gca, 'YTick',[], 'XTick', [])
%     t_string = sprintf('Time = %i sec', round(i*integral_time));
%     text( hpccrop(4)-20, 10, t_string)
    
%     subplot_tight(2,3,2:3,[.05 .05])
%     s = cat(3, spks1(:,inds), spks1(:,inds)+slowspks1(:,inds), spks1(:,inds)+slowspks1(:,inds));
%     s = cat(3, spks1(1:nsegs,inds), spks1(1:nsegs,inds), spks1(1:nsegs,inds));
%     image(s)

%     s = spks1(1:nsegs,inds);
%     stacked_traces(s, 2, {'Color', [0, .1, .3], 'LineWidth', .1}, [.9 1 1])
%     hold on
%     plot([backwin+1, backwin+1], [-5 nsegs+10], 'r-')
%     axis([i-backwin, i+backwin*5, 0, nsegs*1.1])

    
    df = Y2sub(:,:,i)-bg2;    
    subplot(dfax2)
    imagesc(df, [-accscale accscale])
    axis image off
    crop = [max(acccrop(2)-i/cropspd, 1), min(acccrop(4)+i/cropspd, min(w2,h2)),...
            max(acccrop(1)-i/cropspd, 1), min(acccrop(3)+i/cropspd,  min(w2,h2))];
    axis(crop)

    subplot(contax2); cla
    df(df<-accscale)=-accscale;
    df(df> accscale)= accscale;
    df = normalize_matrix(df)*.5 + .25;
    df = df/2;
    image(cat(3, df+a2+b2, df+b2, df+a2+b2))
    axis image
    axis(crop)
    set(gca, 'YTick',[], 'XTick', [], 'Color', 'none')

    
%     subplot_tight(2,3,5:6,[.05 .05])
%     s = cat(3, spks2(1:nsegs,inds)+slowspks2(1:nsegs,inds), spks2(1:nsegs,inds), spks2(1:nsegs,inds)+slowspks2(1:nsegs,inds));
%     image(s)

%     s = spks2(1:nsegs,inds);
%     stacked_traces(s, 2, {'Color', [.3, 0, .1], 'LineWidth', .1}, [1 .9 1])
%     hold on
%     plot([backwin+1, backwin+1], [-5 nsegs+10], 'r-')
    set(ax1, 'YTick',[0:20:nsegs], 'XTick', [i], 'XTickLabel', round(i*integral_time), 'Color', 'none',...
        'XLim',[i-backwin, i+backwin*5], 'YLim',[-5, nsegs*1.2])
    set(ax2, 'YTick',[0:20:nsegs], 'XTick', [i], 'XTickLabel', round(i*integral_time), 'Color', 'none',...
        'XLim',[i-backwin, i+backwin*5], 'YLim',[-5, nsegs*1.2])

    drawnow
    if makevid==true
        temp = getframe(gcf);
        v.writeVideo(temp.cdata)
    end
    pause(.01)

end
if makevid==true
v.close()
end












%%
Ycw = imread_big("D:\Sample Data\cortical_window_1P\megha_cw_mouse\mouse_cw2_rec1\msCam_MC.tiff");
[h, w, f] = size(Ycw);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gt = ceil(linspace(1, f, f)/4);

Ysub = zeros(h,w,max(gt));

for i =1:max(gt)
    Ysub(:,:,i) = squeeze(mean(Ycw(:,:,i==gt), 3));
end
% Y1n=Y1sub;
% Y2n=Y2sub;
Ysub = Ysub - quantile(Ysub(:), .05);
Ysub(Ysub<0) = 0;
Ysub = Ysub ./ quantile(Ysub(:), .9999);
Ysub(Ysub>1) = 1;

bg = nanmean(Ysub, 3);
%%
figure(4); clf;
set(gcf, 'Color', 'w')
colormap magma
hpcscale = .1;

makevid=true;
if makevid==true
v = VideoWriter('C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\_visualization\output\cw_dff_slow.avi');
v.Quality = 100;
v.FrameRate = 20;
v.open();
end
k = ones(30,30); k = k./(sum(k(:)));
k2 = ones(500,500); k2 = k2./(sum(k2(:)));
fbg = conv2(bg, k2, 'same');
x = [513   450   352   265   204   172   196   307   430   471   443   394   269   149    72];
y = [259   269   297   354   417   519   599   610   586   515   434   393   317   221   118];
s = [300   300   300   300   300   350   350   350   400    400   450   450   700   850   1000];
a=[4,3];a=a/max(a);
n = 2*max(gt)/3;

ii = [linspace(1, n, length(x))];
xi = interp1(ii, x, linspace(1,n,n), 'pchip');
yi = interp1(ii, y, linspace(1,n,n), 'pchip');
si = interp1(ii, s, linspace(1,n,n), 'pchip');
sii = conv(si, ones(30,1)/30, 'valid');
si(15:end-15) = sii;
    figure(4); clf;
ax = subplot_tight(1,1,1, [0,0]);
for i = 1:1:max(gt)
    %%
    im1 = Ysub(:,:,i);
    df = im1-bg;
    
    dffbg = df - conv2(df, k, 'same');
%     df = imrotate(df, -8);
%     df = normalize_matrix(df);
%     df = df - .5;
%     imagesc(im1, [.1 1]); axis image off
    subplot(ax); cla
    imagesc(im1 + dffbg - fbg, [-.2 .5]); axis image off
    if i < n
        crop = [xi(i), xi(i)+si(i)*a(1), yi(i), yi(i)+si(i)*a(2)];
    else
        crop = [xi(end), xi(end)+si(end)*a(1), yi(end), yi(end)+si(end)*a(2)];
    end
    axis(crop)
    drawnow
    if makevid==true
        temp = getframe(gcf);
        v.writeVideo(temp.cdata)
    else
        pause(.01)        
    end

end
if makevid==true
v.close()
end








