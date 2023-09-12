clear 
topdir = 'C:\Users\gjb326\Documents\HISTOLOGY\DarrylConfocal_PKMZ\';
animals = {'21627' '21628' '22911A' '22911B' '24500'  '24501'  '24502'  '24503'};
mag = '40x'; % or '20x'
region = 'HIPP'; % or '20x'
prefix = 'Garrett Images ACC-HIPP';



dapi = cell(length(animals), 1);
pkc = cell(length(animals), 1);
figure(77); clf
for aa = 1:length(animals)
    f1 = sprintf('%s\\%s %s %s\\MERGED\\%s_%s %s %s MERGED_ch00.tif', topdir, animals{aa}, mag, region, prefix, animals{aa}, mag, region);
    f2 = sprintf('%s\\%s %s %s\\MERGED\\%s_%s %s %s MERGED_ch02.tif', topdir, animals{aa}, mag, region, prefix, animals{aa}, mag, region);
    
    im1 = imread(f1);
    subplot(2,length(animals), aa)
    imagesc(im1)
    dapi{aa} = im1;
    im2 = imread(f2);
    subplot(2,length(animals), aa+length(animals))
    imagesc(im2)
    pkc{aa} = im2;
    drawnow
    
end
mag_adjust = str2double(mag(1:2))/20;
numpoints = 800*mag_adjust;
width = 30*mag_adjust;
ksize = 3*mag_adjust;
bg_val = 20;%.15; % to ignore holes/bubbles if included on line
f1 = figure(1); clf

%%
% for aa = 1:length(animals)
%     im1 = single(dapi{aa});
%     im1 = conv2(im1, ones(ksize,ksize)./(ksize^2), 'same');
%     im2 = single(pkc{aa});
%     im2 = conv2(im2, ones(ksize,ksize)./(ksize^2), 'same');
%     
%     temp = nanmean(cat(3,im1, im2),3);
%     temp(temp<bg_val) = 400;
%     subplot(2,4,aa); cla; hold on
%     imagesc(temp);
% 
%     drawnow
% end
%%
ca1x = NaN(length(animals), 2); ca1y = NaN(length(animals), 2);
ca3x = NaN(length(animals), 2); ca3y = NaN(length(animals), 2);

ca1vals_647 = NaN(length(animals), numpoints, width*2 + 1);
ca1vals_dapi = NaN(length(animals), numpoints, width*2 + 1);
ca3vals_647 = NaN(length(animals), numpoints, width*2 + 1);
ca3vals_dapi = NaN(length(animals), numpoints, width*2 + 1);
saveVars = {'ca1*' 'ca3*' 'ctx*' 'numpoints', 'width', 'mag_adjust', 'ksize'};
% f2 = figure(2);
%%  
for aa = aa:length(animals)
%%
im1 = single(dapi{aa});
im1 = conv2(im1, ones(ksize,ksize)./(ksize^2), 'same');
im2 = single(pkc{aa});
im2 = conv2(im2, ones(ksize,ksize)./(ksize^2), 'same');

figure(f1); clf
subplot(3,1,3)
title('CA1!')
[ca1vals_647(aa,:,:), ca1vals_dapi(aa,:,:), ca1x(aa, :), ca1y(aa, :)] = get_im_vals(im1, im2, numpoints, width, f1);
title('CA3!')
[ca3vals_647(aa,:,:), ca3vals_dapi(aa,:,:), ca3x(aa, :), ca3y(aa, :)] = get_im_vals(im1, im2, numpoints, width);
% title('ALVEUS!')
% [alveus_647(aa,:,:), alveus_dapi(aa,:,:), alvx(aa, :), alvy(aa, :)] = get_im_vals(im1, im2, numpoints, width);
% get_im_vals(im1, im2, numpoints, width, f1);
% get_im_vals(im1, im2, numpoints, width);

end
%%
alvx = NaN(length(animals), 2); alvy = NaN(length(animals), 2);
alveus_647 = NaN(length(animals), numpoints, width*2 + 1);
alveus_dapi = NaN(length(animals), numpoints, width*2 + 1);

ctxx = NaN(length(animals), 2); ctxy = NaN(length(animals), 2);
ctx_647 = NaN(length(animals), numpoints, width*2 + 1);
ctx_dapi = NaN(length(animals), numpoints, width*2 + 1);

thalx = NaN(length(animals), 2); thaly = NaN(length(animals), 2);
thal_647 = NaN(length(animals), numpoints, width*2 + 1);
thal_dapi = NaN(length(animals), numpoints, width*2 + 1);

for aa = 1:length(animals)
%%
im1 = single(dapi{aa});
im1 = conv2(im1, ones(ksize,ksize)./(ksize^2), 'same');
im2 = single(pkc{aa});
im2 = conv2(im2, ones(ksize,ksize)./(ksize^2), 'same');

figure(f1); clf
subplot(3,1,3)
title('CTX!')
[ctx_647(aa,:,:), ctx_dapi(aa,:,:), ctxx(aa, :), ctxy(aa, :)] = get_im_vals(im1, im2, numpoints, width, f1);
title('THAL!')
[thal_647(aa,:,:), thal_dapi(aa,:,:), thalx(aa, :), thaly(aa, :)] = get_im_vals(im1, im2, numpoints, width);
title('ALVEUS!')
[alveus_647(aa,:,:), alveus_dapi(aa,:,:), alvx(aa, :), alvy(aa, :)] = get_im_vals(im1, im2, numpoints, width);
% get_im_vals(im1, im2, numpoints, width, f1);
% get_im_vals(im1, im2, numpoints, width);

end
%%
% save('C:\Users\gjb326\Documents\HISTOLOGY\DarrylConfocal_PKMZ\dw_ver2.mat', saveVars{:});
if false
figure(2); clf
im11 = im1; im22 = im2;
ca1vals_647 = NaN(numpoints, width*2 + 1);
ca1vals_dapi = NaN(numpoints, width*2 + 1);

im11 = normalize_matrix(im11);
im22 = normalize_matrix(im22);
im11 = imrotate(im11, -6);
im22 = imrotate(im22, -6);
ksize = 5;
% hsmall = fspecial('disk', 50);
% bg2 = imfilter(im22, hsmall);
% bg1 = imfilter(im11, hsmall);
% bg2 = imgaussfilt(im22, ksize);
% bg1 = imgaussfilt(im11, ksize);
% bg1 = conv2(im11, ones(ksize,ksize)./(ksize^2), 'same');
bg2 = conv2(im22, ones(ksize,ksize)./(ksize^2), 'same');
% im11 = double(im11)-bg1;
% im22 = double(im22)-bg2;

im = zeros(size(im11,1), size(im11,2), 3);
im(:,:,3) = im11./2 + im22;
im(:,:,2) = im11./2;
im(:,:,1) = im22./1.3;

image(im*2); axis image
[px, py] = ginput(2);
px = round(px); py = round(py);
% axis([min(px) max(px) min(py) max(py)])
imnew = im(min(py):max(py), min(px):max(px),:);
bg3 = 1.4-bg2(min(py):max(py), min(px):max(px),:);
bg3 = cat(3, bg3, 1.3.*bg3, bg3.*.65);
% axis([min(px) max(px) min(py) max(py)])
imnew2 = 2.7*imnew.*bg3;
image([imnew; imnew2])
% imwrite(imnew2, 'C:\Users\gjb326\Documents\HISTOLOGY\test.jpeg')
% imwrite(imnew2, 'C:\Users\gjb326\Documents\HISTOLOGY\test.tif')
end
%%
% base_647 = nanmean(nanmean(thal_647,3),2)*ones(1, numpoints);
% base_dapi = nanmean(nanmean(thal_dapi,3), 2)*ones(1, numpoints);
% % ctx  = nanmean(nanmean(ctx_647,3),2);
% % alv  = nanmean(nanmean(alveus_647,3),2);
% % thal = nanmean(nanmean(thal_647,3),2);
% % base_647 = alv*ones(1, numpoints);


% base_647 = nanmean(nanmean(alveus_647,3),2)*ones(1, numpoints);
% base_647 = nanmean(nanmean(alveus_647,3),2)*ones(1, numpoints);
% base_dapi = nanmean(nanmean(alveus_dapi,3), 2)*ones(1, numpoints);
% base_ca1_dapi = nanmean(nanmean(ca1vals_647(:,1:100,:),3), 2)*ones(1, numpoints);
% base_ca3_dapi = nanmean(nanmean(ca3vals_647(:,1:100,:),3), 2)*ones(1, numpoints);

% ca1_647 = squeeze(nanmedian(ca1vals_647,3))./base_647;
% ca1_dapi = squeeze(nanmedian(ca1vals_dapi,3))./base_ca1_dapi;
% ca3_647 = squeeze(nanmedian(ca3vals_647,3))./base_647;
% ca3_dapi = squeeze(nanmedian(ca3vals_dapi,3))./base_ca3_dapi;
ca1_647 = squeeze(nanmedian(ca1vals_647,3));
ca1_dapi = squeeze(nanmedian(ca1vals_dapi,3));
ca3_647 = squeeze(nanmedian(ca3vals_647,3));
ca3_dapi = squeeze(nanmedian(ca3vals_dapi,3));

% ca3_647 = normalize_rows(ca3_647);
% ca3_dapi = normalize_rows(ca3_dapi);

figure(11); clf
subplot(1,2,1)
plot(ca1_647(1:4,:)', 'k'); hold on; plot(ca1_647(5:end,:)', 'r');
subplot(1,2,2)
plot(ca3_647(1:4,:)', 'k'); hold on; plot(ca3_647(5:end,:)', 'r');

figure(12); clf
subplot(1,2,1)
plot(ca1_dapi(1:4,:)', 'k'); hold on; plot(ca1_dapi(5:end,:)', 'b');
subplot(1,2,2)
plot(ca3_dapi(1:4,:)', 'k'); hold on; plot(ca3_dapi(5:end,:)', 'b');

figure(9); clf; 
subplot(2,2,1); hold on; 
title('CA1')
% plot(nanmean(ca1_dapi(1:4,:),1), 'k:');
% plot(nanmean(ca1_dapi(5:8,:),1), 'r:')
% plot(nanmean(ca1_647(1:4,:),1), 'k');
% plot(nanmean(ca1_647(5:8,:),1), 'r')
plot(ca1_647(1:4,:)', 'k');
plot(ca1_647(5:8,:)', 'r');


subplot(2,2,2); hold on; 
title('CA3')
% plot(nanmean(ca3_dapi(1:4,:),1), 'k:');
% plot(nanmean(ca3_dapi(5:8,:),1), 'r:')
% plot(nanmean(ca3_647(1:4,:),1), 'k');
% plot(nanmean(ca3_647(5:8,:),1), 'r')
plot(ca3_647(1:4,:)', 'k');
plot(ca3_647(5:8,:)', 'r');

subplot(2,2,3); hold on; 
title('CA1')
nsubs = 4; 
trained = 5:8;
untrained = 1:4;
shadedErrorBar(1:numpoints, nanmedian(ca1_647(untrained,:),1), nanstd(ca1_647(untrained,:),[],1)./sqrt(length(untrained)), 'lineprops', 'k');
shadedErrorBar(1:numpoints, nanmedian(ca1_647(trained,:),1), nanstd(ca1_647(trained,:),[],1)./sqrt(length(trained)), 'lineprops', 'r');

subplot(2,2,4); hold on; 
title('CA3')
shadedErrorBar(1:numpoints, nanmedian(ca3_647(untrained,:),1), nanstd(ca3_647(untrained,:),[],1)./sqrt(length(untrained)), 'lineprops', 'k');
shadedErrorBar(1:numpoints, nanmedian(ca3_647(trained,:),1), nanstd(ca3_647(trained,:),[],1)./sqrt(length(trained)), 'lineprops', 'r');




function [pvals_647, pvals_dapi, x, y] = get_im_vals(im1, im2, numpoints, width, h)
im1_og = im1; im2_og = im2;
%%
pvals_647 = NaN(numpoints, width*2 + 1);
pvals_dapi = NaN(numpoints, width*2 + 1);
im1 = im1_og; im2 = im2_og;

im1 = normalize_matrix(im1);
im2 = normalize_matrix(im2);


im = zeros(size(im1,1), size(im1,2), 3);
im(:,:,3) = im1./2 + im2;
im(:,:,2) = im1./2;
im(:,:,1) = im2./1.3;
%%

%
if nargin>4
figure(h);
subplot(3,1,1:2)
image(im*2); axis image
[px, py] = ginput(2);
px = round(px); py = round(py);
axis([min(px) max(px) min(py) max(py)])
else
subplot(3,1,1:2)
end
% image(im(min(py):max(py), min(px):max(px),:)*2)

[x, y] = ginput(2);

dx = x(1)-x(2);
dy = y(1)-y(2);
ang = (atan(dx/dy));
xdiff = sin(ang+pi/2)*width;
ydiff = cos(ang+pi/2)*width;

xx = linspace(x(1), x(2), numpoints);
yy = linspace(y(1), y(2), numpoints);

%%
x1 = round(xx-xdiff);
x2 = round(xx+xdiff);
y1 = round(yy-ydiff);
y2 = round(yy+ydiff);
hold on
plot(x, y, 'w', 'linewidth', 1)
plot(x1, y1,'r:', 'linewidth', 2)
plot(x2, y2,'r:', 'linewidth', 2)
for i = 1:numpoints
    a1 = round(x1(i));
    a2 = round(x2(i));
    a = ( round(linspace(a1, a2, width*2)) );
    b1 = round(y1(i));
    b2 = round(y2(i));
    b = ( round(linspace(b1, b2, width*2)) );
    valid = ( a>0 & a<size(im2,2) ) & ( b>0 & b<size(im2,1) );
    a = a(valid);
    b = b(valid);
    c = sqrt(diff(a).^2 + diff(b).^2);
    a = a(c>0);
    b = b(c>0);
%     plot([a1, a2], [b1, b2], 'r', 'LineWidth', .5)
%     scatter(a,b,'ro')
%     im2(a,b) = 1000;
    for j = 1:length(a)
        pvals_647( i, j) = im2_og(b(j),a(j));
        pvals_dapi( i, j) = im1_og(b(j),a(j));
    end
end
%%
% bg_val = 20;%.15; % to ignore holes/bubbles if included on line
% 
% temp = nanmean(cat(3,im1_og, im2_og),3);
% temp(temp<bg_val) = NaN;
% [h,w,~] = size(im);
% im2 = reshape(im(:,:,1), [h*w,1]);
% im2(isnan(temp(:))) = 1;
% im2 = reshape(im2, [h, w]);
% im(:,:,1) = im2;
% 
% subplot(3,1,1:2); cla; hold on
% image(im);
% pvals_dapi(pvals_dapi<bg_val) = NaN;
% pvals_647(pvals_647<bg_val)   = NaN;

subplot(3,1,3); cla; hold on
plot(nanmean(pvals_dapi,2), 'b')
plot(nanmean(pvals_647,2), 'r')

end




















