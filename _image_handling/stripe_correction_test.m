% Y = load_tiffstack("C:\Users\gjb326\Desktop\test_stripes\test.tiff");
Yog = load_tiffstack("C:\Users\gjb326\Desktop\test_stripes\test2.tiff");
[h,w,nf] = size(Yog);
Ysub = Yog(200:2:500, 200:2:400,:);
Yc = reshape(Ysub, [size(Ysub,1)*size(Ysub,2), nf]);
Y = permute(Yog, [2 1 3]);
% [h,w,nf] = size(Y);
Yt = reshape(Y, [h*w, nf]);
c = corr(double(Yc));
bads = find(c(:,1)<.95);
% for i = 1:size(Y,3)
%     
% end
%%

for ii = 1:length(bads)
    %%
    fnum = bads(ii)
refFrame = squeeze(Y(:,:,858));
frame = squeeze(Y(:,:,fnum));
refLine1 = squeeze(Yt(:,800));
refLine2 = squeeze(Yt(:,858));
frame_line = squeeze(Yt(:,fnum));
frame_line_copy = frame_line;

x1 = 169;
y1 = 418;
x2 = 448;
y2 = 431;

idx1 = w*(y1-1) + x1;
idx2 = w*(y2-1) + x2;
buff_size = idx2-idx1;
idx_odds = idx1 + (buff_size)*(2*linspace(-100,100, 201));
idx_odds = idx_odds( (idx_odds > 0) & (idx_odds < (h*w)) );

idx_evs = idx1 + (buff_size)*(2*linspace(-100,100, 201)) - (buff_size);
idx_evs = idx_evs( (idx_evs > 0) & (idx_evs < (h*w)) );

% for i = 1:length(idx_odds)
%     a = idx_odds(i);
%     b = a+buff_size;
%     if b>(h*w)
%         bdiff = b-(h*w);
%         frame_line(a:b) = 0;
%         frame_line(1:bdiff) = 0;
%     else
%         frame_line(a:b) = 0;
%     end    
% end
for i = 1:length(idx_evs)
    % values from prior stripe
    a2diff = 0;
    a1 = idx_evs(i);
    a2 = a1+buff_size;
    % index to stripe
    b2diff = 0;
    if i+1>length(idx_evs)
    b1 = idx_evs(1);
    b2 = b1+buff_size;
    else
    b1 = idx_evs(i+1);
    b2 = b1+buff_size;
    end
    
    if a2>(h*w)
        a2diff = a2-(h*w);
        a2 = (h*w);
%         vals = [frame_line_copy(a1:a2); frame_line_copy(1:a2diff)];
        inds = [a1:a2, 1:a2diff];
    else
%         vals = [frame_line_copy(a1:a2)];
        inds = [a1:a2];
    end
    if b2>(h*w)
        b2diff = b2-(h*w);
        b2 = (h*w);
%         inds = [b1:b2, 1:b2diff];
        vals = [frame_line_copy(b1:b2); frame_line_copy(1:b2diff)];
    else
%         inds = [b1:b2];
        vals = [frame_line_copy(b1:b2)];
    end
    frame_line(inds) = vals;
end
% Yt2(idx1:idx2) = 0;
yyy2 = reshape(frame_line', [h,w]);
% figure; imagesc(yyy2);
% drawnow
% pause(1)
Yog(:,:,fnum) = yyy2';
end
f = Fast_Tiff_Write("C:\Users\gjb326\Desktop\test_stripes\test3.tiff");
for i = 1:1000
f.WriteIMG(Yog(:,:,i))
end
f.close()
% figure; imshowpair(yyy2, YY1)


%%
figure; hold on; 
% plot(conv(refLine1, ones(100,1)./100, 'same')); 
bg = uint8(median(Yt,2));
plot(abs(diff(frame_line_copy-bg)))














