


% px = ms.room.svm_decoding.pred_x;
% py = ms.room.svm_decoding.pred_y;
% rx = ms.room.svm_decoding.x;
% ry = ms.room.svm_decoding.y;
px = ms.room.pfield_decode.xdecoded;
py = ms.room.pfield_decode.ydecoded;
rx = ms.room.pfield_decode.xbin;
ry = ms.room.pfield_decode.ybin;
t = ms.room.svm_decoding.t;

figure; hold on;
histogram(ms.room.svm_decoding.pred_err)
histogram(ms.room.pfield_decode.decode_dist)
histogram(ms.room.svm_decoding.rand_err_median)

bads = ms.room.svm_decoding.pred_err > nanmedian(ms.room.svm_decoding.pred_err);

figure; histogram2(px(bads), py(bads), ms.params.pos_bins, ms.params.pos_bins)
figure; histogram2(px(~bads), py(~bads), ms.params.pos_bins, ms.params.pos_bins)
figure; histogram2(rx(bads), ry(bads), ms.params.pos_bins, ms.params.pos_bins)
figure; histogram2(rx(~bads), ry(~bads), ms.params.pos_bins, ms.params.pos_bins)

% mt = ms.timestamps./1000;
% framei = interp1(mt, ms.frameNum, t, 'nearest')
%%
win = 10;
clr = jet(win + 1);
for i = 1+win:2:length(t)-win
    fs = i-win:i;
    figure(1); clf; hold on;
    plot(rx, ry, 'Color', [.5 .5 .5])
    plot(rx(fs), ry(fs), 'k-')
    plot(px(fs), py(fs), 'r-')
    scatter(px(fs), py(fs), 60, clr, '.')
    plot(rx(i), ry(i), 'kx')
%     axis([-30 30 -30 30])
    axis([-1 21 -1 21])
    drawnow
    
end