function show_mscam_dff(fname, frames2show)
Y = imread_big(fname);
[h,w,nf] = size(Y);
minY = min(Y,[],3);
if nargin<2
    frames2show = 1:4:nf;
end
figure;
imshow(minY);
[px,py,resp]=ginput(2);

px = round(sort(px)); py = round(sort(py));

Ysub = Y(py(1):py(2), px(1):px(2), frames2show);
Ydff = Y(py(1):py(2), px(1):px(2), frames2show) - minY(py(1):py(2), px(1):px(2));
for i=1:length(frames2show)
    figure(1); clf;
    subplot_tight(2,1,1, [0,0]);
    imagesc(squeeze(Ysub(:,:,i)), [0 190]); axis image off;
    subplot_tight(2,1,2, [0,0]);
    imagesc(squeeze(Ydff(:,:,i)), [-5 30]); axis image off; drawnow;
end


end