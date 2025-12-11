function [lastx, lasty, lastfnum] = fernando_helper_click_mini(im, cropX, cropY, fnum)

figure(101); clf;
set(gcf,'Color', [.7 1 1])
im2 = uint8(im(cropY(1):cropY(2), cropX(1):cropX(2),:));
image(im2); axis image
title('Click miniscope LED')
[lastx, lasty] = ginput(1);
lastx = round(lastx);
lasty = round(lasty);
lastfnum = fnum-3;
% ledVal = mean(mean(mean(im2(lasty-2:lasty+2, lastx-2:lastx+2, :))));

end