function [cropX, cropY, lastfnum] = fernando_helper_click_behavROI(im)

    figure(101); clf;
    set(gcf,'Color', [.7 .7 1])
    image(uint8(im)); axis image
    title('Crop the behavior valid area, click 2 opposing corners')
    [cropX, cropY] = ginput(2);
    cropX = round([min(cropX) max(cropX)]);
    cropY = round([min(cropY) max(cropY)]);
    lastfnum=-1;

end
