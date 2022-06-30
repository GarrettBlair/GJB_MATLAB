figure(1); clf

frames = 1:2:size(C,2);

YY = uint8(fullA*C(:,frames));
YR = uint8(zeros(size(YY)));

ddir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\2022_06_18\16_44_33\MiniLFOV';
fname = sprintf('%s/msCam_MC.tiff', ddir);
for i = frames
%     t = imread_big(fname, i);
    YR(:,i) = reshape(imread(fname, i), [dims(1)*dims(2), 1]);
    
end
%%
for i = frames
    temp = YY(:,i);
    im = reshape(temp, dims);
    temp = YR(:,i);
    imr = reshape(temp, dims);
    subplot_tight(1,2,1, [0 0])
    imagesc(im, [-4 10]);
    axis image off
    subplot_tight(1,2,2, [0 0])
    imagesc(imr, [30 250]);
    axis image off
    colormap gray
    drawnow
    
end