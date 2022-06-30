clear 
% fname = 'C:\Users\gjb326\Desktop\Sample Data\Hipp15\sub_MC.tiff';
fname = 'C:\Users\gjb326\Desktop\Sample Data\Hipp18\sub_MC.tiff';
% finfo = imfinfo(fname);

% tic
% TIFFStack appears slightly slower and more complicated variable output
% m1 = TIFFStack(fname);
% toc


tic
m2 = imread_big(fname);
toc

%%
maxf = max(m2, [], 3);
meanf = uint8(mean(m2, 3));
minf = uint8(mean(m2, 3));

figure(1) % ; imagesc(maxf - meanf)

[h, w, numF] = size(m2);

for i = 1:2:numF
    f = m2(:,:,i);
    fs = uint8((f - minf) + 128);
    imshow([f fs]);
    drawnow
end