function hisogram_rgb(im, bins)
if nargin < 2
    bins = [0:.01:2];
end
figure; 
hold on; 
rgb='rgb'; 
for i = 1:3
    histogram(im(:,:,i), bins, 'FaceColor', rgb(i), 'EdgeAlpha', 0); 
end