% uiopen('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\2023_01_11\15_56_35_TR5\MiniLFOV\Plot Values.csv',1)
%%
mcx = PlotValues1.Mean(1:end);
fname = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Acc19947\2023_01_11\15_56_35_TR5\MiniLFOV\msCam_MC.tiff';
[ugp, ugn] = id_bad_frames(mcx);
% if any(intersect(ugn, ugp))
%     x = intersect(ugn, ugp);
%     ugn = ugn(ugn~=x);
%     ugp = ugp(ugp~=x);
% end
%%
show_bad_frames(ugp, ugn, fname)

%%
function [ugp, ugn] = id_bad_frames(mcx)
d = abs(diff(mcx));
h_thresh = median(d) + 10*std(d);
l_thresh = median(d) + 2*std(d);

bad_ind = find(d>h_thresh);
gi_prev = NaN(length(bad_ind),1);
gi_next = NaN(length(bad_ind),1);
for i = 1:length(bad_ind)
    bi = bad_ind(i);
    gi_prev(i) = find(d(1:bi)<= l_thresh, 1, 'last');
    gi_next(i) = find(d(bi:end)<= l_thresh, 1, 'first') + bi - 1;
    
end
    
    
figure(6); clf; 
hold on
plot(d)
plot(bad_ind, d(bad_ind), 'ro')
ugp = unique(gi_prev);
ugn = unique(gi_next);
plot(ugp, d(ugp), 'go')
plot(ugp, d(ugp), 'gx')
plot(ugn, d(ugn), 'bo')
end
%%
function show_bad_frames(ugp, ugn, fname)
for i = 1:length(ugp)
    span = ugp(i):ugn(i);
    figure(7); clf
    for j = 1:length(span)
        f = imread(fname, span(j));
        subplot_tight(1, length(span), j, [0 0])
        imshow(f)
        axis image tight
        title(span(j))
    end
    drawnow;
    [~, ~, resp] = ginput(1);
    if resp == 113
        return
    end
        
end
end