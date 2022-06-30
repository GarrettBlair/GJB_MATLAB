anums = 8;
% anums = [6 7 8 9 30 31 32 34 35 36 37];
c = 0;
ddir = 'C:\Users\gjb326\Desktop\Sample Data\MiniData Hipp Blair\';

tt = [.5 1 2 4 8 16 32];
nspks = normalize_matrix(S(:,1:2048));
[nsegs, nsamples] = size(nspks);
sumspks = NaN(length(tt), nsegs, nsamples);
pvcorr_t = NaN(length(tt), nsamples, nsamples);
pvcorr_p = NaN(length(tt), nsamples, nsamples);
pvcorr_z = NaN(length(tt), nsamples, nsamples);

for ttsub = 1:length(tt)
    %%
time_res = round((1/11) * 11 * tt(ttsub)); % 1/11 is the scope sampling rate after half frame downsampling

subplot_c = ceil(sqrt(length(anums)));
subplot_r = floor(sqrt(length(anums)));
% for cc = 1:length(anums)
    
fname = sprintf('%sHipp%d\\caiman_cnmfe_out.mat', ddir, anums(cc));
load(fname)

% nspks = normalize_matrix(S);
smooth_nspks = conv2(1, ones(1,5), nspks, 'same');
smooth_nspks = normalize_matrix(smooth_nspks);
%
nsub = floor(nsamples/time_res)-1;
tvec = 1:time_res:nsamples;
%%%% Calculate similarity by fixed bins
% sumspks = NaN(nsegs, nsub);
% pvcorr_t = NaN(nsub, nsub);
% pvcorr_p = NaN(nsub, nsub);
% for i = 1:nsub
%     sumspks(:, i) = sum(nspks(:, tvec(i):tvec(i+1)), 2);
%     sumspks(:, i) = sum(smooth_nspks(:, tvec(i):tvec(i+1)), 2);
% end
% for i = 1:nsub
% for j = i:nsub
%    [pvcorr_t(i, j),  pvcorr_p(i, j)] = corr(sumspks(:, i), sumspks(:, j), 'type', 'Kendall');
%    pvcorr_t(j, i) = pvcorr_t(i, j);
%    pvcorr_p(j, i) = pvcorr_p(i, j);
% end
% end
%%%% Calculate similarity by sliding bins
for i = 1:nsamples-time_res
%     sumspks(:, i) = sum(nspks(:, i:i+time_res), 2);
    sumspks(ttsub, :, i) = sum(smooth_nspks(:, i:i+time_res), 2);
end
for i = 1:nsamples
for j = i:nsamples
   [pvcorr_t(ttsub, i, j),  pvcorr_p(ttsub, i, j)] = corr(sumspks(ttsub, :, i)', sumspks(ttsub, :, j)', 'type', 'Kendall');
   pvcorr_t(ttsub, j, i) = pvcorr_t(ttsub, i, j);
   pvcorr_p(ttsub, j, i) = pvcorr_p(ttsub, i, j);
end
end

fakepos = ones(size(smooth_nspks,2), 1);
p_neighbors_vec=[0.075/15 0.075];
rand('seed',0)
ReducedDataArray=DimentionalityReduction_Ver1(1*(smooth_nspks>0), p_neighbors_vec);
v2=ReducedDataArray{3}';

options.dims = 3; %1:3 %range of dimensions to calculate typ 2:20 for dim calculations, 3 for plot
options.landmarks = 1:size(smooth_nspks,2); %recommended value
nNeighbours = 5;
D = L2_distance(smooth_nspks, smooth_nspks, 1);
[Y, R, E] = IsomapII(D, 'k', nNeighbours, options); 
v1 = Y.coords{3,1};


% pcorr_z = pvcorr_t;
temp = squeeze(pvcorr_p(ttsub,:,:));
temp(find(eye(size(temp))==1)) = NaN;
pvcorr_z(ttsub, :, :) = temp;
end
%%
figure(1);clf
% subplot_tight(3,2,[1 2]);
% imagesc(nspks); axis tight
% subplot_tight(3,2,[3 4]);
% imagesc(sumspks); axis tight
% c = c+1;
for i = 1:length(tt)
subplot_tight(2, length(tt), i, [.03 .01]);
pvc = squeeze(pvcorr_t(i,:,:));
imagesc(pvc, [-.1 .5]); colormap viridis; colorbar
% yyaxis('right')
% plot(10*sum(smooth_nspks,1) + nsamples); 
axis image off
set(gca, 'YDir', 'normal')
title(sprintf('Kendall Tau, t=%2.1f', tt(i)))
% set(gca, 'XTick', [0:60*11:nsub], 'XTickLabel', [0:60*11:nsub]*11)
subplot_tight(2, length(tt), i+length(tt), [.03 .01]);
pvc = squeeze(pvcorr_z(i,:,:));
imagesc(-1.*log10(pvc), [-2 10]); colormap viridis; axis square; colorbar
set(gca, 'YDir', 'normal')
title('Kendall Tau Sig')
axis image off
drawnow
end
% yyyyyyy
%%
fname_out = sprintf('%sHipp%d\\popvec_corr_multi_integration_times.mat', ddir, anums(cc));
% save(fname_out, 'pvcorr_t', 'pvcorr_p', 'time_res','sumspks', 'nspks', 'smooth_nspks');
% load(fname_out, 'pvcorr_t', 'pvcorr_p', 'time_res','sumspks', 'nspks', 'smooth_nspks');
% end
%%

im_file1 = 'C:\Users\gjb326\Desktop\Sample Data\MiniData Hipp Blair\Hipp8\sub_MC.tiff';
im_file2 = 'C:\Users\gjb326\Desktop\Sample Data\MiniData Hipp Blair\Hipp8\sub_MC_dff.tiff';
prop_f = sum(nspks,1);
prop_f(1) = prop_f(2);
prop_f = prop_f./max(prop_f);
t = (1/11)*(1:nsamples)/60;
    figure(2); clf
    set(gcf, 'Color', [.2 .2 .2], 'Position', [59 110 1039 846])
%     v = VideoWriter('C:\Users\gjb326\Desktop\Sample Data\MiniData Hipp Blair\Hipp8\pcov_vis.avi', 'Uncompressed AVI');
    v = VideoWriter('C:\Users\gjb326\Desktop\Sample Data\MiniData Hipp Blair\Hipp8\pcov_vis2.avi');
    v.Quality = 80;
%     v.VideoCompressionMethod = 100;
    v.open();
for i = 1:2:nsamples
    %%
    f1 = imread(im_file1, i);
    f2 = imread(im_file2, i);
    
%     subplot_tight(1,3,1);
    figure(2); clf
    subplot_tight(2,2,1, [.03 .001]);
    imshow(f1);
    
    subplot_tight(2,2,3, [.03 .001]);
    imshow(f2);
    
    subplot_tight(5, 2, 2, [0.03 .01]);
    cla;
    hold on
    plot(t, prop_f, 'k'); 
    plot([t(i) t(i)], [-.1 1], 'r-'); 
    set(gca, 'YColor', 'w', 'XColor', 'w')
    xlabel('Time (min)')
    axis tight
    
    subplot_tight(6, 2, [6:2:12], [0.03 .01])
    cla; 
    hold on
    imagesc(pvcorr_t, [-.1 .5]); colormap viridis;
    set(gca, 'YDir', 'normal')
    plot([i i], [0 nsamples], 'r-'); 
    plot([0 nsamples], [i i], 'r-'); 
    axis tight square off;
    
    drawnow
    temp = getframe(gcf);
    c = temp.cdata;
    v.writeVideo(c)
end
v.close();











