toyrat_dir  = 'F:\Fernando\ToyRat cellreg\';
movgrid_dir = 'F:\Fernando\MovGrid cellreg\';
animals = ['A', 'B', 'C', 'D', 'E', 'F', 'G'];

retrieval_sess = 2;
names = {'toyrat_spatial_footprints_01', 'toyrat_spatial_footprints_02', 'fearretrieval_spatial_footprints_02'};
for aLoop = 2 % 1:length(animals)-1
    
    fnames = cell(3,1);
    fnames{1} = [toyrat_dir animals(aLoop) '\spatial_footprints_01.mat'];
    fnames{2} = [toyrat_dir animals(aLoop) '\spatial_footprints_02.mat'];
    fnames{3} = [movgrid_dir animals(aLoop) '\aligned_data_struct.mat'];
    temp = load(fnames{1});
    s1 = temp.A;
    [nsegs1, h1, w1] = size(s1);
    max_pixels = 0;
    [s1, bads] = norm_A_conts(s1, .4, max_pixels);
    %
    temp = load(fnames{2});
    s2 = temp.A;
    [nsegs2, h2, w2] = size(s2);
    [s2, bads] = norm_A_conts(s2, .4, max_pixels);

    cr = load(fnames{3});
    
    s3 = cr.aligned_data_struct.spatial_footprints{retrieval_sess};
    [s3, bads] = norm_A_conts(s3, .65, max_pixels);
    nsegs2show=[];
    [contours_shifted, projections_shifted, down_right_corrections, scaling_corrections] = manual_contour_matching_footprints_only({s1,s2,s3}, [], [], nsegs2show, []);
        fprintf('\n__________________\n\t%s', animals(aLoop))
    %%
    for i = 1:length(contours_shifted)
        aligned_fn = sprintf('%s%s\\aligned_spatial_footprints_0%d.mat', toyrat_dir,  animals(aLoop), i);
        sessname = names{i};
        filename = fnames{i};
        
        fprintf('\n%s ', sessname)
        fprintf('\n [y,x][%d, %d],   %1.2f', down_right_corrections(i,1), down_right_corrections(i,2), scaling_corrections(i))
        fprintf('\n__________________')
        
        
        contours = contours_shifted{i};
        vars2save = {'contours', 'sessname', 'filename'};
        save(aligned_fn, vars2save{:})
    end

end

%%


% A [-2 9]
% B [ 8 -2] maybe [x3 y0]
% C [ 4 3]
% D [ 1 -8]
% E [ 4 3]
% F [-2 -3]
% G [-8 10]
%%
% temp = load('C:\Users\gjb326\Downloads\spatial_footprints_01.mat');
% s2 = temp.A;




function [s1, bads] = norm_A_conts(s1, normthresh, npx)
[nsegs1, h1, w1] = size(s1);
for i = 1:nsegs1
    c = squeeze(s1(i, :, :));
    if any(c(:)>0)
        c = c-min(c(:));
        c = c./max(c(:));
        s1(i,:,:) = c;
    end
end
bads = true(nsegs1,1);
for i = 1:nsegs1
    c = squeeze(s1(i, :, :));
    if any(c(:)>0)
        bads(i) = false;
        if npx>0
        [~,ord] = sort(c(:), 'descend');
        c(ord(npx+1:end)) = 0;
        end
        c = c-min(c(:));
        c = c./max(c(:));
        if npx==0
            c(c<normthresh) = 0;
        end
        s1(i,:,:) = c;
    end
end
if sum(bads)>0
    warning('bads = %d', sum(bads))
end

end

function [yOffset, xOffset] = xcorr_gjb(img1, img2)
% Perform normalized cross-correlation
correlation = normxcorr2(img1, img2);

% Find peak in correlation matrix
[maxCorr, maxIndex] = max(abs(correlation(:)));
[yPeak, xPeak] = ind2sub(size(correlation), maxIndex);

% Calculate offset between images
yOffset = yPeak - size(img1, 1);
xOffset = xPeak - size(img1, 2);

% Display results
fprintf('Best alignment offset: x = %d, y = %d\n', xOffset, yOffset);

% Optional: visualize the matching location
figure;
im = cat(3, img1, img1*0, img2);
image(im);
hold on;
rectangle('Position', [xOffset, yOffset, size(img1,2), size(img1,1)], ...
          'EdgeColor','r', 'LineWidth', 2);
title('Best match location in larger image');

end