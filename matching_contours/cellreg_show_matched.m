clear 

% ddir = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_HPC\cellreg_subIL\";
% crfile = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_HPC\cellreg_subIL\cellRegistered_20251014_165122.mat";
ddir = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_ACC\cellreg_subIL\";
crfile = "D:\GarrettBlair\APA\HPCACC34990\matching_contours\manual_alignment_ACC\cellreg_subIL\cellRegistered_20251014_165409.mat";

% ddir = "D:\GarrettBlair\APA\HPCACC24500\matching_contours\manual_alignment_ACC\cellreg\";
% crfile = "D:\GarrettBlair\APA\HPCACC24500\matching_contours\manual_alignment_ACC\cellreg\cellRegistered_20231116_151132.mat";

% ddir = "D:\GarrettBlair\APA\HPCACC24500\matching_contours\manual_alignment_ACC\cellreg\";
% crfile = "D:\GarrettBlair\APA\HPCACC24500\matching_contours\manual_alignment_HPC\cellreg\cellRegistered_20231116_152054.mat";

temp = load(crfile);
cmap = temp.cell_registered_struct.cell_to_index_map;
ns = size(cmap, 2);
%%

figure; 
subplot(1,3,1)
[o, l] = pop_overlap(cmap);
imagesc(cmap)
subplot(1,3,2)
imagesc(o, [0, 1])
colorbar
subplot(1,3,3)
shadedErrorBar(1:ns, nanmean(l,2), nanstd(l,[],2))
ylim([0 1])
colormap magma
%%
[~,ny,nx] = size(temp.cell_registered_struct.spatial_footprints_corrected{1});
close all
contour_thresh = .65;
for i = 2% 1:ns
%     shared = cmap(:,j)>0 & cmap(:,i)>0;
    c_im = zeros(ny,nx,3);
    idx=0;
    for j = [i-1:i+1]
        idx=idx+1;
        if j>0 && j<=ns
            if j == i
                if j==1
                    shared = (cmap(:,i+1)>0 & cmap(:,i)>0);
                elseif j==ns
                    shared = (cmap(:,i-1)>0 & cmap(:,i)>0);
                else
                    shared = (cmap(:,i+1)>0 & cmap(:,i)>0) | (cmap(:,i-1)>0 & cmap(:,i)>0);
                end
            else
                shared = cmap(:,j)>0 & cmap(:,i)>0;
            end
            segs = cmap(shared,j);
            cs = temp.cell_registered_struct.spatial_footprints_corrected{j}(segs,:,:);
            for k = 1:size(cs,1); 
                ccs = squeeze(cs(k,:,:));
                ccs = ccs-min(ccs(:));
                ccs = ccs./max(ccs(:));
                ccs(ccs<=contour_thresh) = 0;
                cs(k,:,:) = ccs;
            end
                
%             cs(cs<=.4)=0;
            c_im(:,:, mod(j,3)+1) = sum(cs,1).*1;
        end
    end    
    figure(i); 
    image(c_im)
    title(i)
end
%%
[~,ny,nx] = size(temp.cell_registered_struct.spatial_footprints_corrected{1});
sessns = 1:3;
contour_thresh = .5;
ns=length(sessns);
shared = sum(cmap(:,sessns)>0, 2) == ns; %  & cmap(:,i)>0

sharedvals = [.4 .4 .4];
othervals  = [.4 .4 .4];
figure(1001); clf
set(gcf, 'Position', [200, 200, 1000, 450])
for i = 1:ns
    %     shared = cmap(:,j)>0 & cmap(:,i)>0;
    c_im = zeros(ny,nx,3);
    c_im = zeros(ny,nx,3);
    segs1 = cmap(shared,i);
    segs2 = cmap(shared==false & cmap(:,i)>0, i);
    im1 = seg_im(temp.cell_registered_struct.spatial_footprints_corrected{i}(segs1,:,:), contour_thresh);

    im2 = seg_im(temp.cell_registered_struct.spatial_footprints_corrected{i}(segs2,:,:), contour_thresh);
    
    c_im(:,:, i) = im1*.95;
    for ii = 1:3; c_im(:,:, ii) = squeeze(c_im(:,:, ii)) + squeeze(im2*.4); end
    for ii = 1:3
        if i~=ii
            c_im(:,:, ii) = squeeze(c_im(:,:, ii)) + squeeze(im1*.2); 
        end
    end
    subplot_tight(1,ns,i, [.1,.01])
    image(c_im)
    axis image off
    title(sprintf('sess: %d,  %d of %d', sessns(i), sum(shared), sum(cmap(:,i)>0)))
    
end
disp(sum(shared))
%%









function im = seg_im(cs, contour_thresh)

for k = 1:size(cs,1)
    ccs = squeeze(cs(k,:,:));
    ccs = ccs-min(ccs(:));
    ccs = ccs./max(ccs(:));
    ccs(ccs<=contour_thresh) = 0;
    cs(k,:,:) = ccs;
end
im = sum(cs,1).*1;

end