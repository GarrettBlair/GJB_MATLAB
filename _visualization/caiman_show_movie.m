function caiman_show_movie(spatial, activity, dims, bg_image)

% spatial = fullA;%spatial;
% activity = C(:,1:4:end);
% bg_image = minFrame./50;
Y = spatial*activity;

for i = 1:size(activity,2)
    f = reshape(Y(:,i), dims);
    imagesc(f + bg_image, [0 10])
    drawnow; 
end