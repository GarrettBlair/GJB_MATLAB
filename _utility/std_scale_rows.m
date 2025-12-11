function [mat] = std_scale_rows(mat)
dim = 2;
mean_mat = nanmean(mat, dim)*ones(1, size(mat,dim));
mat = mat-mean_mat;

min_mat = nanstd(mat, [], dim)*ones(1, size(mat,dim));
mat = mat./min_mat;
% mat = sqrt(mat);%./min_mat;
% max_mat = max(mat, [], dim)*ones(1, size(mat,dim));
% mat = mat./max_mat;
% mat = round(mat*10000)/10000;

end
