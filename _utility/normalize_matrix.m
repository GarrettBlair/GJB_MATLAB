function [mat] = normalize_matrix(mat)

min_mat = min(mat(:));
mat = mat - min_mat;
max_mat = max(mat(:));
mat = mat./max_mat;
mat = round(mat*10000)/10000;

end