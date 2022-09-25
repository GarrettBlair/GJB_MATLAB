function [scrambled_mat] = scramble_mat(mat, dim)
if ~ismatrix(mat)
    error('matrix arg (#1) must be of size <= 2!');
end
n = size(mat, dim);
[~, ord] = sort(rand(n,1));
if dim == 1
scrambled_mat = mat(ord, :);
elseif dim == 2
scrambled_mat = mat(:, ord);   
else
    error('dimension arg (#2) must be 1 or 2!');
end

