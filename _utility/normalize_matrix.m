function [mat] = normalize_matrix(mat)

if ismatrix(mat)
    min_mat = min(mat(:));
    mat = mat - min_mat;
    max_mat = max(mat(:));
    mat = mat./max_mat;
    mat = round(mat*10000)/10000;
elseif ndims(mat)==3
    for i = 1:size(mat,1)
        m = mat(i,:,:);
        min_mat = min(m(:));
        m = m - min_mat;
        max_mat = max(m(:));
        m = m./max_mat;
        mat(i,:,:) = round(m*10000)/10000;
    end
end

end