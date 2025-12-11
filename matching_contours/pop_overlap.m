function [overlap_mat, lin_overlap] = pop_overlap(cellmap)
%%
[~, nsess] = size(cellmap);
overlap_mat = NaN(nsess);

for i = 1:nsess-1
    for j = i+1:nsess
        a = sum(cellmap(:,i)>0 & cellmap(:,j)>0);
        b = sum(cellmap(:,i)>0 | cellmap(:,j)>0);
        overlap_mat(i,j) = a/b;
        overlap_mat(j,i) = a/b;
    end
end

lin_overlap = NaN(nsess, nsess);
for sep = 1:nsess-1
    idx = 0;
    for i = 1:nsess-sep
        for j = i+sep
            idx = idx+1;
            lin_overlap(sep, idx) = overlap_mat(i,j);
        end
    end
end