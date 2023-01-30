function [overlap_mat] = pop_overlap(cellmap)
%%
[~, nsess] = size(cellmap);
overlap_mat = NaN(nsess);

for i = 1:nsess
    for j = 1:nsess
        if i~=j
        a = sum(cellmap(:,i)>0 & cellmap(:,j)>0);
        b = sum(cellmap(:,i)>0 | cellmap(:,j)>0);
        overlap_mat(i,j) = a/b;
        end
    end
end