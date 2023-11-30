function splits = splits_equal_occupancy(var_binned, counts, nsplits)
splits = NaN(length(var_binned),1);
odd_ind = [];
split_equal = 1;
% ux = unique(counts);
for i = 1:length(counts)
    %     ind = find(var_binned==ux(i));
    ind = find(var_binned==i);
    %         [~, ranord] = sort(rand(length(ind),1));
    %             ind = ind(ranord);
    if mod(length(ind), 2) == 1
        % this is to keep the last split value from always geting one extra
        % when there's an odd number
        odd_ind = ind(end);
        ind = ind(1:end-1);
    else
        odd_ind = [];
    end
    if ~isempty(ind)
        xi = floor(linspace(0, length(ind), nsplits+1));
        for j = 1:nsplits
            splits(ind(xi(j)+1:xi(j+1))) = j;
        end
    end
    if ~isempty(odd_ind)
        splits(odd_ind) = split_equal;
        split_equal = mod(split_equal, nsplits)+1;
    end
end