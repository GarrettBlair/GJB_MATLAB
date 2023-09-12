function [inds] = time_to_indx(t, ts)
inds = NaN(size(ts));

for i = 1:length(ts)
    inds(i) = find(min(abs(t-ts(i))) == abs(t-ts(i)), 1, 'first');
end