function [peth] = make_peth_TTL(data, inds, pwidth)
%%
n = length(inds);
peth = NaN(n, 2*pwidth + 1);
for i = 1:length(inds)
    a = inds(i)-pwidth;
    b = inds(i)+pwidth;
    % only accepts peths that fit
    if a>0 && b<=length(data)
        peth(i, :) = data(a:b);
    end
end
    