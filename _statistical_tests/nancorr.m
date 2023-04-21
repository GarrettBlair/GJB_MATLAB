function [corr_val, pval] = nancorr(a,b)
a1 = NaN(length(a), 2);
a1(:,1) = a;
a1(:,2) = b;

nanind = any(isnan(a1), 2);
a1 = a1(~nanind, :);

[corr_val, pval] = corr(a1(:,1), a1(:,2));