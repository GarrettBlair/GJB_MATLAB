function [corr_val, pval] = nancorr(a,b)
% performs a corr on two signals after removing any row with nans between either
if nargin==2
    a1 = NaN(length(a), 2);
    a1(:,1) = a;
    a1(:,2) = b;
    
    nanrow = any(isnan(a1), 2);
    a1 = a1(~nanrow, :);
    
    [corr_val, pval] = corr(a1(:,1), a1(:,2));
elseif nargin==1
    if size(a,2)<=1
        error('First arg must have more than 1 column')
    else
        nanrow = any(isnan(a), 2);
        a1 = a(~nanrow, :);
        [corr_val, pval] = corr(a1);
    end
end