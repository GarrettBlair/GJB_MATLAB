function [jitter] = gb_rand_jitter(y, s)
%%
% make a distribution of random points based on each column data.
% data in col y will get greater variance the closer it is to the mean to
% aid in visualization

if isempty(s)
    s = 1;
end
if isrow(y)
    y = y';
end
jitter = NaN(size(y));
for i = 1:size(y,2)
    ysub = y(:,i);
    j_cent = rand(size(ysub))-.5;
    
%     goodind = ~isnan(ysub);
    meanx = nanmean(ysub);
    devx = abs(ysub-meanx);
    devx - devx-min(devx);
    devx = 1 - devx./max(devx);
%     j = 1 - ((devx-meanx)/maxx)
    jitter = (j_cent.*devx)./s;
    
    
%     
%     devx = (devx-min(devx))./(1+max(devx-min(devx)));
%     % devx = 1- devx./(max(devx)+1);
%     % jitter = j_cent./(s./(devx + .001*nanmin(x)));
%     jitter(:, i) = j_cent./((devx.^3 + 1)*s);
end
