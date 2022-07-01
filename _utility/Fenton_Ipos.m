function [i]=Fenton_Ipos(spks, window_size, countmap, secmap)
s = normalize_rows(spks);
sumspks = bin_spks(s, window_size, false);
nspk = ceil(normalize_rows(sumspks)*4);
nspk(isnan(nspk))=0;

% k=find(c>0 & s>0); % find pixels in both maps with positive values
% 
% c=c(k); % keep positive values only
% s=s(k); % keep positive values only
% 
% fx=c./s;        % firing rate at each location 
% px=s/sum(s(:)); % probability of location
% 
% fmean=sum(px.*fx); % mean firing rate
% 
% infoRate=sum(fx.*px.*log2(fx))-fmean*log2(fmean);
% infoPerSpike=infoRate/fmean;