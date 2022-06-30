function [infoPerSpike, infoRate]=infoSkaggs(countmap, secmap)

% infoSkaggs computes Shannon mutual information rates according to Skaggs 1992.
%
% Usage:
%   infoPerSpike=infoSkaggs(countmap, secmap)
%   [infoPerSpike, infoRate]=infoSkaggs(countmap, secmap)
%
% Input:
%   countmap is the total spike counts at various locations.
%   secmap is the total time (sec) spent at various locations.
%
% Output:
%   infoPerSpike is the information (bits) per spike.
%   infoRate is the information (bits) per sec.

c=countmap;
s=secmap;

if ~(size(c,1)==size(s,1))
 c=c';
end

k=find(c>0 & s>0); % find pixels in both maps with positive values

c=c(k); % keep positive values only
s=s(k); % keep positive values only

fx=c./s;        % firing rate at each location 
px=s/sum(s(:)); % probability of location

fmean=sum(px.*fx); % mean firing rate

infoRate=sum(fx.*px.*log2(fx))-fmean*log2(fmean);
infoPerSpike=infoRate/fmean;