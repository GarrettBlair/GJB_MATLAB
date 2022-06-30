
function [vmap, idex] = make_vmap_1D(behvar,TS,binedges,intvar)

%compute occupancy map (histogram of time spent at binned values) for a 1D behavioral variable

%behvar = values of sampled behavioral variable
%TS = timemstamps of each value in behvar
%binedges = vector of bin edges for the visit map
%intvar = interval variable for filtering occupancy times

%vmap = occupancy map
%indices of elements in behvar that passed through the interval filter = occupancy map

vmap = zeros(1,length(binedges)-1); %visit map
idex = [];

for i=1:size(intvar(:,1)) %loop through intervals
    
    rdex = find((TS>=intvar(i,1)) & (TS<=intvar(i,2)));
    if ~isempty(rdex)
        [xcount, xbin] = histc(behvar(rdex),binedges); 
        vmap = vmap + (xcount(1:(end-1)))';
    end
    if length(intvar(:,1))>1 %only return speed filtered values if there is more than one interval
        idex = [idex; rdex];
    end
end

vmap(vmap == 0) = NaN;