
function smap = make_smap_1D(behspk,binedges)

%compute spike count map (histogram of spike counts at binned values) for a 1D behavioral variable

%behspk = values of sampled behavioral variable at each spike occurrence
%binedges = vector of bin edges for the spike count map

%smap = spike count map

smap = zeros(1,length(binedges)-1); %spike counts

smap = histc(behspk,binedges); 

smap = smap(1:(end-1));



