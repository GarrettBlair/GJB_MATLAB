
function smap = make_Ca_map_1D(posx,bin_edges,Cavec)

%nbins = number of spatial bins along each dimension of the spatial map (x and y dimension have same bin resolution)
%mapspk = spike time stamps (these are assumed to already be speed filtered
%posx,posy = x,y values of position data
%TS = timestamps of positobn data

[~, bin] = histc(posx,bin_edges); 
smap = zeros(1, length(bin_edges)-1);
for i = 1:length(bin_edges)-1
        smap(i) = sum(Cavec(bin == i));
end
