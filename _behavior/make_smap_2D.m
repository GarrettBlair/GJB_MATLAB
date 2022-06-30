
function smap = make_smap_2D(posx,posy,bins,nbins,mapspk,TS)

%nbins = number of spatial bins along each dimension of the spatial map (x and y dimension have same bin resolution)
%mapspk = spike time stamps (these are assumed to already be speed filtered
%posx,posy = x,y values of position data
%TS = timestamps of positobn data

%smap = spike count map
% % % posx = 100*((posx - min(posx))/(2*max(posx)));
% % % posy = 100*((posy - min(posy))/(2*max(posy)));
% % % 
spkx = interp1(TS, posx, mapspk);
spky = interp1(TS, posy, mapspk);
        
spbinsize = 100/nbins;
smap = zeros(nbins); %spike counts

[xcount, xbin] = histc(spkx,bins); 
for i=1:nbins
    if ~isempty(spky(find(i==xbin)))
        [ycount, ybins] = histc(spky(find(i==xbin)),bins);
        smap(:,i) = ycount(nbins:-1:1);
    else
        smap(i,1:nbins) = 0;
    end
end


