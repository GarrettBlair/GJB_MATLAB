
function vmap = make_vmap_2D(posx,posy,bins,nbins,minsamples)

%posx,posy = x,y values of position samples
%nbins = number of spatial bins along each dimension of the spatial map (x and y dimension have same bin resolution)
%minsamples = occupancy inclusion threshold (in samples) for spatial bins; 
%          bins occuppied for less than mintime will be excluded from the rate map

%vmap = speed-filtered occupancy map

% % % spbinsize = 100/nbins;
% % % vmap = NaN(nbins); %visit map
% % % posx = 100*((posx - min(posx))/(2*max(posx))); % normalize position
% % % posy = 100*((posy - min(posy))/(2*max(posy)));
% vmap = NaN(nbins,nbins);
% vmap = zeros(nbins,nbins);
% 
% [xcount, xbin] = histc(posx,bins); 
% for x=1:nbins
%     if ~isempty(posy(find(x==xbin)))
%         [ycount, ybins] = histc(posy(find(x==xbin)),bins);
%         vmap(:,x) = ycount(nbins:-1:1);
%     end
% end
vmap = histcounts2(posx, posy, bins, bins);
vmap(vmap<minsamples)=NaN;
vmap = vmap/(nansum(vmap(:)));
% vmap(find(vmap(:)<minsamples))=0;
