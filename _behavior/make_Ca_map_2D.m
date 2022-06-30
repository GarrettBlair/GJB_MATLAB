
function smap = make_Ca_map_2D(posx,posy,bin_edges,Cavec)

%nbins = number of spatial bins along each dimension of the spatial map (x and y dimension have same bin resolution)
%mapspk = spike time stamps (these are assumed to already be speed filtered
%posx,posy = x,y values of position data
%Don't use speed filtered, but do use speed filtered traces
%TS = timestamps of positobn data

if length(posx) ~= length(Cavec)
    fprintf('Dimensionality missmatch!\n')
    return
end
% % [~, xbin] = histc(posx,bin_edges); 
% % [~, ybin] = histc(posy,bin_edges); 
% % smap = zeros(length(bin_edges)-1);
% % for i = 1:length(bin_edges)-1
% %     for j = 1:length(bin_edges)-1
% %         smap(j,i) = sum(Cavec(xbin == i & ybin == j));
% %     end
% % end
% % smap = smap(end:-1:1, :);
[~, ~, ~, xbin, ybin] = histcounts2(posx, posy, bin_edges, bin_edges); 
smap = zeros(length(bin_edges)-1);
for i = 1:length(bin_edges)-1
    for j = 1:length(bin_edges)-1
        smap(j,i) = sum(Cavec(xbin == i & ybin == j));
    end
end
