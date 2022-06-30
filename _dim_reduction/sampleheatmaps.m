figure(1);
create_heatmap(ncells, nbins, sparsity)

function create_heatmap(ncells, nbins, sparsity)

M = zeros(ncells, nbins);
centerbin = randi(nbins, [ncells, 1]);

for i = 1:ncells
    v = zeros(1, nbins);
    v(centerbin(i)) = 1;
    n = (rand(1, nbins))/sparsity;
M(i, :) = v+n;  
end

end

% function decimate_heatmap(M, ndec)
% for d = 1:ndec
%     for i = 1:ncells
%         
%     end
% end
% 
% end


function normsort_plot(M)
[mmax,pos] = max(M, [], 2);
a =  ones(1, size(M,2));
mnorm = mmax*a;

MM = M./mnorm;
[~,pos] = sort(pos);

MM = MM(pos,:);
imagesc(MM)
end