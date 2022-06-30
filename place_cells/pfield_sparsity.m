function sparsity = pfield_sparsity(pfields)
%% Calculates the sparsity of a 1D place field, ranges from 0 (uniform) to 1 (impulse-like)
% a la Ravassard, Mehta et al. 2013
map_dims = ndims(pfields);
if map_dims == 2
%     map_type = '1D';
    [nsegs, nbins1] = size(pfields);
    pfields_alt = pfields;
elseif map_dims == 3
%     map_type = '2D';
    [nsegs, nbins1, nbins2] = size(pfields);
    pfields_alt = reshape(pfields, [nsegs, nbins1*nbins2]);
else
    error('Unkown map type!')
%     return
end
nanbins = squeeze(sum(isnan(pfields), 1))==nsegs; % find bins that are always nan, i.e. never visited

pfields_alt(isnan(pfields_alt)) = 0;
sparsity = NaN(size(pfields_alt, 1), 1);
nbins = sum(~nanbins(:)); % size(pfields_alt, 2);
for i = 1:nsegs
    p1 = pfields_alt(i,:);
    a = sum(p1).^2;
    b = nbins*sum(p1.^2);
    sparsity(i) = (1 - a./b)*(nbins/(nbins-1));
    try 
        p0 = squeeze(pfields(i,:,:));
    catch
        p0 = squeeze(pfields(i,:));
    end
%     subplot(1,2,1)
%     imagesc(p0)
%     subplot(122)
%     plot(p1)
%     title(sparsity(i))
%     drawnow
%     input('');
end