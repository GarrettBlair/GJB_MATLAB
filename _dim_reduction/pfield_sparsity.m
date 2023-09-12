function sparsity = pfield_sparsity(pfields, occupancy)
%% Calculates the sparsity of a 1D place field, ranges from 0 (uniform) to 1 (impulse-like)
% a la Ravassard, Mehta et al. 2013 by GJB
pfields(isnan(pfields)) = 0;
sparsity = NaN(size(pfields, 1), 1);
for i = 1:size(pfields)
    if ndims(pfields)==2
    p1 = squeeze(pfields(i,:));
%     L = size(pfields, 2);
    L = sum(occupancy(:)>0);
    elseif ndims(pfields)==3
    p1 = squeeze(pfields(i,:,:));
    p1 = p1(:);
    L = sum(occupancy(:)>0); % size(pfields, 2) * size(pfields, 3);
    end
    a = sum(p1).^2;
    b = L*sum(p1.^2);
    sparsity(i) = (1 - a./b)*(L/(L-1));
%     plot(p1)
%     title(sparsity(i))
%     drawnow
%     input('');
end