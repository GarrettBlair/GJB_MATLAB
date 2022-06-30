function sparsity = pfield_sparsity(pfields)
%% Calculates the sparsity of a 1D place field, ranges from 0 (uniform) to 1 (impulse-like)
% a la Ravassard, Mehta et al. 2013 by GJB
pfields(isnan(pfields)) = 0;
sparsity = NaN(size(pfields, 1), 1);
L = size(pfields, 2);
for i = 1:size(pfields)
    p1 = pfields(i,:);
    
    a = sum(p1).^2;
    b = L*sum(p1.^2);
    sparsity(i) = (1 - a./b)*(L/(L-1));
%     plot(p1)
%     title(sparsity(i))
%     drawnow
%     input('');
end