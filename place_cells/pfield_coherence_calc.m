function [coh, coh_z] = pfield_coherence_calc(pfields, occupancy)

[nsegs, nbins2, nbins1] = size(pfields);
neighbor_av_r = pfields*NaN;
for c1 = 1:nbins1
    for c2 = 1:nbins2
        rind = c1+[-1  0  1 -1 1 -1 0 1];
        rind = rind(rind>0 & rind<=nbins1);
        cind = c2+[-1 -1 -1  0 0  1 1 1];
        cind = cind(cind>0 & cind<=nbins2);
%         b = nanmean(nanmean(squeeze(pfields(:,cind,rind)), 2), 3);
        b = nanmean(nanmean((pfields(:,cind,rind)), 2), 3);
        neighbor_av_r(:,c2,c1) = b;
    end
end
%%
coh = NaN(nsegs,1);
for i = 1:nsegs
    p1 = squeeze(pfields(i,:,:));
    p2 = squeeze(neighbor_av_r(i,:,:));
    p2(isnan(occupancy)) = NaN;
%     figure(8675309);
%     subplot(121); imagesc(p1);
%     subplot(122); imagesc(p2)
%     drawnow
    pcheck = p1.*p2;
    valid = ~isnan(pcheck(:));
    if any(valid(:))
        coh(i) = corr(p1(valid), p2(valid));
    else
        drawnow
    end
    
end
% Z transform r values (Spearman to Fisher) : https://www.statisticshowto.com/fisher-z/
coh_z = .5*( log(1+coh) - log(1-coh) );
