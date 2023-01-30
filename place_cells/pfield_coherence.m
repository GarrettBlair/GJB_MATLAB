%%
% load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files\2022_09_14_H17_52_22_WTR8_@placecells.mat')
pfields_r = ms.room.pfields;
pfields_a = ms.room.pfields;
[nsegs, nbins1, nbins2] = size(pfields_r);
coh_room = NaN(nsegs,1);
coh_arena= NaN(nsegs,1);
rad = 1;
neighbor_av_r = pfields_r*NaN;
neighbor_av_a = pfields_a*NaN;
for c1 = 2:nbins1-1
    for c2 = 2:nbins2-1
%     a = squeeze(pfields(:,c1,c2));
    rind = c1+[-1  0  1 -1 1 -1 0 1];
    cind = c2+[-1 -1 -1  0 0  1 1 1];
    b = nanmean(nanmean(squeeze(pfields_r(:,cind,rind)), 2), 3);
    neighbor_av_r(:,c2,c1) = b;
    
    b = nanmean(nanmean(squeeze(pfields_a(:,cind,rind)), 2), 3);
    neighbor_av_a(:,c2,c1) = b;
    end
end
%%
for i = 1:nsegs
            p1 = squeeze(pfields_r(i,:,:));
            p2 = squeeze(neighbor_av_r(i,:,:));
            pcheck = p1.*p2;
            valid = ~isnan(pcheck(:));
            coh_room(i) = corr(p1(valid), p2(valid));
            
%             figure(2); clf
%             subplot(221)
%             imagesc(p1);
%             title(sprintf('Room coh- %0.2f', coh_room(i)))
%             subplot(222)
%             imagesc(p2);
            p1 = squeeze(pfields_a(i,:,:));
            p2 = squeeze(neighbor_av_a(i,:,:));
            pcheck = p1.*p2;
            valid = ~isnan(pcheck(:));
            coh_arena(i) = corr(p1(valid), p2(valid));
            
%             subplot(223)
%             imagesc(p1);
%             title(sprintf('Arena coh- %0.2f', coh_arena(i)))
%             subplot(224)
%             imagesc(p2);
%             drawnow
%             pause(1)
            
end
% Z transform r values (Spearman to Fisher) : https://www.statisticshowto.com/fisher-z/
coh_arena_z = .5*( log(1+coh_arena) - log(1-coh_arena) );
coh_room_z = .5*( log(1+coh_room) - log(1-coh_room) );
