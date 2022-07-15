%% function to evaluate each cell grpahically
nc = 5; nr = 4;

[nsegs, nframes] = size(ms.neuron.C);

contours = gbContours(ms.neuron.A, ms.neuron.dims, [], .6);
craw = normalize_rows(ms.neuron.C+ms.neuron.YrA);
% spks = squeeze(sum(ms.neuron.S_matw,1));
spks = normalize_rows(ms.neuron.S_matw);
caiman_c = normalize_rows(ms.neuron.C);
win = 50;
dims = double(ms.neuron.dims);
t = ms.timestamps./1000;
%%
kern = gb_kernel(4, 'step3');
cspks = spks;
for j = 1:1
for i = 1:nsegs
    cspks(i,:) = conv(cspks(i,:), kern, 'same');
end
end
cspks = normalize_rows(cspks);

spd = ms.speed_epochs';
figure(1)
set(gcf, 'Position', [558         137        1133         782]); clf;
r = [1  .4  .4];
g = [.2 .8 .2];
b = [.4 .4  1];
for i = 1:nsegs
    %%
    figure(1); clf
    
    c = squeeze(contours(i,:,:)>0)/4;
    c0 = squeeze(sum(contours(setdiff(1:nsegs, i),:,:)>0,1))/4;
    bg = ms.neuron.meanFrame./500;
    cim = cat(3, c0+bg, c+bg, bg);
    [~, rx] = max(sum(c,1));
    rxx = rx-win:rx+win;
    rxx = rxx(rxx>0);
    rxx = rxx(rxx<=dims(2));
    [~, ry] = max(sum(c,2));
    ryy = ry-win:ry+win;
    ryy = ryy(ryy>0);
    ryy = ryy(ryy<=dims(1));
    cim_sub = cim(ryy, rxx, :);
    cim_sub = cim_sub./max(max(max(cim_sub)));
    
    t_raw  = craw(i,:);
    t_spks = spks(i,:);
    t_c = cspks(i,:);
    spdspk = t_spks>0 & spd;
    stlspk = t_spks>0 & ~spd;
%     s = t_spks>0;
    spd_size = round(t_spks(spdspk)*100);
    stl_size = round(t_spks(stlspk)*100);
    
    room_map    = squeeze(ms.room.pfields(i,:,:));
    room_map_s  = squeeze(ms.room.pfields_smooth(i,:,:));
    arena_map    = squeeze(ms.arena.pfields(i,:,:));
    arena_map_s  = squeeze(ms.arena.pfields_smooth(i,:,:));
    
%     pnr = ms.neuron.pnr(ryy, rxx);
    pnr = ms.neuron.pnr;
%     pnr = pnr./max(pnr(:))/2;
%     cim2 = cat(3, c0+pnr, c+pnr, pnr);
%     cim2_sub = cim2(ryy, rxx, :);
%     cim2_sub = cim2_sub./max(max(max(cim2_sub)));
%     pnr = pnr./max(pnr);
    %   CONTOUR   %
    subplot_tight(nr, nc, 1) 
    image(cim_sub)
    axis image off
    %   ROOM VMAP   %
    subplot_tight(nr, nc, 2);
    imagesc(ms.room.vmap, 'AlphaData', ms.room.pfield_alpha)
    axis image off
    %   ROOM POS AND SPIKES   %
    
    subplot_tight(nr, nc, 3); hold on
    plot(ms.room.x, ms.room.y*-1, 'Color', r)
    scatter(ms.room.x(stlspk), ms.room.y(stlspk)*-1, stl_size, 'MarkerFaceColor', [.4 .4 .4], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75)
    scatter(ms.room.x(spdspk), ms.room.y(spdspk)*-1, spd_size, 'MarkerFaceColor', r/2, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .85)
    axis square
    %   ROOM PMAP   %
    subplot_tight(nr, nc, 4);
    imagesc(room_map, 'AlphaData', ms.room.pfield_alpha)
    axis image off
    %   ROOM PMAP SMOOTH   %
    subplot_tight(nr, nc, 5);
    imagesc(room_map_s, 'AlphaData', ms.room.pfield_alpha)
    axis image off
    %   CONTOUR   %
    subplot_tight(nr, nc, 1+nc) 
    imagesc(pnr)
%     image(cim2_sub)
%     colorbar
    axis image off
    %   ARENA VMAP   %
    subplot_tight(nr, nc, 2+nc); hold on
    imagesc(ms.arena.vmap, 'AlphaData', ms.arena.pfield_alpha)
    axis image off
    %   ARENA POS AND SPIKES  %
    subplot_tight(nr, nc, 3+nc); hold on
    plot(ms.arena.x, ms.arena.y*-1, 'Color', b)
%     scatter(ms.arena.x(s), ms.arena.y(s)*-1, s_size, 'MarkerFaceColor', b/3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75)
    scatter(ms.arena.x(stlspk), ms.arena.y(stlspk)*-1, stl_size, 'MarkerFaceColor', [.4 .4 .4], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75)
    scatter(ms.arena.x(spdspk), ms.arena.y(spdspk)*-1, spd_size, 'MarkerFaceColor', b/2, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .85)
    axis square
    %   ARENA PMAP  %
    subplot_tight(nr, nc, 4+nc);
    imagesc(arena_map, 'AlphaData', ms.arena.pfield_alpha)
    axis image off
    
    %   ARENA PMAP SMOOTH   %
    subplot_tight(nr, nc, 5+nc);
    imagesc(arena_map_s, 'AlphaData', ms.arena.pfield_alpha)
    axis image off
    
    %   SPEED EPOCHS   %
%     subplot_tight(nr, 1, nr-1)
%     hold on
%     plot(ms.speed_epochs*80, 'Color', g, 'LineWidth', 2)
%     plot(ms.arena.speed, 'k')
%     axis tight
%     ylim([0 80])
    
    %   TRACE AND SPIKES   %
    subplot_tight(nr, 1, nr-1:nr)
    hold on
    colormap viridis
%     imagesc([cat(1, spd, spd)], -)
    imagesc([cat(1, ms.arena.speed_smooth', ms.arena.speed_smooth', ms.arena.speed_smooth')], [0 20])
%     t_c = t_c + 1;
    plot(1.5*(t_c+.75), 'Color', 'r', 'LineWidth', 2)
%     plot(t_spks*1.5, 'Color', r, 'LineWidth', 2)
    plot(t_raw+1, 'Color', [.8 .8 .8])
    axis tight off
    colorbar
    text(1.1*nframes, 2.5, 'speed (clipped)', 'Rotation', 270)
%     input('')
end