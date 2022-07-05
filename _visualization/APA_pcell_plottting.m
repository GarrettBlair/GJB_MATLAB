%% function to evaluate each cell grpahically
nc = 5; nr = 4;

[nsegs, nframes] = size(ms.neuron.C);

contours = gbContours(ms.neuron.A, ms.neuron.dims, [], .6);
craw = normalize_rows(ms.neuron.C+ms.neuron.YrA);
spks = normalize_rows(ms.neuron.S);
caiman_c = normalize_rows(ms.neuron.C);
win = 50;
dims = double(ms.neuron.dims);
t = ms.timestamps./1000;

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
    t_c = caiman_c(i,:);
    s = t_spks>0;
    s_size = round(t_spks(s)*100);
    
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
    scatter(ms.room.x(s), ms.room.y(s)*-1, s_size, 'MarkerFaceColor', r/3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75)
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
    scatter(ms.arena.x(s), ms.arena.y(s)*-1, s_size, 'MarkerFaceColor', b/3, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75)
    %   ARENA PMAP  %
    subplot_tight(nr, nc, 4+nc);
    imagesc(arena_map, 'AlphaData', ms.arena.pfield_alpha)
    axis image off
    
    %   ARENA PMAP SMOOTH   %
    subplot_tight(nr, nc, 5+nc);
    imagesc(arena_map_s, 'AlphaData', ms.arena.pfield_alpha)
    axis image off
    
    %   SPEED EPOCHS   %
    subplot_tight(nr, 1, nr-1)
    hold on
    plot(ms.speed_epochs*80, 'Color', g, 'LineWidth', 2)
    plot(ms.arena.speed, 'k')
    axis tight
    ylim([0 80])
    
    %   TRACE AND SPIKES   %
    subplot_tight(nr, 1, nr)
    hold on
    plot(t_c*1.75, 'Color', g, 'LineWidth', 2)
%     plot(t_spks*1.5, 'Color', r, 'LineWidth', 2)
    plot(t_raw, 'k')
    axis tight
    input('')
end