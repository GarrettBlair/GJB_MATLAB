function [sfep, room_inferred_ipos, arena_inferred_ipos] = SFEP_map_quantification(ms, smoothing_size, plotting, verbose)
%% Generate SFEP maps for real data
r = ms.room; a = ms.arena;

ipos = abs(r.momentary_pos_info) - abs(a.momentary_pos_info);
ipos_mean= nanmean(ipos,1);

ipos_mean_zero = ipos_mean;
ipos_mean_zero(isnan(ipos_mean_zero)) = 0;

smoothing_kern = gausswin(smoothing_size*2+1);
smoothing_kern = smoothing_kern*smoothing_kern';
smoothing_kern = smoothing_kern./(sum(smoothing_kern(:)));

[room_sfep_map, room_sfep_map_smooth, room_sfep_coherence, room_sfep_bits, ~] = sfep_map(r.svm_decoding.x, r.svm_decoding.y,...
    ipos_mean_zero, ms.params.pos_bins, ms.params.pos_bins, r.vmap, smoothing_kern);
[arena_sfep_map, arena_sfep_map_smooth, arena_sfep_coherence, arena_sfep_bits, ~] = sfep_map(a.svm_decoding.x, a.svm_decoding.y,...
    ipos_mean_zero, ms.params.pos_bins, ms.params.pos_bins, a.vmap, smoothing_kern);

if verbose==true
fprintf('\nREAL ROOM bits:%2.2f   coh:%2.2f', room_sfep_bits, room_sfep_coherence)
fprintf('\nREAL ARENA bits:%2.2f   coh:%2.2f', arena_sfep_bits, arena_sfep_coherence)
end


%% Generate null ipos driven only by position
if ~isfield(ms.room, 'inferred_ipos') && ~isfield(ms.arena, 'inferred_ipos')
    if verbose==true
        fprintf('\n\t! Generating inferred ipos for comparison...')
    end
    spks = normalize_rows(ms.spks);
    % pfr=ms.room.pfields_smooth;
    % pfa=ms.arena.pfields_smooth;
    pfr=ms.room.pfields;
    pfa=ms.arena.pfields;
    
    spksroom  = NaN(size(spks));
    spksarena = NaN(size(spks));
    spksroom_sub  = NaN(size(spks));
    spksarena_sub = NaN(size(spks));
    
    [~, ~, ~, xbr, ybr] = histcounts2(ms.room.x, ms.room.y, ms.params.pos_bins, ms.params.pos_bins);
    [~, ~, ~, xba, yba] = histcounts2(ms.room.x, ms.room.y, ms.params.pos_bins, ms.params.pos_bins);
    xbr(xbr==0) = interp1(find(xbr~=0), xbr(xbr~=0), find(xbr==0), 'nearest');
    ybr(ybr==0) = interp1(find(ybr~=0), ybr(ybr~=0), find(ybr==0), 'nearest');
    xba(xba==0) = interp1(find(xba~=0), xba(xba~=0), find(xba==0), 'nearest');
    yba(yba==0) = interp1(find(yba~=0), yba(yba~=0), find(yba==0), 'nearest');
    % infer the firing rate from the place fields
    for i = 1:length(xbr)
        spksroom(:,i)  = pfr(:, ybr(i), xbr(i));
        spksarena(:,i) = pfa(:, yba(i), xba(i));
    end
    spksroom = normalize_rows(spksroom);
    spksarena = normalize_rows(spksarena);
    % equalize for the total number of spikes by downsampling (choosing the
    % n largest values)
    quant_level = 10; % how finely to quantize the data, higher more fine. should be >1 integer
    for i = 1:size(spksroom)
        s = spks(i,:);
        s(isnan(s)) = 0;
        ns = nansum(s>0);
        quant_s = round(s*quant_level);
        us = unique(quant_s(quant_s>0));
        
        aa = spksroom(i,:); 
        aa(isnan(aa)) = 0;
        quant_aa = round(aa*quant_level);
        aaa = aa*0; 
        
        bb = spksarena(i,:); 
        bb(isnan(bb)) = 0;
        quant_bb = round(bb*quant_level);
        bbb = bb*0; 
        
        for j = 1:length(us)
            nsu = sum(quant_s == us(j));
            xx = find(abs(quant_aa - us(j)) == min(abs(quant_aa - us(j))));
            while length(xx) < nsu
                quant_aa(xx) = NaN;
                tempxx = find(abs(quant_aa - us(j)) == min(abs(quant_aa - us(j))));
                xx = [xx, tempxx];
            end
            [~,ord] = sort(rand(size(xx)));
            aaa(xx(ord(1:nsu))) = us(j)/quant_level;

            xx = find(abs(quant_bb - us(j)) == min(abs(quant_bb - us(j))));
            while length(xx) < nsu
                quant_bb(xx) = NaN;
                tempxx = find(abs(quant_bb - us(j)) == min(abs(quant_bb - us(j))));
                xx = [xx, tempxx];
            end
            [~,ord] = sort(rand(size(xx)));
            bbb(xx(ord(1:nsu))) = us(j)/quant_level;
        end
        spksroom_sub(i,:) = aaa;
        spksarena_sub(i,:) = bbb;
%         figure(10);clf; hold on; plot(aaa+2); plot(bbb+1); plot(s)
%         title(sprintf('s=%3.1f   aaa=%3.1f   bbb=%3.1f\ns=%d   aaa=%d   bbb=%d', nansum(s), nansum(aaa), nansum(bbb), nansum(s>0), nansum(aaa>0), nansum(bbb>0)))
%         input('')
    end
    spksroom_sub = normalize_rows(spksroom_sub);
    spksarena_sub = normalize_rows(spksarena_sub);
    
    arenams = ms; arenams.spks = spksarena_sub;
    roomms = ms; roomms.spks = spksroom_sub;
    arenams.params.skip_ensemble = true;
    roomms.params.skip_ensemble = true;
    [arena_inferred_ipos] = Fenton_ipos(arenams, .25, 'arena', arenams.params);
    [room_inferred_ipos] = Fenton_ipos(roomms, .25, 'room', roomms.params);
    if verbose==true
        fprintf(' done!\n')
    end
    ms.arena.inferred_ipos = arena_inferred_ipos;
    ms.room.inferred_ipos = room_inferred_ipos;
else
    if verbose==true
        fprintf('\n\tUsing previous inferred ipos\n')
    end
    arena_inferred_ipos = ms.arena.inferred_ipos;
    room_inferred_ipos  = ms.room.inferred_ipos;
end
%% Generate SFEP maps for inferred null data
inferred_av_ipos = nanmean( abs(room_inferred_ipos) - abs(arena_inferred_ipos), 1);


inferred_av_ipos_zero = inferred_av_ipos;
inferred_av_ipos_zero(isnan(inferred_av_ipos_zero)) = 0;

[room_sfep_map_inf, room_sfep_map_inf_smooth, room_sfep_inf_coherence, room_sfep_inf_bits, ~] = sfep_map(r.svm_decoding.x, r.svm_decoding.y,...
    inferred_av_ipos_zero, ms.params.pos_bins, ms.params.pos_bins, r.vmap, smoothing_kern);

[arena_sfep_map_inf, arena_sfep_map_inf_smooth, arena_sfep_inf_coherence, arena_sfep_inf_bits, ~] = sfep_map(a.svm_decoding.x, a.svm_decoding.y,...
    inferred_av_ipos_zero, ms.params.pos_bins, ms.params.pos_bins, a.vmap, smoothing_kern);

if verbose==true
fprintf('\nINFER ROOM bits:%2.2f   coh:%2.2f', room_sfep_inf_bits, room_sfep_inf_coherence)
fprintf('\nINFER ARENA bits:%2.2f   coh:%2.2f', arena_sfep_inf_bits, arena_sfep_inf_coherence)
end

if plotting==true
% figure(9400); clf; 
% hold on
% plot(inferred_av_ipos, 'k');
% axis tight
% % ylim([-2.5 2.5])
% % yyaxis('right')
% plot(ipos_mean, 'b');
% axis tight

figure(9399); clf;   colormap magma
subplot_tight(2,3,1); i = room_sfep_map_smooth; 
ii = imagesc(i, [0 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('room sfep');
subplot_tight(2,3,2); i = room_sfep_map_inf_smooth;
ii = imagesc(i, [0 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('room control')
subplot_tight(2,3,3); i = room_sfep_map_smooth - room_sfep_map_inf_smooth;
ii = imagesc(i, [-1 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('room sfep diff')

subplot_tight(2,3,4); i = arena_sfep_map_smooth;
ii = imagesc(i, [0 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('arena sfep')
subplot_tight(2,3,5); i = arena_sfep_map_inf_smooth;
ii = imagesc(i, [0 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('arena control')
subplot_tight(2,3,6); i = arena_sfep_map_smooth - arena_sfep_map_inf_smooth;
ii = imagesc(i, [-1 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('arena sfep diff')

figure(9398); clf;  colormap magma
subplot_tight(2,3,1); i = room_sfep_map ;
ii = imagesc(i, [0 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('room sfep'); 
subplot_tight(2,3,2); i = room_sfep_map_inf ;
ii = imagesc(i, [0 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('room control')
subplot_tight(2,3,3); i = room_sfep_map  - room_sfep_map_inf ;
ii = imagesc(i, [-1 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('room sfep diff')

subplot_tight(2,3,4); i = arena_sfep_map ;
ii = imagesc(i, [0 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('arena sfep')
subplot_tight(2,3,5); i = arena_sfep_map_inf ;
ii = imagesc(i, [0 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('arena control')
subplot_tight(2,3,6); i = arena_sfep_map  - arena_sfep_map_inf ;
ii = imagesc(i, [-1 1]); axis image off; ii.AlphaData = double(~isnan(i)); set(gca, 'YDir', 'normal'); axis off image; colorbar; title('arena sfep diff')
end

% room_diff = nanmean(nanmean(room_sfep_map  - room_sfep_map_inf));
% arena_diff = nanmean(nanmean(arena_sfep_map - arena_sfep_map_inf));
% room_absdiff = nanmean(nanmean(abs(room_sfep_map  - room_sfep_map_inf)));
% arena_absdiff = nanmean(nanmean(abs(arena_sfep_map - arena_sfep_map_inf)));
room_diff = nanmean(nanmean(room_sfep_map_smooth  - room_sfep_map_inf_smooth));
arena_diff = nanmean(nanmean(arena_sfep_map_smooth - arena_sfep_map_inf_smooth));
room_absdiff = nanmean(nanmean(abs(room_sfep_map_smooth  - room_sfep_map_inf_smooth)));
arena_absdiff = nanmean(nanmean(abs(arena_sfep_map_smooth - arena_sfep_map_inf_smooth)));
% room_diff = nanmean(nanmean(room_sfep_map  - room_sfep_map_inf));
% arena_diff = nanmean(nanmean(arena_sfep_map - arena_sfep_map_inf));
% room_absdiff = nanmean(nanmean(abs(room_sfep_map  - room_sfep_map_inf)));
% arena_absdiff = nanmean(nanmean(abs(arena_sfep_map - arena_sfep_map_inf)));
room_corr = nancorr(room_sfep_map_smooth(:), room_sfep_map_inf_smooth(:));
arena_corr = nancorr(arena_sfep_map_smooth(:), arena_sfep_map_inf_smooth(:));
% room_corr = nancorr(room_sfep_map(:), room_sfep_map_inf(:));
% arena_corr= nancorr(arena_sfep_map(:), arena_sfep_map_inf(:));



sfep = [];
% maps
sfep.room_vmap              = r.vmap;
sfep.arena_vmap             = a.vmap;
sfep.room_map               = room_sfep_map;
sfep.room_map_smooth        = room_sfep_map_smooth;
sfep.arena_map              = arena_sfep_map;
sfep.arena_map_smooth       = arena_sfep_map_smooth;
sfep.infer.room_map         = room_sfep_map_inf;
sfep.infer.room_map_smooth  = room_sfep_map_inf_smooth;
sfep.infer.arena_map        = arena_sfep_map_inf;
sfep.infer.arena_map_smooth = arena_sfep_map_inf_smooth;

% metrics
% corr doesn't appear to work well
% [sfep.room_corr(1), sfep.room_corr(2)] = nancorr(room_sfep_map_smooth(:), room_sfep_map_inf_smooth(:));
% [sfep.arena_corr(1), sfep.arena_corr(2)] = nancorr(arena_sfep_map_smooth(:), arena_sfep_map_inf_smooth(:));
sfep.room_diff              = room_diff;
sfep.arena_diff             = arena_diff;
sfep.room_absdiff           = room_absdiff;
sfep.arena_absdiff          = arena_absdiff;
sfep.room_corr              = room_corr;
sfep.arena_corr             = arena_corr;

if verbose==true
fprintf('\n\tROOM diff: %2.2f', room_diff)
fprintf('\n\tARENA diff: %2.2f', arena_diff)
fprintf('\n\n')
end

sfep.room_coh               = room_sfep_coherence;
sfep.room_bits              = room_sfep_bits;
sfep.arena_coh              = arena_sfep_coherence;
sfep.arena_bits             = arena_sfep_bits;
sfep.infer.room_coh         = room_sfep_inf_coherence;
sfep.infer.room_bits        = room_sfep_inf_bits;
sfep.infer.arena_coh        = arena_sfep_inf_coherence;
sfep.infer.arena_bits       = arena_sfep_inf_bits;


end
%%
function [sfep_map, sfep_map_smooth, sfep_coherence, sfep_bits, sfep_bitsRate] = sfep_map(x,y,sfepsig, xbins, ybins, vmap, smoothing_kern)
% [sfep_map, ~] = make_occupancymap_2D(x, y, sfepsig, xbins, ybins);
[sfep_map, map_counts] = sfep_prob_map_2D(x, y, sfepsig, xbins, ybins);
% sfep_map = sfep_map./vmap;
sfep_map = sfep_map;
sfep_map_smooth = sfep_map;
sfep_map_smooth(isnan(sfep_map_smooth))=0;
sfep_map_smooth = conv2(sfep_map_smooth, smoothing_kern, 'same');
sfep_map_smooth(isnan(vmap)) = NaN;

sfep_map(isnan(vmap)) = 0;
[~, coherence] = pfield_coherence_calc(permute(cat(3, sfep_map, sfep_map), [3,1,2]), vmap);
sfep_coherence = coherence(1);
[sfep_bits, sfep_bitsRate] = infoSkaggs(normalize_matrix(sfep_map), vmap);
sfep_map(isnan(vmap)) = NaN;

sfep_map = normalize_matrix(sfep_map);
sfep_map_smooth = normalize_matrix(sfep_map_smooth);

end
