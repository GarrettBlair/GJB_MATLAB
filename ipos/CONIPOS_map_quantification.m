function [sfep, room_inferred_ipos, arena_inferred_ipos] = CONIPOS_map_quantification(ms, smoothing_size, plotting, verbose)
%% Generate SFEP maps for real data
r = ms.room; a = ms.arena;

ipos_mean = abs(r.conjoint_ipos_av) - abs(a.conjoint_ipos_av);
% ipos_mean= nanmean(ipos,1);

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
    fprintf('\nREAL ROOM bis:%2.2f   coh:%2.2f', room_sfep_bits, room_sfep_coherence)
    fprintf('\nREAL ARENA bis:%2.2f   coh:%2.2f', arena_sfep_bits, arena_sfep_coherence)
end


%% Generate null ipos driven only by position
if ~isfield(ms.room, 'inferred_ipos') && ~isfield(ms.arena, 'inferred_ipos')
    if verbose==true
        fprintf('\n\t! Generating inferred ipos for comparison...')
    end
    spks = ms.spks;
    % pfr=ms.room.pfields_smooth;
    % pfa=ms.arena.pfields_smooth;
    pfr=ms.room.pfields;
    pfa=ms.arena.pfields;
    
    spksroom  = NaN(size(spks));
    spksarena = NaN(size(spks));
    
    [~, ~, ~, xbr, ybr] = histcounts2(ms.room.x, ms.room.y, ms.params.pos_bins, ms.params.pos_bins);
    [~, ~, ~, xba, yba] = histcounts2(ms.room.x, ms.room.y, ms.params.pos_bins, ms.params.pos_bins);
    xbr(xbr==0) = interp1(find(xbr~=0), xbr(xbr~=0), find(xbr==0), 'nearest');
    ybr(ybr==0) = interp1(find(ybr~=0), ybr(ybr~=0), find(ybr==0), 'nearest');
    xba(xba==0) = interp1(find(xba~=0), xba(xba~=0), find(xba==0), 'nearest');
    yba(yba==0) = interp1(find(yba~=0), yba(yba~=0), find(yba==0), 'nearest');
    
    for i = 1:length(xbr)
        spksroom(:,i)  = pfr(:, ybr(i), xbr(i));
        spksarena(:,i) = pfa(:, yba(i), xba(i));
    end
    spksroom = normalize_rows(spksroom);
    spksarena = normalize_rows(spksarena);
    
    arenams = ms; arenams.spks = spksarena;
    roomms = ms; roomms.spks = spksroom;
    arenams.params.skip_ensemble = false;
    roomms.params.skip_ensemble = false;
    [~, ~, ~, arena_inferred_ipos] = Fenton_ipos(arenams, .25, 'arena', arenams.params);
    [~, ~, ~, room_inferred_ipos] = Fenton_ipos(roomms, .25, 'room', roomms.params);
    if verbose==true
        fprintf(' done!\n')
    end
    ms.arena.inferred_ipos = arena_inferred_ipos;
    ms.room.inferred_ipos = room_inferred_ipos;
else
    if verbose==true
        fprintf('\n\tUsing previous inferred ipos\n')
        arena_inferred_ipos = ms.arena.inferred_ipos;
        room_inferred_ipos  = ms.room.inferred_ipos;
    end
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
    fprintf('\nINFER ROOM bis:%2.2f   coh:%2.2f', room_sfep_inf_bits, room_sfep_inf_coherence)
    fprintf('\nINFER ARENA bis:%2.2f   coh:%2.2f', arena_sfep_inf_bits, arena_sfep_inf_coherence)
end

if plotting==true
    figure;
    subplot_tight(2,3,1); imagesc(room_sfep_map_smooth); axis image off; title('room sfep')
    subplot_tight(2,3,2); imagesc(room_sfep_map_inf_smooth); axis image off; title('room control')
    subplot_tight(2,3,3); imagesc(room_sfep_map_smooth-room_sfep_map_inf_smooth); axis image off; title('room sfep diff')
    subplot_tight(2,3,4); imagesc(arena_sfep_map_smooth); axis image off; title('arena sfep')
    subplot_tight(2,3,5); imagesc(arena_sfep_map_inf_smooth); axis image off; title('arena control')
    subplot_tight(2,3,6); imagesc(arena_sfep_map_smooth-arena_sfep_map_inf_smooth); axis image off; title('arena sfep diff')
    figure;
    subplot_tight(2,3,1); imagesc(room_sfep_map); axis image off; title('room sfep')
    subplot_tight(2,3,2); imagesc(room_sfep_map_inf); axis image off; title('room control')
    subplot_tight(2,3,3); imagesc(room_sfep_map-room_sfep_map_inf); axis image off; title('room sfep diff')
    subplot_tight(2,3,4); imagesc(arena_sfep_map); axis image off; title('arena sfep')
    subplot_tight(2,3,5); imagesc(arena_sfep_map_inf); axis image off; title('arena control')
    subplot_tight(2,3,6); imagesc(arena_sfep_map-arena_sfep_map_inf); axis image off; title('arena sfep diff')
end

room_diff = nanmean(nanmean(room_sfep_map  - room_sfep_map_inf));
arena_diff = nanmean(nanmean(arena_sfep_map - arena_sfep_map_inf));
room_absdiff = nanmean(nanmean(abs(room_sfep_map  - room_sfep_map_inf)));
arena_absdiff = nanmean(nanmean(abs(arena_sfep_map - arena_sfep_map_inf)));

if verbose==true
    fprintf('\n\tROOM  DIFF: %2.2f', room_diff)
    fprintf('\n\tARENA DIFF: %2.2f', arena_diff)
    fprintf('\n\n')
end


sfep = [];
% maps
sfep.room_map               = room_sfep_map;
sfep.room_map_smooth        = room_sfep_map_smooth;
sfep.arena_map              = arena_sfep_map;
sfep.arena_map_smooth       = arena_sfep_map_smooth;
sfep.infer.room_map         = room_sfep_map_inf;
sfep.infer.room_map_smooth  = room_sfep_map_inf_smooth;
sfep.infer.arena_map        = arena_sfep_map_inf;
sfep.infer.arena_map_smooth = arena_sfep_map_inf_smooth;

% metrics
sfep.room_diff              = room_diff;
sfep.arena_diff             = arena_diff;
sfep.room_absdiff           = room_absdiff;
sfep.arena_absdiff          = arena_absdiff;

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
    [sfep_map, ~] = make_occupancymap_2D(x, y, sfepsig, xbins, ybins);
    sfep_map = sfep_map./vmap;
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
