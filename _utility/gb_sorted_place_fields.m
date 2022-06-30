function [ordered_segs, sorted_spks, sorted_colored_spks, sorted_pfields, sorted_pos, good_pcells] = gb_sorted_place_fields(ms, idx, prob_thresh)

if isempty(idx)
    idx = 1:size(ms.neuron.A,2);
end
if ~exist('prob_thresh', 'var') || isempty(prob_thresh)
    prob_thresh = 1;
end
% max_firing = [(ms.SL.spike_prob(idx,:)) (ms.SR.spike_prob(idx,:))...
%     (ms.LL.spike_prob(idx,:)) (ms.LR.spike_prob(idx,:))];
max_firing = [(ms.SL.info.info_prob(idx,:)) (ms.SR.info.info_prob(idx,:))...
    (ms.LL.info.info_prob(idx,:)) (ms.LR.info.info_prob(idx,:))];


% max_firing_loc = [max(ms.SL.field_center{end}(idx,:), [], 2) max(ms.SR.field_center{end}(idx,:), [], 2)...
%     max(ms.LL.field_center{end}(idx,:), [], 2) max(ms.LR.field_center{end}(idx,:), [], 2)];
max_firing_loc = [max(ms.SL.field_center(idx,:), [], 2) max(ms.SR.field_center(idx,:), [], 2)...
    max(ms.LL.field_center(idx,:), [], 2) max(ms.LR.field_center(idx,:), [], 2)];
[min_prob, maxpos] = min(max_firing, [], 2);

good_pcells = min_prob>prob_thresh;
% maxpos(good_pcells) = 0;
strs = {'SL', 'SR', 'LL', 'LR'};
% strs = {'SL', 'SR'};

%%
ordered_segs = [];
if isfield(ms.neuron, 'C_raw')
    act = ms.neuron.C_raw;
else
    act = ms.neuron.YrA+ms.neuron.C;
end
act(act<0) = 0;
act = normalize_matrix(act(idx,:));
spks = 1*(ms.neuron.ispks(idx,:)>0);


sorted_colored_spks = zeros(size(act,1), size(act,2), 3);
sorted_pos = zeros(size(act,1), 1, 3);
prev_segs = 1;
cmap = jet(8);
cmap = cmap(2:7,:);
% cmap = 1-zeros(4,3);%[0 1 0;0 0 1;0 .7 0;0 0 .7];
pathcmap = [0 0 .9; .2 .9 .2; .2 .2 .6; .2 .6 .2];
sorted_pfields = [];

for i = 1:length(strs)
    segs = find(maxpos == i);% & min_prob<=.05);
    [~, ord] = sort(max_firing_loc(segs, i), 'descend');
    sorted_pos(prev_segs:prev_segs+length(segs)-1,:) =ones(length(segs), 1)*pathcmap(i,:);
    
    for jj = 1:3
        sorted_colored_spks(prev_segs:prev_segs+length(segs)-1, :, jj) =  spks(segs(ord), :)*cmap(i,jj);
    end
    prev_segs = length(segs) + prev_segs;
    [~, ord] = sort(max_firing_loc(segs, i), 'ascend');
    ordered_segs = cat(1, ordered_segs, segs(ord));
    
    
%     temp = eval(sprintf('ms.%s.smoothed_place_fields{end}(idx,:)', strs{i}));
% %     temp = eval(sprintf('ms.%s.place_fields{end}', strs{i}));
%     sorted_pfields = cat(2, sorted_pfields, temp);
%     temp = temp(segs, :);
%     resize_field= NaN(size(temp,1), bestfieldsize);
%     for j = 1:size(temp,1)
%         xx = interp1(1:size(temp,2), temp(j,:), linspace(1, size(temp,2), bestfieldsize));
%         resize_field(j,:) = conv(xx, ones(5,1)./5, 'same');
%     end
%     norm_fields = cat(1, norm_fields, resize_field);
end
% [~, max_resize_pos] = max(norm_fields, [], 2);
% [~, max_resize_ord] = sort(max_resize_pos);
% norm_fields = norm_fields(max_resize_ord, :);
sorted_spks = spks(ordered_segs,:);

% sorted_pfields = normalize_matrix(sorted_pfields, 1);
