function [neuron] = remove_segments(neuron, bad_seg_idx, verbose)
if ~exist('verbose', 'var') && ~isempty(verbose)
    verbose = false;
end
if ~isempty(bad_seg_idx)
    good_idx = setdiff([1:size(neuron.A,2)], bad_seg_idx);
    
    fields = fieldnames(neuron);
    [segs, frames] = size(neuron.C);
    %     dims = size(neuron.C);
    for i = 1:size(fields,1)
        temp = eval(sprintf('neuron.%s', fields{i}));
        size1 = size(temp);
        if all(size(temp) == [segs, frames])
            can_remove = true;
            temp = temp(good_idx, :);
        elseif all(size(temp) == [neuron.dims(1)*neuron.dims(2), segs])
            can_remove = true;
            temp = temp(:, good_idx);
        elseif all(size(temp) == [1, segs])
            can_remove = true;
            temp = temp(1, good_idx);
        else
            can_remove = false;
        end
        size2 = size(temp);
        if can_remove
            if verbose
                fprintf('Changed\t- %s:  \n\twas: %d x %d  \n\tnow: %d x %d\n', fields{i}, size1(1), size1(2), size2(1), size2(2))
            end
            eval(sprintf('neuron.%s = temp;', fields{i}));
        else
            if verbose
                fprintf('Kept\t- %s:  \n\tsize: %d x %d\n', fields{i}, size1(1), size1(2))
            end
        end
    end
    neuron.idx_components = good_idx;
    if ~isrow(bad_seg_idx)
        bad_seg_idx = bad_seg_idx';
    end
    neuron.idx_components_bad = unique([neuron.idx_components_bad bad_seg_idx]);
    neuron.idx_components = good_idx;
end
%     good_idx = setdiff([1:size(neuron.A,2)], bad_seg_idx);
%     neuron.A         = neuron.A(:, good_idx);
%     if isfield(neuron, 'fullA')
%         neuron.fullA = neuron.fullA(:, good_idx);
%     end
%     neuron.C         = neuron.C(good_idx, :);
%     neuron.YrA       = neuron.YrA(good_idx, :);
%     neuron.C_raw     = neuron.C + neuron.YrA;
%     neuron.S         = neuron.S(good_idx, :);
%     neuron.sn        = neuron.sn(good_idx);
%     neuron.SNR_comp  = neuron.SNR_comp(good_idx);
%     neuron.r_values  = neuron.r_values(good_idx);
%     neuron.lam       = neuron.lam(good_idx);
%     neuron.idx_components = good_idx;
%     neuron.idx_components_bad = [neuron.idx_components_bad bad_seg_idx'];
% end
end