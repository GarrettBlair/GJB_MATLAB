function [neuron] = caiman_downsample(neuron, frames2keep, verbose)
if ~exist('verbose', 'var') && ~isempty(verbose)
    verbose = false;
end
good_idx = frames2keep;%setdiff([1:size(neuron.A,2)], bad_seg_idx);
[numsegs, numframes] = size(neuron.C);
if length(good_idx)~=length(numframes)
    fields = fieldnames(neuron);
    %     dims = size(neuron.C);
    for i = 1:size(fields,1)
        temp = eval(sprintf('neuron.%s', fields{i}));
        size1 = size(temp);
        if all(size(temp) == [numsegs, numframes])
            can_remove = true;
            temp = temp(:, good_idx);
%         elseif all(size(temp) == [neuron.dims(1)*neuron.dims(2), segs])
%             can_remove = true;
%             temp = temp(:, good_idx);
        elseif all(size(temp) == [1, numframes])
            can_remove = true;
            temp = temp(1, good_idx);
        elseif all(size(temp) == [numframes, 1])
            can_remove = true;
            temp = temp(good_idx, 1);
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
    neuron.good_frames_kept = frames2keep;
elseif verbose == true
    fprintf('All frames kept, no ds required\n')
end
end