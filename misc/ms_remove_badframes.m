function [ms] = ms_remove_badframes(ms, badframes, verbose)

if ~exist('verbose', 'var') && ~isempty(verbose)
    verbose = false;
end
try
frameVec = ms.frameNum;
catch
frameVec = 1:size(ms.S,2);
end
if ~isempty(badframes)
    if islogical(badframes)
        good_idx = find(~badframes);
    elseif isscalar(badframes)
        good_idx = setdiff(frameVec, badframes);
    end
    frames = length(frameVec);
    fields = fieldnames(ms);
    for i = 1:size(fields,1)
        temp = eval(sprintf('ms.%s', fields{i}));
        size1 = size(temp);
        if size(temp, 2) == frames
            can_remove = true;
            temp = temp(:, good_idx);
        elseif size(temp, 1) == frames
            can_remove = true;
            temp = temp(good_idx, :);
        elseif all(size(temp) == [1, frames])
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
            eval(sprintf('ms.%s = temp;', fields{i}));
        else
            if verbose
                fprintf('Kept\t- %s:  \n\tsize: %d x %d\n', fields{i}, size1(1), size1(2))
            end
        end
    end
end
ms.bad_frames_removed = true;
end