function [ms] = downsample_ms_struct(ms, ds)
% if ~exist('verbose', 'var') && ~isempty(verbose)
%     verbose = false;
% end
numFrames = length(ms.frameNum);
if isscalar(ds)
    ds_frames = 1:ds:numFrames;
elseif isvector(ds) && length(ds)==2
    ds_frames = ds(1):ds(2);
else
    ds_frames = ds;
end
ms = ds_struct(ms, ds_frames, numFrames);

end

function str = ds_struct(str, ds_frames, numFrames)
f = fieldnames(str);
nf = length(f);
for i = 1:nf
    temp = eval(sprintf('str.%s;', f{i}));
    if isstruct(temp)
        temp = ds_struct(temp, ds_frames, numFrames);
    else
        if size(temp,1)==numFrames
            temp = temp(ds_frames, :);
        elseif size(temp,2)==numFrames
            temp = temp(:, ds_frames);
        end
    end
    eval(sprintf('str.%s = temp;', f{i}));
    
end
end
