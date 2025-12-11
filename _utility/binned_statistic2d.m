function [h_map, binx, biny] = binned_statistic2d(x, y, signal, binsx, binsy, method)
% BINNED_STATISTIC Compute 2D binned statistics of a signal.
%
%   h_map = BINNED_STATISTIC(x, y, signal, binsx, binsy, method)
%
%   This function bins the (x, y) coordinate data using the specified
%   bin edges (`binsx` and `binsy`) and computes the specified
%   statistic over the corresponding `signal` values in each bin.
%
%   Inputs:
%       x      - Vector of x-coordinates (same length as `signal`)
%       y      - Vector of y-coordinates (same length as `signal`)
%       signal - Vector of signal values to be aggregated per bin
%       binsx  - Vector of bin edges along the x-axis
%       binsy  - Vector of bin edges along the y-axis
%       method - String specifying the statistic to compute:
%                'sum', 'mean', 'median', 'nansum', 'nanmean', 'nanmedian'
%
%   Output:
%       h_map  - 2D matrix of aggregated values per (y, x) bin,
%                size is (length(binsy)-1, length(binsx)-1)
%                Bins with no data are filled with NaN.
%
%   Example:
%       x = randn(1000,1);
%       y = randn(1000,1);
%       s = rand(1000,1);
%       binsx = -3:0.5:3;
%       binsy = -3:0.5:3;
%       H = binned_statistic2d(x, y, s, binsx, binsy, 'mean');
%       imagesc(binsx, binsy, H); colorbar; axis xy;
%
%   See also: histcounts2, accumarray, binned_statistic1d
%   GJB 5/16/2025

% Bin the data into 2D indices
[~, ~, ~, binx, biny] = histcounts2(x, y, binsx, binsy);
nbinsx = length(binsx)-1;
nbinsy = length(binsy)-1;
h_map = NaN(nbinsy, nbinsx); % Preallocate output

valid_method = {'sum', 'mean', 'median', 'nansum', 'nanmean', 'nanmedian'};
% Check for invalid method
if ~any(strcmp(method, valid_method))
    warning('Accepted methods are: %s', strjoin(valid_method, ', '));
    error('%s is not a valid method!', method)
end

% Filter out values that fall outside binning range
valid = binx > 0 & biny > 0;
binx = binx(valid);
biny = biny(valid);
signal = signal(valid);

% Convert 2D subscripts to linear indices
lin_idx = sub2ind([nbinsy, nbinsx], biny, binx);

% Use accumarray to compute the statistic
method_func = str2func(method);
try
    stats = accumarray(lin_idx, signal, [nbinsy*nbinsx, 1], method_func, NaN);
catch
    % accumarray may not support 'nanmean' or 'nanmedian' directly in old MATLAB versions
    % Fallback if needed
    if strcmp(method, 'nanmean')
        method_func = @(x) mean(x, 'omitnan');
    elseif strcmp(method, 'nanmedian')
        method_func = @(x) median(x, 'omitnan');
    end
    stats = accumarray(lin_idx, signal, [nbinsy*nbinsx, 1], method_func, NaN);
end

% Reshape to 2D
h_map(:) = stats;

%% Older slower version
% function [h_map] = binned_statistic2d(x, y, signal, binsx, binsy, method)
% [h, ~, ~, binx, biny] = histcounts2(x, y, binsx, binsy);
% ux = unique(binx);
% uy = unique(biny);
% valid_method = {'sum', 'mean', 'median', 'nansum', 'nanmean', 'nanmedian'};
% h_map = NaN(size(h));
% 
% if ~any(strcmp(method, valid_method))
%     warning('Accepted methods are: %s', strjoin(valid_method, ', '));
%     error('%s is not a valid method!', method)
% else
%     method_func = str2func(method);
%     for iy = 1:length(uy)
%         for ix = 1:length(ux)
%             idx = binx == ux(ix) & biny == uy(iy);
%             if any(idx)
%                 h_map(uy(iy), ux(ix)) = method_func(signal(idx));
%             end
%         end
%     end
% end
% end


end