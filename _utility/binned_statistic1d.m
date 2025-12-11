function [h_map, binx] = binned_statistic1d(x, signal, binsx, method)
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
%       signal - Vector of signal values to be aggregated per bin
%       binsx  - Vector of bin edges along the x-axis
%       method - String specifying the statistic to compute:
%                'sum', 'mean', 'median', 'nansum', 'nanmean', 'nanmedian'
%
%   Output:
%       h_map  - 1D matrix of aggregated values per x bin,
%                size is (length(binsx)-1)
%                Bins with no data are filled with NaN.
%
%   Example:
%       x = randn(1000,1);
%       s = rand(1000,1);
%       binsx = -3:0.5:3;
%       H = binned_statistic1d(x, s, binsx, 'mean');
%       plot(binsx, H);
%
%   See also: histcounts2, accumarray
%   GJB 5/16/2025

% Bin the data into 2D indices
[~, ~, binx] = histcounts(x, binsx);
nbinsx = length(binsx)-1;
h_map = NaN(nbinsx,1); % Preallocate output

valid_method = {'sum', 'mean', 'median', 'nansum', 'nanmean', 'nanmedian', 'mode', 'max', 'nanstd'};
% Check for invalid method
if ~any(strcmp(method, valid_method))
    warning('Accepted methods are: %s', strjoin(valid_method, ', '));
    error('%s is not a valid method!', method)
end

% Filter out values that fall outside binning range
valid = binx > 0;
binx = binx(valid);
signal = signal(valid);
[nc,nr] = size(binx);
if nc==1
    binx = binx';
end
% Convert 2D subscripts to linear indices
lin_idx = binx; % sub2ind([nbinsy, nbinsx], biny, binx);

% Use accumarray to compute the statistic
method_func = str2func(method);
try
    stats = accumarray(binx, signal, [nbinsy*nbinsx, 1], method_func, NaN);
catch
    % accumarray may not support 'nanmean' or 'nanmedian' directly in old MATLAB versions
    % Fallback if needed
    if strcmp(method, 'nanmean')
        method_func = @(x) mean(x, 'omitnan');
    elseif strcmp(method, 'nanmedian')
        method_func = @(x) median(x, 'omitnan');
    end
    stats = accumarray(binx, signal, [nbinsx, 1], method_func, NaN);
end

% Recast to h_map
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