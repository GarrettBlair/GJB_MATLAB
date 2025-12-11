function [arr_z] = zscore_nan(arr, dim)
%ZSCORE_NAN_DIM Z-score normalization along any dimension, ignoring NaNs.
%
%   Z = ZSCORE_NAN_DIM(X, DIM) standardizes X by subtracting the mean and 
%   dividing by the standard deviation along dimension DIM, while omitting NaNs.
%
%   Inputs:
%       X   - Input N-dimensional array
%       DIM - Dimension along which to perform z-score normalization
%
%   Output:
%       Z   - Z-score normalized array (same size as X)
%
%   Example:
%       X = randn(5, 4, 3); X(2,:,1) = NaN;
%       Z = zscore_nan_dim(X, 1);

    % Compute mean and std along specified dimension, omitting NaNs
    mu = mean(arr, dim, 'omitnan');
    sigma = std(arr, 0, dim, 'omitnan');

    % Expand size to match input for broadcasting
    sz = ones(1, ndims(arr));
    sz(dim) = size(arr, dim);

    mu = repmat(mu, sz);
    sigma = repmat(sigma, sz);

    % Perform z-score normalization
    arr_z = (arr - mu) ./ sigma;

    % Avoid division by zero: set std==0 to NaN in result
    arr_z(sigma == 0) = NaN;
end
