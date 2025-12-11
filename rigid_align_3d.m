function [A_transformed, R, t] = rigid_align_3d(A, B)
% RIGID_ALIGN_3D Aligns two 3D matrices using only rotation and translation.
%
%   [A_aligned, R, t] = RIGID_ALIGN_3D(A, B)
%
%   Inputs:
%       A - M x N x 3 matrix of 3D points to align (source)
%       B - M x N x 3 matrix of 3D reference points (target)
%
%   Outputs:
%       A_aligned - Aligned version of A with same size as A
%       R         - 3x3 rotation matrix
%       t         - 1x3 translation vector
%
%   Requires A and B to be the same size and contain corresponding points.

    % Check input size
    if ~isequal(size(A), size(B))
        error('A and B must be the same size.');
    end

    % Reshape to M*N x 3
    A_pts = reshape(A, [], 3);
    B_pts = reshape(B, [], 3);

    % Remove NaNs (if any)
    valid = all(~isnan(A_pts) & ~isnan(B_pts), 2);
    A_pts = A_pts(valid, :);
    B_pts = B_pts(valid, :);

    % Subtract centroids
    centroid_A = mean(A_pts, 1);
    centroid_B = mean(B_pts, 1);
    A_centered = A_pts - centroid_A;
    B_centered = B_pts - centroid_B;

    % Compute optimal rotation using SVD (Kabsch algorithm)
    H = A_centered' * B_centered;
    [U, ~, V] = svd(H);
    R = V * U';

    % Correct for reflection
    if det(R) < 0
        V(:,3) = -V(:,3);
        R = V * U';
    end

    % Compute translation
    t = centroid_B - centroid_A * R;

    % Apply transformation
    A_transformed = A_pts * R + t;

%     % Reshape back to original shape
%     A_aligned = reshape(NaN(size(A,1)*size(A,2), 3), size(A,1), size(A,2), 3);
%     A_aligned(repmat(reshape(valid, size(A,1), size(A,2)), 1, 1, 3)) = reshape(A_transformed', [], 1);
end