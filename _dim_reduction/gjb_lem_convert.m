function ReducedDataArray=DimentionalityReduction_Ver1(data,p_neighbors_vec)
%Inputs:

%activity_mat_active_frames: N-by-T activity matrix, N is number of neurons
%and T is number of time frames (includes only active frames)

%p_neighbors_vec: 1-by-2 vector of portion of neighbors that are considered close. First
%entry for the first iteraion and second for second iteration.

%Output:
%ReducedDataArray: Array of the data. First Cell original data, second cell data after first 
%iteration of dimentionality reduction, and third cell after second iteration of dimentionality reduction. 
%%

do_isomap = false;

ndims = 10;
p_neighbors1=p_neighbors_vec(1);
p_neighbors2=p_neighbors_vec(2);

% data = activity_mat_active_frames';
[N, nc] = size(data);

knn    = ceil(p_neighbors1*N); 
% knn
sigma2 = 25;

m = size(data,1);
dt = squareform(pdist(data));
[srtdDt,srtdIdx] = sort(dt,'ascend');
dt = srtdDt(1:knn+1,:);
nidx = srtdIdx(1:knn+1,:);

tempW  = ones(size(dt));

% build weight matrix
i = repmat(1:m,knn+1,1);
W = sparse(i(:),double(nidx(:)),tempW(:),m,m); 
W = max(W,W'); % for undirected graph.

% The original normalized graph Laplacian, non-corrected for density
ld = diag(sum(W,2).^(-1/2));
DO = ld*W*ld;
DO = max(DO,DO');

% get eigenvectors
clearvars -except DO knn p_neighbors2 v nc ndims data p_neighbors_vec do_isomap
% [v,d] = eigs(DO,10,'la'); % GJB get all eigs
% v_varexplained = cumsum(diag(d))./sum(diag(d));
% ndim_half = find(v_varexplained>=.5,1);
if do_isomap
    fprintf('ISOMAP\n')
    options.dims = 1:ndims;
    options.verbose = true;
    options.display = true;
    [Y, R, E] = IsomapII(DO, 'k', knn, options); 
    v = Y.coords{ndims}';
else
    [v,d] = eigs(DO,ndims,'la');
end

data2 = [v(:,1:ndims)];
% data2 = [v(:,1:ndim_half)];
N                = size(data2,1);
knn    = ceil(p_neighbors2*N); % each patch will only look at its knn nearest neighbors in R^d
% knn
m                = size(data2,1);
dt               = squareform(pdist(data2));
[srtdDt,srtdIdx] = sort(dt,'ascend');
dt               = srtdDt(1:knn+1,:);
nidx             = srtdIdx(1:knn+1,:);

tempW  = ones(size(dt));

i = repmat(1:m,knn+1,1);
W = sparse(i(:),double(nidx(:)),tempW(:),m,m); 
W = max(W,W');

ld = diag(sum(W,2).^(-1/2));
DO = ld*W*ld;
DO = max(DO,DO');
%%%%%%%%%%%%%% GJB changed to 4 fom 10 dims
clearvars -except DO p_neighbors2 v nc ndims do_isomap

if do_isomap
    fprintf('ISOMAP\n')
    options.dims = 1:ndims;
    options.verbose = true;
    options.display = true;
    [Y, R, E] = IsomapII(DO, 'k', knn, options); 
    v2 = Y.coords{ndims}';
else
    [v2,dd] = eigs(DO,ndims,'la');
end
% [v2,~] = eigs(DO,4,'la');
v2 = v2(:, 1:ndims);
% v2_varexplained = cumsum(diag(dd))./sum(diag(dd));
ReducedDataArray{1}=[]; % data;
ReducedDataArray{2}=v;
ReducedDataArray{3}=v2;

% varianceExplained{1}=[]; % data;
% varianceExplained{2}=v_varexplained;
% varianceExplained{3}=v2_varexplained;