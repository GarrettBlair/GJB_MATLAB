load swiss_roll_data

D = L2_distance(X_data(:,1:10000), X_data(:,1:10000), 1);
options.dims = 1:10;
% [Y, R, E] = Isomap(D, 'k', 7, options);
% 
% 
% D = L2_distance(X_data(:,1:10000), X_data(:,1:10000), 1);
% options.dims = 1:10;
% options.landmarks = 1:50;
% [Y, R, E] = IsomapII(D, 'k', 7, options);

clear
%% IsoMap example %%
%%%%%%%%%%%%%%%%%%%%

% ddir = 'C:\Users\gjb326\Desktop\Sample Data\MiniData Hipp Blair\';
% fname = sprintf('%sHipp%d\\caiman_cnmfe_out.mat', ddir, 8);
fname = 'C:\Users\gjb326\Desktop\Sample Data\MiniData Hipp Blair\Hipp8\Hipp8_linear4.mat';
load(fname);
%%
% nspks = double(normalize_matrix(S(:,1:2048)));
normspks = double(normalize_matrix(ms.neuron.S));

% nspks = double(normalize_matrix(ms.neuron.C_raw));
% s = ms.neuron.ispks; s(s<0) = 0;
% normspks = double(normalize_matrix(s));
[nsegs, nsamples] = size(normspks);
% kern = [zeros(1,4), ones(1,5)]; kern = kern./sum(kern(:));
kern = [1]; kern = kern./sum(kern(:));
smooth_spks = conv2(1, kern, normspks, 'same');
smooth_spks = normalize_matrix(smooth_spks);
% smooth_spks = smooth_spks(1:2:end,1:2:8000);
[pco_tau, pco_prob] = Fenton_pco(smooth_spks, 16, true);
figure(46); clf; imagesc(pco_tau); colorbar

options.dims = 1:3;%1:60; %1:3 %range of dimensions to calculate typ 2:20 for dim calculations, 3 for plot
options.landmarks = 1:size(smooth_spks,2); %recommended value 1:size(datalength);
options.verbose = true;
options.overlay = true;
options.display = true;
nNeighbours = 3;
D = L2_distance(smooth_spks, smooth_spks, 1);
[Y, R, E] = IsomapII(D, 'k', nNeighbours, options); 
v1 = Y.coords{3,1};

%% Ziv Laplacian example  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fakepos = ones(size(smooth_spks,2), 1);
% p_neighbors_vec=[0.075/15 0.075];
p_neighbors_vec=[0.075/15 0.075];

% ActivityThreshold=1;
% downsample_flag=2;
% ActvitySmoothingSize=2;
% NumOfBins=24;
% [all_trials_activity_smoothed_mat position_active_frames activity_mat_active_frames legitimacy_vec]=...
%     SmoothingAndThresholdingData(position_per_frame,neuronal_activity_mat,downsample_flag,ActvitySmoothingSize,ActivityThreshold);
% [all_trials_activity_smoothed_mat position_active_frames activity_mat_active_frames legitimacy_vec]=...
%     SmoothingAndThresholdingData(fakepos, smooth_nspks, 0, 0, 0);
rand('seed',0)
ReducedDataArray=DimentionalityReduction_Ver1(1*(smooth_spks>0), p_neighbors_vec);

v2=ReducedDataArray{3}';

%% MIND example %%
%%%%%%%%%%%%%%%%%%
ms = [];
ms.t(:,1) = (1:size(smooth_spks,2))/11;
ms.x(:,1) = ones(size(smooth_spks,2), 1);
ms.f = smooth_spks';

parameters = [];
parameters.normalizepTF = true;
parameters.dt = 1;
parameters.pca.n = .95;   % how many principal components to retain (or fraction of variance)
parameters.lm.n = 5000;    % how many landmark points
parameters.prune_lm_by_time = false; % prune landmarks so they are separated in time, if requested
% initialize ppca forest
parameters.dim_criterion = 10; % guess
parameters.ndir = 3; % guess
parameters.min_leaf_pts = 10; % guess
parameters.ntrees = 10; % guess 10?
parameters.verbose = true; % logical
%
parameters.rwd.type = 'discrete';    % type of random walk ('continuous' or 'discrete')
parameters.rwd.sym = 'avg';    % how to symmetrize global distances ('avg' or 'min')
parameters.rwd.all_geo = true;    % if true, all distances will be geodesic distances. otherwise,
    % distances between connected points will be ideal local distances and
    % distances between non-connected points will be filled in with
    % geodesic distances.
parameters.rwd.d = 3;% dimensionality of space in which random walk is performed (only used for continuous random walk)
parameters.rwd.var_scale = .05;
    % variance of diffusion kernel, expressed as fraction of maximum
    % possible variance (only used for continuous random walk; shouldn't
    % matter much)

dat = mindAsFunction(ms,parameters);
% saves output in struct 'dat' with fields/sub-structs:
%   data: original data
%   pca: global pca
%   forest: ppca forest for estimating transition probabilities
%   lm: landmark points
%   tp: transition probabilities
%   rwd: random walk distances
    
parameters.embed.d = 2; % guess?
parameters.embed.mode = 'mds'; % 'mds' or 'cmds' or 'hierarchical'; guess?
parameters.learnmapping = true; % guess?
parameters.mapping.mode = 'gp'; % 'lle' or 'wn' or 'tlle' or 'gp'; guess?
% parameters.mapping.k = 5; % LLE mode; choices of k for k nearest neighbors; guess?
% parameters.mapping.lambda = .5; % LLE mode; % choices of regularization parameter; guess?
% parameters.mapping.nfolds_lle = 5; % LLE mode; idk guess?
% parameters.mapping.h; % WN w-n kernel mode; choices of h
% parameters.tmargin_sec = 10; % guess?
% % TLLE uses same params as LLE but is "lle with time-blocked cv"
[dat, embed] = embedAsFunction(dat, parameters);

v3 = dat.embed.f2m.y(:,1:2)';
v3 = cat(1, v3, v3(2,:));




%%
h = figure(22); clf
nsamples = 10000;
downSample = 4;
trail = 20;
draw_history = false;

y_all = cell(3,1);
y_all{1} = v1;
y_all{2} = v2;
% y_all{3} = v3;

plot_dim_by_time(y_all, nsamples, downSample, trail, h, draw_history)



%%  




function plot_dim_by_time(y, nsamples, downSample, trailingSamples, figHandle, keep_history)
nb = trailingSamples;
num_data = length(y);
for i = nb+1:downSample:nsamples
    %%
figure(figHandle); clf
    for nn = 1:num_data
        yy = y{nn};
        if ~isempty(yy)
%         figure(figHandle);
        subplot_tight(1, num_data, nn, [.03 .06])
        cla; hold on
        %     plot3(yy(1,:),yy(2,:),yy(3,:), '-', 'Color', [.7 .7 .7], 'LineWidth',2)
        scatter3(yy(1,:),yy(2,:),yy(3,:), 5, 'k.', 'MarkerEdgeColor', [.3 .3 .3])
        plot3(yy(1,i-nb:i),yy(2,i-nb:i),yy(3,i-nb:i), 'r.-', 'LineWidth',2)
        scatter3(yy(1,i-nb:i),yy(2,i-nb:i),yy(3,i-nb:i), 15, 'r.', 'MarkerEdgeColor', [.8 .8 .8])
        if keep_history
        plot3(yy(1,1:i-nb),yy(2,1:i-nb),yy(3,1:i-nb), '-', 'Color', [.9 .6 .6], 'LineWidth', .02)
        end
        set(gca, 'View', [90+i/5, 36])
        %     axis([-10 10 -10 10 -10 10])
        axis tight
        end
    end
        drawnow
end
end
