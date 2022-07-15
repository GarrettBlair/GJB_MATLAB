%%
params = [];
params.arena_radius             = 40; % in cm
params.pos_bins                 = -40:4:40; % in cm, x and y
params.yaw_bin                  = -pi:pi/8:pi;
params.behav_smoothing_interval = .5; % in seconds, length of smoothing kernel
params.occupancy_thresh         = 1; % in seconds, minimum time in each bin to be counted for place map
params.pfield_kernel_radius     = 3; % kernel ends up being [n*2 + 1] in bins
% params.speed_thresh             = 5; % speed thresh in cm/sec
params.num_partitions           = 2;
params.max_spd_thresh           = 100;
params.min_spd_thresh           = 5;
params.min_samples              = 10;

params.rotate_behav             = true;
params.nan_interp               = true;
params.correct_dt               = true; % correct for large jumps in timestamp file when constructing vmap
params.plotting                 = false;

recdir = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\2022_06_28\16_13_39\';
load([recdir, 'MiniLFOV\deconv_sweep.mat'])
load("C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\processed_files\2022_06_28___16_13_39_ms_placecells_data.mat")
% [ms, behav, params] = APA_generate_data_struct(recdir, params);
% [ms] = APA_generate_placmaps(ms, params);

%%
smin = -50:5:50;
smin = -50:5:0;

smat = [];
smat_weighted = [];
for i = 1:length(smin)
    if smin(i)==0
    varName = sprintf('smin_None');
    smin_weight = 1;
    elseif smin(i)<0
    varName = sprintf('smin_neg%d', abs(smin(i)));
    smin_weight = abs(smin(i));
    else
    varName = sprintf('smin_%d', smin(i));
    smin_weight = abs(smin(i));
    end
    s = eval(varName);
    smat = cat(3, smat, s);
    smat_weighted = cat(3, smat_weighted, s*smin_weight);
    
end
S1 = squeeze(sum(smat,3));
S1 = S1(ms.neuron.idx_components,:);
S2 = squeeze(sum(smat_weighted,3));
S2 = S2(ms.neuron.idx_components,:);
figure(1); clf
stacked_traces(ms.neuron.C+ms.neuron.YrA, 1, {'k-'})
stacked_traces(S1, 1, {'r-', 'LineWidth', 1.5});
stacked_traces(S2, 1, {'b-', 'LineWidth', 1});
% S = normalize_rows(S);
% figure(9); clf; hold on
% for i = 1:size(S,1); t = 1.3*(S(i,:)+i); plot(t,'k'); end