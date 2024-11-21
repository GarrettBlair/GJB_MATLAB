%% USED PARAMS FOR MOUSE PKMZ APA IMAGING EXPTS
this_params_ver = 'v2.1';

params = [];
params.arena_radius             = 20; % in cm
params.arena_center             = [ 127.5, 127.5]; % pixel center of behav cam, [x,y]
params.pixpercm                 = 6.25; % pixel radius of behav cam env
params.behav_fps                = 30;
params.behav_smoothing_interval = .25; % in seconds, length of smoothing kernel


params.pos_bins                 = [-22 -17.5:2.5:17.5 22]; % in cm, x and y
params.yaw_bins                 = -pi:pi/8:pi;
params.rho_bins                 = [0 7 14 21];
params.occupancy_thresh         = 1; % in seconds, minimum time in each bin to be counted for place map
params.pfield_kernel_radius     = 3; % kernel ends up being [n*2 + 1] in bins
params.smin_vals                = -50:5:-10; % smin values used to create the 'deconv_sweep.mat'
% params.bayesian_random_shuffle  = 500;
params.num_random_shuffle_pcell = 100; 
params.num_random_shuffle_decode= 10; 
% params.speed_thresh             = 5; % speed thresh in cm/sec
params.num_partitions           = 2;
% params.max_spd_thresh           = 100;
% params.min_spd_thresh           = 5;
params.max_spd_thresh           = 100;
params.min_spd_thresh           = -1; % if goes below this, counts as not moving
params.min_samples              = 10;
params.ipos_int_time            = .25; % seconds, binning time for computing momentary spatial information
params.x_fold_training          = 5; % decoding split train/test partition number
params.decode_bins_theta        = -pi:pi/12:pi; % polar angle binning
params.decode_bins_rho          = linspace(0,22,4); % anglar distance binning
% params.rotate_behav             = true;
params.nan_interp               = true;
params.remove_all_bad_caiman_segs = false;
params.remove_segs_outside      = true;
params.correct_dt               = true; % correct for large jumps in timestamp file when constructing vmap
params.plotting                 = false;
params.skip_contour_bounding    = false;
params.cameraName               = {'HPC_miniscope1'};%, 'ACC_miniscope2'}; % will default to 'MiniLFOV'
% params.reuse_contour_crop       = 'Crop_params.mat'; % use the previous ms file contour crop, unless empty
params.reuse_contour_crop       = 'bounding_box.mat'; % use the previous ms file contour crop, unless empty

params.parfor_progbar           = false; % whether to show the parfor progress figure
params.polar_decode             = false; % whether to use polar coords for svm decoding
params.skip_ensemble            = false; % flag to skip ensemble ipos calc






