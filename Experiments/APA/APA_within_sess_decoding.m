function [decodeStruct_room, decodeStruct_arena] = APA_within_sess_decoding(ms, decode_params)
% pos_bins = [-45, -36:4:36, 45];
% occ_thresh = 1; % minimum occupancy in time (sec) to be valid fo mapping error
% x_fold_training = 10;
% timeres = .25;

pos_bins = decode_params.pos_bins;%[-45, -36:4:36, 45];
occ_thresh = decode_params.occupancy_thresh; % minimum occupancy in time (sec) to be valid fo mapping error
x_fold_training = decode_params.x_fold_training;
timeres = decode_params.ipos_int_time;

% t = ms.timestamps./1000;
% spks = ms.spks; % normalize_rows(ms.neuron.S_matw);
% spks(isnan(spks)) = 0;
% spks_bin = bin_spks_time(spks, timeres, t, false);
% spks_bin(isnan(spks_bin))=0;

tempparams.pos_bins = pos_bins; % in cm, x and y
tempparams.occupancy_thresh = occ_thresh; % minimum in seconds
tempparams.skip_ensemble = true; %can skip ensemble prob calc
%%
% % % spd = ms.speed_epochs; % [spd; spd(end)];
% % % is_moving = average_spks_time(spd', timeres, t, false, 'mean')>=.5;
[room_struct] = bin_session_data(ms, timeres, 'room', decode_params);
is_moving = speed_calc_gb(room_struct.x, room_struct.y) > decode_params.min_spd_thresh;


spks_bin = room_struct.spks_bin; 

ref_struct = room_struct;

x = ref_struct.x; 
y = ref_struct.y; 
x(isnan(x)) = interp1(find(~isnan(x)), x(~isnan(x)), find(isnan(x)), 'nearest', 'extrap');
y(isnan(y)) = interp1(find(~isnan(y)), y(~isnan(y)), find(isnan(y)), 'nearest', 'extrap');
occ_map = ref_struct.vmap; 
[pos_angle, rth] = cart2pol(x,y);

room_struct.abins          = decode_params.decode_bins_theta;%-pi:pi/12:pi; % range for binning angular position
room_struct.rbins          = decode_params.decode_bins_rho;%0:15:45; % range fo binning polar distance
room_struct.abin_center    = room_struct.abins(2:end) - abs(diff(room_struct.abins))/2;
room_struct.rbin_center    = room_struct.rbins(2:end) - abs(diff(room_struct.rbins))/2;
room_struct.theta          = pos_angle;
room_struct.r              = rth;
room_struct.train_flag     = is_moving;
room_struct.num_random_shuffle  = decode_params.num_random_shuffle_decode;
room_struct.parfor_progbar = decode_params.parfor_progbar;

fprintf('\n\tDecoding room position  ')
decodeStruct_room = svm_decode_sub(room_struct, {pos_angle, rth}, {room_struct.abins room_struct.rbins}, spks_bin, x_fold_training, 'eucl_dist_real');
decodeStruct_room.spks = spks_bin;

[sum_all, counts_all]       = make_occupancymap_2D(x,  y,  decodeStruct_room.pred_err, pos_bins, pos_bins);
av_err_map = sum_all./counts_all; 
av_err_map(occ_map <= occ_thresh) = NaN;
decodeStruct_room.error_map     = av_err_map;

% [sum_all, ~]       = make_occupancymap_2D(x,  y,  decodeStruct_room.rand_err, pos_bins, pos_bins);
% av_err_map = sum_all./counts_all; 
% av_err_map(occ_map <= occ_thresh) = NaN;
% decodeStruct_room.rand_error_map     = av_err_map;

%%
[arena_struct] = bin_session_data(ms, timeres, 'arena', decode_params);
% [~, arena_struct, ~] = Fenton_ipos(ms, timeres, 'arena', tempparams);
ref_struct = arena_struct;

x = ref_struct.x; 
y = ref_struct.y; 
occ_map = ref_struct.vmap; 
[pos_angle, rth] = cart2pol(x,y);

arena_struct.abins          = decode_params.decode_bins_theta;%-pi:pi/12:pi; % range for binning angular position
arena_struct.rbins          = decode_params.decode_bins_rho;%0:15:45; % range fo binning polar distance
arena_struct.abin_center    = arena_struct.abins(2:end) - abs(diff(arena_struct.abins))/2;
arena_struct.rbin_center    = arena_struct.rbins(2:end) - abs(diff(arena_struct.rbins))/2;
arena_struct.theta          = pos_angle;
arena_struct.r              = rth;
arena_struct.train_flag     = is_moving;
arena_struct.num_random_shuffle  = decode_params.num_random_shuffle_decode;
arena_struct.parfor_progbar = decode_params.parfor_progbar;

fprintf('\n\tDecoding arena position ')
decodeStruct_arena = svm_decode_sub(arena_struct, {pos_angle, rth}, {arena_struct.abins arena_struct.rbins}, spks_bin, x_fold_training, 'eucl_dist_real');
decodeStruct_arena.spks = spks_bin;

[sum_all, counts_all]       = make_occupancymap_2D(x,  y,  decodeStruct_arena.pred_err, pos_bins, pos_bins);
av_err_map = sum_all./counts_all; 
av_err_map(occ_map <= occ_thresh) = NaN;
decodeStruct_arena.error_map     = av_err_map;

