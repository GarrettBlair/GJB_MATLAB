function [decodeStruct_room, decodeStruct_arena] = APA_within_sess_decoding(ms, decode_params)
% pos_bins = [-45, -36:4:36, 45];
% occ_thresh = 1; % minimum occupancy in time (sec) to be valid fo mapping error
% x_fold_training = 10;
timeres = decode_params.ipos_int_time;

[room_struct] = bin_session_data(ms, timeres, 'room', decode_params);
[arena_struct] = bin_session_data(ms, timeres, 'arena', decode_params);

is_moving = speed_calc_gb(arena_struct.x, arena_struct.y) > decode_params.min_spd_thresh;
fprintf('\n\tDecoding room position ')
[room_struct, decodeStruct_room] = decode_sub(room_struct, decode_params, is_moving);

fprintf('\n\tDecoding arena position ')
[arena_struct, decodeStruct_arena] = decode_sub(arena_struct, decode_params, is_moving);

end

function [frame_struct, decode_struct] = decode_sub(frame_struct, decode_params, is_moving)
pos_bins = decode_params.pos_bins;%[-45, -36:4:36, 45];
occ_thresh = decode_params.occupancy_thresh; % minimum occupancy in time (sec) to be valid fo mapping error
x_fold_training = decode_params.x_fold_training;
ispolar = decode_params.polar_decode;

spks_bin = frame_struct.spks_bin;
x = frame_struct.x;
y = frame_struct.y;
x(isnan(x)) = interp1(find(~isnan(x)), x(~isnan(x)), find(isnan(x)), 'nearest', 'extrap');
y(isnan(y)) = interp1(find(~isnan(y)), y(~isnan(y)), find(isnan(y)), 'nearest', 'extrap');
occ_map = frame_struct.vmap;

if ispolar==true
    [varDim1, varDim2] = cart2pol(x,y);
    frame_struct.bins1          = decode_params.decode_bins_theta;
    frame_struct.bins2          = decode_params.decode_bins_rho;
    frame_struct.theta          = varDim1;
    frame_struct.r              = varDim2;
else % cartesian
    varDim1 = x;
    varDim2 = y;
    frame_struct.bins1          = pos_bins;
    frame_struct.bins2          = pos_bins;
end

frame_struct.bins1_center           = frame_struct.bins1(2:end) - abs(diff(frame_struct.bins1))/2;
frame_struct.bins2_center           = frame_struct.bins2(2:end) - abs(diff(frame_struct.bins2))/2;
frame_struct.train_flag             = is_moving;
frame_struct.num_random_shuffle     = decode_params.num_random_shuffle_decode;
frame_struct.parfor_progbar         = decode_params.parfor_progbar;

decode_struct = svm_decode_sub(frame_struct, {varDim1, varDim2}, {frame_struct.bins1 frame_struct.bins2},...
    spks_bin, x_fold_training, 'eucl_dist_real');
decode_struct.spks = spks_bin;

[sum_all, counts_all]       = make_occupancymap_2D(x,  y,  decode_struct.pred_err, pos_bins, pos_bins);
av_err_map = sum_all./counts_all;
av_err_map(occ_map <= occ_thresh) = NaN;
decode_struct.error_map     = av_err_map;

end
