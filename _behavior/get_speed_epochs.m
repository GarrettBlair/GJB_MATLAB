function [partition_ID, speed_epochs] = get_speed_epochs(spd_vect, params)
%% Split session into n segments with segments of continuous behavior above speed threshold
% Evaluates the speed vector to find continuous samples that have at least
% n samples above max_thresh, then finds when those epochs fall below
% min_thresh (min stretching is optional)
% Inputs: 
% num_parts - number of equal sample partitions
% speed samples above max thresh (high speed)
num_parts       = params.num_partitions;
max_spd_thresh  = params.max_spd_thresh;
min_spd_thresh  = params.min_spd_thresh;
min_samples     = params.min_samples;

%%
good_speed = spd_vect >= min_spd_thresh; 

spd_change_start = find([diff([0; good_speed])]);
if spd_change_start > 1
spd_change_end = spd_change_start - 1;
else
spd_change_end = spd_change_start(2:end) - 1;
end
good_speed_start = spd_change_start(good_speed(spd_change_start));
good_speed_end = spd_change_end(good_speed(spd_change_end));
if length(good_speed_start) ~= length(good_speed_end)
    good_speed_start = good_speed_start(1:min([length(good_speed_start) length(good_speed_end)]));
    good_speed_end = good_speed_end(1:min([length(good_speed_start) length(good_speed_end)]));    
end
speed_length = good_speed_end - good_speed_start;

good_length = speed_length(speed_length>=min_samples);
good_length_start = good_speed_start(speed_length>=min_samples);
good_length_end = good_speed_end(speed_length>=min_samples);
% Refine indices to where the rat slows down below min_spd
if exist('min_spd_thresh', 'var')~=0 && ~isempty(max_spd_thresh)
    stop_spd = spd_vect <= max_spd_thresh; % speed samples below min thresh (high speed)
    ginit = 1;
    stretch_start = [];
    stretch_end  = [];
    str_inds = 1;
    for i = 1:length(good_length_start)
        if ginit <= good_length_start(i)
            g1 = good_length_start(i);
            g2 = good_length_end(i);
            s1 = ginit + find(stop_spd(ginit:g1), 1, 'last');
            if isempty(s1)
                s1 = g1;
            end
            s2 = g2 + find(stop_spd(g2:end), 1, 'first') - 2;
            
            stretch_start(str_inds, 1) = s1;
            if ~isempty(s2)
                stretch_end(str_inds, 1) = s2;
            else
                 stretch_end(str_inds, 1) = length(spd_vect);
            end
            ginit = stretch_end(str_inds, 1)+1;
            str_inds = str_inds+1;
        end
    end
    good_length_start = stretch_start;
    good_length_end = stretch_end;
    good_length = good_length_end-good_length_start;
end

parts = floor(sum(good_length)/num_parts);
length_cumsum = cumsum(good_length);

part_start(1) = 1;
part_end(1) = find(length_cumsum>=parts, 1, 'first');
for i = 2:num_parts
    part_start(i) = part_end(i-1)+1;
    part_end(i) = find(length_cumsum>=parts*i, 1, 'first');    
end

% randomly redistribute the labels throughout the session
% [~, rand_ord] = sort(rand(length(good_length_start), 1));
% good_length_start = good_length_start(rand_ord);
% good_length_end   = good_length_end(rand_ord);
% no need garrett sigh


partition_ID = zeros(length(spd_vect), 1);
num_samples  = zeros(num_parts, 1);
for i = 1:num_parts
    t1 = good_length_start(part_start(i):part_end(i));
    t2 = good_length_end(part_start(i):part_end(i));
    for j = 1:length(t1)
        partition_ID(t1(j):t2(j)) = i;
        num_samples(i) = num_samples(i) + length(t1(j):t2(j));
    end
end
speed_epochs = partition_ID>0;


end
