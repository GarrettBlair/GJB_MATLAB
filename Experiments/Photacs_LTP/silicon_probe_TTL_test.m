recording_folder = 'C:\Users\gjb326\Documents\Open Ephys\TTL_test\2022-04-19_14-20-28';
json_file = '\Record Node 108\experiment1\recording1\structure.oebin'; % ".oebin" is the JSON file they refer to in docs

fname = strcat(recording_folder, json_file);

input_index = 1; % DAQ board input, shouldn't need to change

% load continuous recording data
A = load_open_ephys_binary(fname, 'continuous', input_index);

% all channels should have the same bitvolt value
bv = A.Header.channels(1).bit_volts;

D = bv*A.Data; 
fs = double(A.Header.sample_rate);
timestamps = double(A.Timestamps);
timestamps = (timestamps-timestamps(1))/fs; % convert sample num to second
%% TTL
raw_ttl = D(33,:); % 33 is the first TTL index (1-32 are the data channels for this headstage)
% normalize ttl
ttl = raw_ttl - min(raw_ttl); ttl = ttl./max(ttl);
ttl_thresh = .5; % value that declares ttl as "on", between 0 and 1. Lower will decrease latency, higher increase
% binarize ttl
ttl = ttl>ttl_thresh;

% find the change indices
delta_ttl = [0 abs(diff(ttl))>0]; % zero padding accounts for diff changing the size
% sign of the change, onset (positive) vs offset (negative)
sign_ttl = [0 (diff(ttl))]; % zero padding accounts for diff changing the size

ttl_onset  = find(delta_ttl==1 & sign_ttl>0);
ttl_offset = find(delta_ttl==1 & sign_ttl<0);

% concatenate indices to include only onsets with a subsequent offset
if length(ttl_onset) < length(ttl_offset)
    valid_ttl = find(ttl_offset > ttl_onset(1));
    ttl_offset = ttl_offset(valid_ttl);
elseif length(ttl_onset) > length(ttl_offset)
    valid_ttl = length(ttl_offset);
    ttl_onset = ttl_onset(1:valid_ttl);
end
num_ttl = length(ttl_onset);

ttl_duration = mean(timestamps(ttl_offset) - timestamps(ttl_onset));
ttl_interval = mean(timestamps(ttl_onset(2:end)) - timestamps(ttl_onset(1:end-1)));
ttl_std = std(timestamps(ttl_onset(2:end)) - timestamps(ttl_onset(1:end-1)));
ttl_freq = 1/ttl_interval;
if ttl_std>1
    warning('Inconsistent timing')
end

% plot them to check
figure(1); clf; hold on
plot(timestamps, ttl); 
% plot(delta_ttl*.5);
% plot(sign_ttl*.2); 
scatter(timestamps(ttl_onset), ones(num_ttl,1), 'go') % show onsets
scatter(timestamps(ttl_offset), ones(num_ttl,1), 'ro') % show offsets
legend('TTL', 'on', 'off')
ylim([ -.5 1.5])
%% Plot data 
inds = 1:size(D,2); % plot all samples
figure(2); clf;
hold on
ch2plt = [1:33];%[2 4 5 7]; % if you only want a couple
for i = 1:length(ch2plt)
    single_channel = D(ch2plt(i),(inds));
    single_channel = single_channel/15000;%250; % arbitrary scaling
    plot(timestamps(inds), single_channel + i);    
end

% plot ttl around data if you want
if num_ttl>0
    for i = 1:num_ttl
       pvect = [timestamps(ttl_onset(i)), -2,  timestamps(ttl_offset(i))-timestamps(ttl_onset(i)), 5+length(ch2plt)];
       rectangle('Position', pvect, 'FaceColor', 'none', 'EdgeColor', 'm')
    end   
end
axis tight

p = make_peth_TTL(D(ch2plt(1),(inds)), ttl_onset, 3e2);
figure; imagesc(p)







