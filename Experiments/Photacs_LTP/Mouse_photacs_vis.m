clear 
recording_folder = 'C:\Users\gjb326\Desktop\Recording Data\GarrettBlair\MOUSE_PHOTACS\2022-04-21_17-36-00\Record Node 103';
json_file = '\experiment1\recording1\structure.oebin'; % ".oebin" is the JSON file they refer to in docs
% json_file = '\experiment1\recording2\structure.oebin'; % ".oebin" is the JSON file they refer to in docs

fname = strcat(recording_folder, json_file);

input_index = 1; % DAQ board input, shouldn't need to change

% load continuous recording data
A = load_open_ephys_binary(fname, 'continuous', input_index);

% all channels should have the same bitvolt value
bv = A.Header.channels(1).bit_volts;
% probe_order = [21 12 22 11 23 10 24 9 13 20 14 19 15 18 16 17];
probe_order = [9:9+16];
D = bv*A.Data; 
fs = double(A.Header.sample_rate);
timestamps = double(A.Timestamps);
timestamps = (timestamps-timestamps(1))/fs; % convert sample num to second
%% TTL
raw_ttl = D(33,:); % 33 is the first TTL index (1-32 are the data channels for this headstage)
D = D(probe_order,:);
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

ttl_duration = median(timestamps(ttl_offset) - timestamps(ttl_onset));
ttl_interval = median(timestamps(ttl_onset(2:end)) - timestamps(ttl_onset(1:end-1)));
ttl_std = std(timestamps(ttl_onset(2:end)) - timestamps(ttl_onset(1:end-1)));
ttl_freq = 1/ttl_interval;
% if ttl_std>1
%     warning('Inconsistent timing')
% end
% 
% % plot them to check
% figure; clf; hold on
% plot(timestamps, ttl); 
% % plot(delta_ttl*.5);
% % plot(sign_ttl*.2); 
% scatter(timestamps(ttl_onset), ones(num_ttl,1), 'go') % show onsets
% scatter(timestamps(ttl_offset), ones(num_ttl,1), 'ro') % show offsets
% legend('TTL', 'on', 'off')
% ylim([ -.5 1.5])
%% Plot data 
inds = 1:size(D,2); % plot all samples
figure; clf;
hold on

ch2plt = [1:16];%[2 4 5 7]; % if you only want a couple

if num_ttl>0
    for i = 1:num_ttl
        pvect = [timestamps(ttl_onset(i)), -2,  timestamps(ttl_offset(i))-timestamps(ttl_onset(i)), 5+length(ch2plt)];
        rectangle('Position', pvect, 'FaceColor', [.4 1 .5], 'EdgeColor', 'none')
    end
end

for i = 1:length(ch2plt)
    single_channel = D(ch2plt(i),(inds));
    single_channel = single_channel/5000;%250; % arbitrary scaling
    plot(timestamps(inds), single_channel + i);
end
single_channel = raw_ttl;
single_channel = single_channel/5000;%250; % arbitrary scaling
plot(timestamps(inds), single_channel + i+1);

% plot ttl around data if you want
axis tight

% p = make_peth_TTL(D(ch2plt(1),(inds)), ttl_onset, 3e2);
% figure; imagesc(p)
%%
inds = [1.749711865164762e+06, 1.963884185272361e+06];
% subD = D(:, inds(1):inds(2));
subD = D;
m = NaN(16);
for i = 1:16
    for j = 1:16
        m(i,j) = corr(subD(i,:)', subD(j,:)');
    end
end
ms = m(1,:);
[~, ord] = sort(ms, 'descend');

figure; imagesc([D(ord,:); D])