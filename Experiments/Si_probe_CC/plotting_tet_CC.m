% recording_folder = 'C:\Users\gjb326\Desktop\Sample Data\Fenton Si probes\RAT2_2022_02_26\aboveCC_CG\2022-02-26_14-48-47';
recording_folder = 'C:\Users\gjb326\Desktop\Sample Data\Fenton Si probes\RAT1_2022_02_25\aboveCC_S1HL\2022-02-25_19-13-59';
recording_folder = 'C:\Users\gjb326\Desktop\Sample Data\Fenton Si probes\RAT1_2022_02_25\aboveCC_S1HL\2022-02-25_19-13-59';
% recording_folder = 'C:\Users\gjb326\Desktop\Sample Data\Fenton Si probes\RAT2_2022_02_26\withinCC_CG\2022-02-26_15-03-02';
% recording_folder = 'C:\Users\gjb326\Desktop\Sample Data\Fenton Si probes\RAT2_2022_02_26\withinCA1\2022-02-26_16-05-47';
json_file = '\Record Node 101\experiment1\recording1\structure.oebin'; % ".oebin" is the JSON file they refer to in docs

fname = strcat(recording_folder, json_file);

input_index = 1; % DAQ board input, shouldn't need to change

% load continuous recording data
A = load_open_ephys_binary(fname, 'continuous', input_index);

% all channels should have the same bitvolt value
bv = A.Header.channels(1).bit_volts;

% NOTE: openephys auto recorded 32 channels, the probe only has 16 (channels 9:24)
D = bv*A.Data(9:24,:); 
fs = double(A.Header.sample_rate);
ts = double(A.Timestamps);
ts = (ts-ts(1))/fs;

%%%%% Channel groupings %%%%% see also "2x2 tet config.jpg"
g1 = [4 2 7 5]; % dorsomedial channels
g2 = [3 1 8 6]; % ventramedial channels
g3 = [12 10 15 13]; % dorsolateral channels
g4 = [11 9 16 14]; % ventralateral channels
channel_space = 25;% channel spacing (within group) 25 micron, hypotenuse
medial_space = 150; % medial spacing 150 micros
dorsal_space = 200; % dorsal spacing 200 micros

all_ch_order = [4 2 7 5 3 1 8 6 12 10 15 13 11 9 16 14];
channel_dev_ml = [0 -1 1 0];
channel_dev_dv = [1 0 0 -1];
group_dev_ml = [0 0 1 1];
group_dev_dv = [0 1 0 1];

%% Plot data 
dec = floor(size(D,2)*.1);
inds = 1:size(D,2); % plot all samples
inds = 1*dec:3*dec; % plot a smaller portion of the data
figure; clf;
hold on
ch2plt = [1:16];%[2 4 5 7]; % if you only want a couple
for i = 1:length(ch2plt)
    single_channel = D(ch2plt(i),(inds));
    single_channel = single_channel/500;%250; % arbitrary scaling
    plot(ts(inds), single_channel + i);    
end
axis tight
% aboveCC_CG recording has some units on channels 2 and 4 
%% Plot 3D data 
if false
    dec = floor(size(D,2)*.01);
    % inds = 1:size(D,2); % plot all samples
    inds = 61*dec:62*dec; % plot a smaller portion of the data
    figure(2); clf;
    hold on
    ch2plt = [1:16];%[2 4 5 7];% % if you only want a couple
    for i = 1:length(ch2plt)
        single_channel = D(ch2plt(i),(inds));
        single_channel = single_channel/3; % arbitrary scaling
        ch_ind = mod(i-1,4)+1;
        group_ind = floor((i-1)/4)+1;
        pos_dv = channel_dev_dv(ch_ind)*channel_space + group_dev_dv(group_ind)*dorsal_space;
        pos_dv = pos_dv + single_channel*0;
        pos_ml = channel_dev_ml(ch_ind)*channel_space + group_dev_ml(group_ind)*medial_space;
        pos_ml = pos_ml + single_channel*0;
        plot3(ts(inds), pos_ml, single_channel+pos_dv, '-');
    end
end
%% Filtering data 
dec = floor(size(D,2)*.1);
inds = 1:size(D,2); % plot all samples
inds = 0*dec+1:3*dec; % plot a smaller portion of the data
figure; clf;
hold on
ch2plt = [1:16];%[2 4 5 7]; % if you only want a couple
for i = 1:length(ch2plt)
    x = bandpass(D(i,:), [300 6000], 30000);
%     x = bandpass(D(i,:), [5 15], 30000);
%     x = bandpass(D(i,:), [50 500], 30000);
    single_channel = x(inds);
    single_channel = single_channel/250; % arbitrary scaling
    plot(ts(inds), single_channel + i);   
    drawnow
end
axis tight
















