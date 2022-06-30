recording_folder = 'C:\Users\gjb326\Desktop\Sample Data\Fenton Si probes\RAT2_2022_02_26\aboveCC_CG\2022-02-26_14-48-47';
json_file = '\Record Node 101\experiment1\recording1\structure.oebin'; % ".oebin" is the JSON file they refer to in docs

fname = strcat(recording_folder, json_file);

input_index = 1; % DAQ board input, shouldn't need to change

% load continuous recording data
A = load_open_ephys_binary(fname, 'continuous', input_index);

% all channels should have the same bitvolt value
bv = A.Header.channels(1).bit_volts;

% NOTE: openephys auto recorded 32 channels, the probe only has 16 (channels 9:24)
D = bv*A.Data(9:24,:); 
ts = A.Timestamps;

%%%%% Channel groupings %%%%% see also "2x2 tet config.jpg"
g1 = [4 2 7 5]; % dorsomedial channels
g2 = [3 1 8 6]; % ventramedial channels
g3 = [12 10 15 13]; % dorsolateral channels
g4 = [11 9 16 14]; % ventralateral channels
% channel spacing (within group) 25 micron, hypotenuse
% medial spacing 150 micros
% dorsal spacing 200 micros

%% Plot data 
dec = floor(size(D,2)*.1);
% inds = 1:size(D,2); % plot all samples
inds = 6*dec:7*dec; % plot a smaller portion of the data
figure;
hold on
ch2plt = [1:16];%[2 4 5 7]; % if you only want a couple
for i = 1:length(ch2plt)
    single_channel = D(ch2plt(i),(inds));
    single_channel = single_channel/250; % arbitrary scaling
    plot(ts(inds), single_channel + i);    
end
% aboveCC_CG recording has some units on channels 2 and 4 