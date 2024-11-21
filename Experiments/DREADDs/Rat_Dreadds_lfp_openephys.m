% C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\Experiments\Si_probe_CC\plotting_si_CC.m
recording_folder = 'J:\ltp\HPC 24521 post DCZ\2024-10-29_16-34-15\'; % date recording file
json_file = 'Record Node 101\experiment1\recording1\structure.oebin'; % ".oebin" is the JSON file they refer to in docs

fname = strcat(recording_folder, json_file);

input_index = 1; % DAQ board input, shouldn't need to change

% load continuous recording data
A = load_open_ephys_binary(fname, 'continuous', input_index);
