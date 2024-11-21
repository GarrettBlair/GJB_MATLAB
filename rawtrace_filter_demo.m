sampling_rate = 30;
low_pass = .1;
std_thresh = 1.5; % standard deviation
plotting   = true; % standard deviation
[dff_filterd, ~] = ca_trace_dff_filter(raw_traces, sampling_rate, low_pass, std_thresh, plotting);
