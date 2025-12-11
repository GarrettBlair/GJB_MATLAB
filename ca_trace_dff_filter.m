function [dff_filt, raw_filt, dff, low] = ca_trace_dff_filter(raw_traces, fs, lowpass_fs, std_thresh, plotting)


if nargin < 5 || isempty(plotting)
    plotting = false;
end
if nargin < 4 || isempty(std_thresh)
    std_thresh = 1.5;
end
if nargin < 3 || isempty(lowpass_fs)
    lowpass_fs = .1;
end
if nargin < 2 || isempty(fs)
    fs = 30;
end

%% pad the traces to remove edge effects
[ns,nt] = size(raw_traces);
npad = round(1/lowpass_fs);
hpad = round(.5/lowpass_fs);
% p1 = mean(raw_traces(:,1:npad), 2);
% p2 = mean(raw_traces(:,end:-1:end-npad), 2);
p1 = raw_traces(:,1);
p2 = raw_traces(:,end);
% valid = cat(2, zeros(1, 1*hpad), ones(1,nt), zeros(1, 3*hpad))==1;
valid = cat(2, zeros(1, 2*hpad), ones(1,nt), zeros(1, 2*hpad))==1;
padded_traces = cat(2, p1*ones(1, npad), raw_traces, p2*ones(1, npad));

% low pass filter the trace
low = real(lowpass(padded_traces', lowpass_fs, fs))';
low = low(:, valid);
% low = low(:, npad:end);
% low(:, 1:hpad) = low(:, hpad+1)*ones(1,hpad);
% high = real(highpass(raw', 10, 30))';

dff = diff(low, 1, 2);
dff = cat(2, dff(:,1), dff);

stddff = real(std_scale_rows(dff));
% stddff(stddff<std_thresh) = 0;
% dff(dff<0) = 0;
% low(low<0) = 0;
% stddff(stddff<2 & dff<0) = 0;
raw_filt = low;
raw_filt(stddff<std_thresh | dff<=0 | low<=0) = 0;
raw_filt(raw_filt<0) = 0;
dff_filt = dff;
dff_filt(raw_filt<=0) = 0;

% only take the positive values
dff_acel = diff(dff_filt, 1, 2);
dff_acel = cat(2, dff_acel(:,1), dff_acel);
dff_filt(dff_acel<0) = 0;

if plotting==true
    nr = normalize_rows(raw_traces);
    s = nr;%raw_traces;
    s(dff_filt==0) = NaN;
    figure;
    try
        stacked_traces(nr(1:1:end,:),  .9, {'k-', 'LineWidth', 1});
        stacked_traces(s(1:1:end,:),   .9, {'m.-', 'LineWidth', 1});
    catch
        subplot(2,1,1)
        imagesc(nr)
        subplot(2,1,2)
        imagesc(s)
    end
    drawnow;

end

