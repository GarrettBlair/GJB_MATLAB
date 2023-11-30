function [dff_filt, raw_filt] = ca_trace_dff_filter(raw_traces, fs, lowpass_fs, std_thresh, plotting)


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
% raw = n.C(1:50:end, :)+n.YrA(1:50:end, :);
% sm = movmedian(raw, [4,4], 2);
% r = real(highpass(sm', .1, 30))';
low = real(lowpass(raw_traces', lowpass_fs, fs))';
% high = real(highpass(raw', 10, 30))';

dff = diff(low, 1, 2);
dff = cat(2, dff(:,1), dff);

stddff = real(std_scale_rows(dff));
% stddff(stddff<2 & dff<0) = 0;
raw_filt = low;
raw_filt(stddff<std_thresh | dff<0 | low<0) = 0;
raw_filt(raw_filt<0) = 0;
dff_filt = dff;
dff_filt(raw_filt<=0) = 0;

if plotting==true
    nr = normalize_rows(raw_traces);
    nd = normalize_rows(dff_filt);
%     s = raw_traces;
%     s(dff_filt==0) = NaN;
%     figure;
%     stacked_traces(raw_traces,  .9, {'k-', 'LineWidth', 1});
%     stacked_traces(s,           .9, {'m-', 'LineWidth', 3});
%     
    figure;
    imagesc(nr + 1*(nd>0))
%     imagesc(nr + nd)
%     colormap(lbmap(256,'RedBlue'))
    colormap redblue
    drawnow
end

