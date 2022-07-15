function stacked_traces(traces, varargin)
% example: stacked_traces(C, 1.5, {'r-', 'LineWidth', 20})
nsegs = size(traces,1);
% traces = normalize_rows(traces);
traces(isnan(traces))=0;
if nargin==1
    scale_t = 1;
    line_props={'k-', 'LineWidth', 1};
elseif nargin==2
    scale_t = varargin{1};
    line_props={'k-', 'LineWidth', 1};
elseif nargin==3
    scale_t = varargin{1};
    line_props=varargin{2};
end

hold on
for i = 1:nsegs
    plot(i-1+traces(i,:)*scale_t, line_props{:})
end
hold off