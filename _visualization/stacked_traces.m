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
elseif nargin==4
    scale_t = varargin{1};
    line_props=varargin{2};
    line_colors = varargin{3};
end

hold on
for i = nsegs:-1:1
    if nargin==4
%     area(i+traces(i,:)*scale_t, 'FaceColor', line_colors(i,:));
%     plot(i+traces(i,:)*scale_t, line_props{:}, 'Color', 'k')
    area(i+traces(i,:)*scale_t, 'FaceColor', 'w');
    plot(i+traces(i,:)*scale_t, line_props{:}, 'Color', line_colors(i,:))
    else
    plot(i+traces(i,:)*scale_t, line_props{:})
    end
end
hold off