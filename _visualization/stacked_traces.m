function stacked_traces(traces, varargin)
% example: stacked_traces(C, 1.5, {'r-', 'LineWidth', 20})
nsegs = size(traces,1);
% traces = normalize_rows(traces);
% traces(isnan(traces))=0;
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
    area_face =varargin{3};
end

hold on
for i = nsegs:-1:1
    t = traces(i,:);
    nanind = isnan(t);
    if any(~isnan(t))
    t(isnan(t))=0;
    t2= i+t*scale_t;
    t2(nanind) = NaN;
    if nargin==4
%     area(i+traces(i,:)*scale_t, 'FaceColor', line_colors(i,:));
%     plot(i+traces(i,:)*scale_t, line_props{:}, 'Color', 'k')
    area(t2, 'FaceColor', area_face);
    end
    plot(t2, line_props{:})
    end
end
axis tight
hold off