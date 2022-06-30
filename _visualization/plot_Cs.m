%% plots all C traces of neuron structure
function traces = plot_Cs(neuron, segs, color_flag, raw_flag, spike_flag)
spc = 1.2;
gain = 1;
if ~exist('segs', 'var') || isempty(segs)
    segs = 1:size(neuron.C,1); % use all segments
end
nsegs = length(segs);

if color_flag == true % use custom color gradient
    ca = .3;
    gg = linspace(ca, 1-ca, nsegs);
    rr = ones(1, nsegs)*.1;
    bb = linspace(1-ca, ca, nsegs);
    colors = [rr' gg' bb'];
    light_colors = [rr'+.2 gg' bb'+.2];
elseif ~exist('color_flag', 'var') || color_flag == false % false or empty == black coloring
    colors = zeros(nsegs, 3);
    light_colors = zeros(nsegs, 3);
else % random coloring
    rng(pi); % consistentcy
    colors = rand(nsegs, 3);
    light_colors = rand(nsegs, 3);
end

if ~exist('raw_flag', 'var')
    raw_flag = false; % use denoised/demixed traces
end

if ~exist('spike_flag', 'var')
    spike_flag = false; % use denoised/demixed traces
end

%% plot traces of neuronC
hold on;
switch raw_flag
    case 0
        for i = nsegs:-1:1
            norm_trace = (neuron.C(segs(i),:)/max(neuron.C(segs(i),:)));
            plot(norm_trace*gain+spc*(i-1), 'Color', colors(i,:), 'LineWidth', 1.5);
            if spike_flag
                spk = find(neuron.S(segs(i),:)>0);
                plot(spk, norm_trace(spk)*gain+spc*(i-1), 'r.')
            end
            traces(i,:) = norm_trace;
        end
        
    case 1
        for i = nsegs:-1:1
            norm_trace = (neuron.C_raw(segs(i),:)/max(abs(neuron.C_raw(segs(i),:))));
%             plot((norm_trace-mean(norm_trace))*gain+spc*(i-1), 'Color', light_colors(i,:), 'LineWidth', .5);
            plot((norm_trace)*gain+spc*(i-1), 'Color', light_colors(i,:), 'LineWidth', 1.5);
            if spike_flag
                spk = find(neuron.S(segs(i),:)>0);
                plot(spk, norm_trace(spk)*gain+spc*(i-1), 'r.')
            end
            traces(i,:) = norm_trace;
        end
end

hold off

% axis([0 length(neuron.C) -.5*spc nsegs+((spc-1)*gain)])
axis tight
set(gca, 'YTick', 0:spc:(nsegs-1)*spc, 'YTickLabel', segs)
