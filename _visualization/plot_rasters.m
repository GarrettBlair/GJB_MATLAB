%% plots all C traces of neuron structure
function plot_rasters(neuron, segs, mark_size)
spc = 1;
if ~exist('segs', 'var')
    segs = 1:size(neuron.C,1); % use all segments
end

if ~exist('mark_size', 'var')
    mark_size = 5;
end

nsegs = length(segs);
%% plot traces of neuronC
xx = [1:length(neuron.C(1,:))];
for i = 1:nsegs
    hold on;
    spks = find(neuron.S(segs(i),:) > 0);
    plot(xx(spks), ones(length(spks),1)+spc*(i-1), 'ks', 'MarkerSize', mark_size, 'MarkerFaceColor', 'k');
    hold off
end
axis([0 length(neuron.C) -1 nsegs+1])
