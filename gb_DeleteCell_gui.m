function gb_DeleteCell_gui(key, axis)

figure(1);
set(gcf,'Position',[6 124 1908 859]);

%commented because this is now assigned elsewhere
%refcell = current_refcell;
frametype = 'max';
if strcmp(frametype,'max')
    %         temp = imresize(ms.maxFrame{:}, 2);
    background = ms.maxFrame{:};
elseif strcmp(frametype,'mean')
    background = ms.meanFrame{:};
else
    background = ms.minFrame{:};
end

segs = 1:size(neuron.C,1);
numsegs = length(segs);
[contours, ~] = contours_only(neuron, ms, [], .4, 0);
centers = NaN(numsegs, 2);
for i = 1:size(contours,1)
    a = squeeze(contours(i,:,:));
    props = regionprops(full(a),'Area','Centroid');
    
    area = 0;
    for j=1:length(props)
        if (props(j).Area > area)
            area = props(j).Area;
            index  = j;
        end
    end
    centers(i,1) = round(props(index).Centroid(1));
    centers(i,2) = round(props(index).Centroid(2));
    
end
refcell = 1;

button = -1000;
last_near = 1;
near_order = numsegs;
currentCell = 1;
cells2delete = [];
allcells = 1:numsegs;
xdim=neuron.options.d1; ydim=neuron.options.d2;

x = centers(currentCell, 1); y = centers(currentCell, 2);
d=sqrt( (x-centers(:,1)).^2 + (y-centers(:,2)).^2 );
[~, near_ord] = sort(d, 'ascend');
neighbors = near_ord(2:6);
while ~(button==113)
    
    selContour = squeeze(contours(currentCell, :, :));
    others = setdiff(allcells, [currentCell cells2delete neighbors']);
    otherContour = squeeze(sum(contours(others, :, :), 1));
    neighborsContour = squeeze(sum(contours(neighbors, :, :), 1));
    if any(cells2delete)
        deleteContour = squeeze(sum(contours(cells2delete, :, :), 1));
    else
        deleteContour = selContour*0;
    end
    
    %%%% REFERENCE PLOT
    %     subplot(6,6,[1:3 7:9 13:15 19:21 25:27 ]); cla
    
    
    ref_im = zeros(xdim, ydim, 3);
    
    ref_im(:, :, 1) = background./255 + deleteContour.*.25 + neighborsContour.*.25;
    ref_im(:, :, 2) = background./255 + selContour.*.8 + neighborsContour.*.25;
    ref_im(:, :, 3) = background./255 + otherContour.*.25;
    %%
    subplot(4,1,[1:3]); cla
    image(ref_im)
    hold on
    scatter(centers(cells2delete,1), centers(cells2delete,2), 5, 'MarkerFaceColor', [.8 .2 0], 'MarkerEdgeColor', 'k')
    scatter(centers(others,1), centers(others,2), 5, 'MarkerFaceColor', [.2 .8 1], 'MarkerEdgeColor', 'k')
    scatter(centers(neighbors,1), centers(neighbors,2), 10, 'MarkerFaceColor', [.6 1 .6], 'MarkerEdgeColor', 'k')
    scatter(centers(currentCell,1), centers(currentCell,2), 10, 'MarkerFaceColor', [.9 .9 .9], 'MarkerEdgeColor', 'k')
    %text(reference_sess.center(cell2index(refcell,refcol), 1)+50, reference_sess.center(cell2index(refcell,refcol), 2)+50, sprintf('%d', cell2index(refcell,refcol)), 'Color', 'r', 'FontSize', 12);
    axis image
    xlabel(['cell ' num2str(currentCell) ],'FontSize',15);
    title('d = delete; q = quit','FontSize',15);
    
    subplot(4,1,4); cla
    plot_Cs(neuron, [currentCell; neighbors], 1, 1, 1);
    
    %% USER INPUT
    
    [x,y,button]=ginput(1);
    scatter(round(x), round(y), 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
    
    %     [x y button]
    if button == 46 % '>'
        button = 29;
    elseif button == 44 % '<'
        button = 28;
    elseif button == 27 % 'esc'
        button = 113;
    elseif button == 127 % 'del'
        button = 100;
    end
    
    switch button
        case 1
            d=sqrt( (x-centers(:,1)).^2 + (y-centers(:,2)).^2 );
            [~, near_ord] = sort(d, 'ascend');
            last_near = 1;
            currentCell=near_ord(last_near);
            
            x = centers(currentCell, 1); y = centers(currentCell, 2);
            d=sqrt( (x-centers(:,1)).^2 + (y-centers(:,2)).^2 );
            [~, near_ord] = sort(d, 'ascend');
            
            fprintf('#%d    ', currentCell)
            
            try
                neighbors = near_ord(last_near+1:last_near+5);
            catch
                neighbors = near_ord(last_near+1:end);
            end
            
        case 116 % t = toggle del cells
            %%
            if any(cells2delete)
                ref_im = zeros(xdim, ydim, 3);
                ref_im(:, :, 1) = background./255 + deleteContour.*.6;
                ref_im(:, :, 2) = background./255;
                ref_im(:, :, 3) = background./255;
                clf;
                image(ref_im)
                hold on
                scatter(centers(cells2delete,1), centers(cells2delete,2), 5, 'MarkerFaceColor', [.8 .2 0], 'MarkerEdgeColor', 'k')
                %             scatter(centers(others,1), centers(others,2), 5, 'MarkerFaceColor', [.2 .8 1], 'MarkerEdgeColor', 'k')
                %             scatter(centers(selcell,1), centers(selcell,2), 7, 'MarkerFaceColor', [.6 1 .6], 'MarkerEdgeColor', 'k')
                axis image
                xlabel(['cell ' num2str(currentCell) ],'FontSize',15);
                title('Toggle (t) or Quit (q) to return to select','FontSize',15);
                button = 10000;
                while ~(button == 116 || button == 113)
                    [~, ~, button] = ginput(1);
                end
            end
        case 29 % right arrow or </, = next nearest cell
            if last_near < numsegs
                last_near = last_near+1;
                currentCell=near_ord(last_near);
                
                x = centers(currentCell, 1); y = centers(currentCell, 2);
                d=sqrt( (x-centers(:,1)).^2 + (y-centers(:,2)).^2 );
                [~, near_ord] = sort(d, 'ascend');
                
                try
                    neighbors = near_ord(last_near+1:last_near+5);
                catch
                    neighbors = near_ord(last_near+1:end);
                end
                fprintf('#%d    ', currentCell)
            end
        case 28 % left arrow or ./> = prev nearest cell
            if last_near > 1
                last_near = last_near-1;
                currentCell=near_ord(last_near);
                
                x = centers(currentCell, 1); y = centers(currentCell, 2);
                d=sqrt( (x-centers(:,1)).^2 + (y-centers(:,2)).^2 );
                [~, near_ord] = sort(d, 'ascend');
                
                try
                    neighbors = near_ord(last_near+1:last_near+5);
                catch
                    neighbors = near_ord(last_near+1:end);
                end
            end
            fprintf('#%d    ', currentCell)
        case 30 % up = next trace
            if last_near < numsegs
                last_near = last_near+1;
                currentCell=near_ord(last_near);
                
                try
                    neighbors = near_ord(last_near+1:last_near+5);
                catch
                    neighbors = near_ord(last_near+1:end);
                end
                fprintf('#%d    ', currentCell)
            end
        case 31 % down arrow = prev trace
            if last_near > 1
                last_near = last_near-1;
                currentCell=near_ord(last_near);
                
                try
                    neighbors = near_ord(last_near+1:last_near+5);
                catch
                    neighbors = near_ord(last_near+1:end);
                end
            end
            fprintf('#%d    ', currentCell)
        case 100 % d = delete cell
            %             neuron.delete(selcell);
            cells2delete = [cells2delete currentCell];
            fprintf('\\nn    #%d Added to delete list! \n cells2delete: %s\n', currentCell, num2str(cells2delete))
        case 101 % e = end cell (last)
            currentCell=numsegs(end);
            fprintf('#%d    ', currentCell)
        case 107 % k = remove cell from delete list
            if ismember(currentCell, cells2delete)
                cells2delete = setdiff(cells2delete, currentCell);
                fprintf('')
            end
        case 115 % s = first cell (last)
            currentCell=numsegs(end);
            fprintf('#%d    ', currentCell)
        case 120 % x = input cell number
            currentCell=input('\n ~~~ Type cell number to jump to: ');
            fprintf('#%d    ', currentCell)
        case 113 % q = quit
            
        otherwise
            fprintf('Unknown button! key list:\nleftClick to choose, rightArrow = next closest, leftArrow = previous closest\n')
    end
    
end
%%
subplot(4,1,[1:3]); cla
image(ref_im)
hold on
scatter(centers(cells2delete,1), centers(cells2delete,2), 5, 'MarkerFaceColor', [.8 .2 0], 'MarkerEdgeColor', 'k')
scatter(centers(others,1), centers(others,2), 5, 'MarkerFaceColor', [.2 .8 1], 'MarkerEdgeColor', 'k')
scatter(centers(neighbors,1), centers(neighbors,2), 10, 'MarkerFaceColor', [.6 1 .6], 'MarkerEdgeColor', 'k')
scatter(centers(currentCell,1), centers(currentCell,2), 10, 'MarkerFaceColor', [.9 .9 .9], 'MarkerEdgeColor', 'k')
%text(reference_sess.center(cell2index(refcell,refcol), 1)+50, reference_sess.center(cell2index(refcell,refcol), 2)+50, sprintf('%d', cell2index(refcell,refcol)), 'Color', 'r', 'FontSize', 12);
axis image
xlabel(['cell ' num2str(currentCell) ],'FontSize',15);
title('d = delete; q = quit','FontSize',15);

subplot(4,1,4); cla
plot_Cs(neuron, [currentCell; neighbors], 1, 1, 1);

%['DONE WITH MATCHING!']

% close(gcf);
function [spatial, ms_spatial] = contours_only(neuron, ms, segs, spatial_thresh, pad)
% ms_spatial? second output
d1s = neuron.options.d1;
d2s = neuron.options.d2;

if ~exist('segs', 'var') || isempty(segs)
    segs = 1:size(neuron.C,1);
end
numSegs = length(segs);
% threshold cell bodies
spatial = zeros(numSegs, d1s, d2s);
% pad = 60;

ms_spatial = zeros(numSegs, pad + ms.height/2, pad + ms.width/2);


maxh = floor(max(ms.hShift(:,1))/2);
minh = floor(abs(min(ms.hShift(:,1)))/2);
maxw = floor(max(ms.wShift(:,1))/2);
minw = floor(abs(min(ms.wShift(:,1)))/2);

nsize = size(reshape(neuron.A(:,1),d1s,d2s));
hdiff = ms.height/2 - (nsize(1) + maxh + minh);
wdiff = ms.width/2 - (nsize(2) + maxw + minw);


neuron.A = full(neuron.A);
for i = 1:length(segs)
    seg = segs(i);
    c = max(neuron.A(neuron.A(:,seg)>0, seg))*spatial_thresh;
    bad = neuron.A(:,seg) < c;
    n = neuron.A(:,seg)/max(neuron.A(:,seg));
    %     nc(i) = sum(n>0);
    %     if
    n(bad) = 0;
    n = reshape(n,d1s,d2s);
    %     n = get_contiguous_contour(n);
    
    %     n = restrict_contour_size(n, pix_limit);
    spatial(seg,:,:) = n;%%%%
    
    ns = [zeros(d1s, maxw + pad/2), n, zeros(d1s, minw + pad/2)]; %
    ns = [zeros(maxh + pad/2, size(ns, 2)); ns; zeros(minh + pad/2, size(ns, 2))];
    if hdiff ~= 0 % pad zeros
        ns = [ns; zeros(1, size(ns,2))];
    end
    if wdiff ~= 0 % pad zeros
        ns = [ns, zeros(size(ns,1), 1)];
    end
    ms_spatial(seg,:,:) = ns;%%%%
end
end

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
end

end