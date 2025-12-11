function sleap_out = correct_sleap_jumps(sleap)
sleap_out = sleap;
%
np = length(sleap_out.node_names);
nf = length(sleap_out.instance_scores);
debug = false;
if debug
close all
figure(101); clf
end

thresh = zeros(np, np);
og_bad = false(nf,1);
for n_thresh = (np-2):-1:0%:np-1
    pd = zeros(nf, np, np);
    mxy_score = sleap_out.instance_scores/length(sleap_out.node_names); % squeeze(min(sleap.point_scores(:,[head_idx, body_idx, tail_idx],:), [], 2));
    for ni = 1:np-1
        for nj = ni+1:np
            a = squeeze(sleap_out.tracks(:, ni, :));
            b = squeeze(sleap_out.tracks(:, nj, :));
            d = sqrt( (a(:,1) - b(:,1)).^2 + (a(:,2) - b(:,2)).^2);
            if n_thresh == np-2
                thresh(ni,nj) = nanmean(d) + 4*nanstd(d);
                thr = thresh(ni,nj);
                og_bad(d>thr) = true;
            else
                thr = thresh(ni,nj);
            end
            pd(:,ni,nj) = d>thr;
            pd(:,nj,ni) = d>thr;
            if debug
                figure(101)
                subplot(np,np, np*(ni-1)+nj)
                hold on
                plot(d+300*n_thresh);
                plot(find(d>thr), d(d>thr)+300*n_thresh, 'ro');
            end
        end
    end
    for ni = 1:np
        bads = sum(squeeze(pd(:, ni, :)),2)>n_thresh;
        y = squeeze(sleap_out.tracks(:, ni, 2));
        x = squeeze(sleap_out.tracks(:, ni, 1));
        y(bads) = NaN;
        x(bads) = NaN;
        isbad = isnan(bads);
        
        y(bads) = interp1(find(~bads), y(~bads), find(bads), 'linear');
        x(bads) = interp1(find(~bads), x(~bads), find(bads), 'linear');
        sleap_out.tracks(:, ni, 2) = y;
        sleap_out.tracks(:, ni, 1) = x;
    end
end
%%
if debug
    og_bad_exten = find(conv(og_bad, ones(5,1), 'same')>0);
    tail = 10;
    x = squeeze(sleap_out.tracks(:, :, 1));
    y = squeeze(sleap_out.tracks(:, :, 2));
    x2 = squeeze(sleap.tracks(:, :, 1));
    y2 = squeeze(sleap.tracks(:, :, 2));
    clrs = plasma(np*2);
    clrs = clrs(floor(np/2):floor(np/2)+np, :);
    figure(); hold on
    mx = mean(x,2);
    my = mean(y,2);
    for i = 1:1:length(og_bad_exten)
        %%
        clf; %hold on
        f = og_bad_exten(i);
        ftail = f-tail:f;
        ftail = ftail(ftail>0 & ftail<=nf);
        for ni = 1:np
            subplot(121); hold on
            plot(x2(ftail, ni), y2(ftail,ni), 'Color', clrs(ni,:))
            scatter(x2(f, ni), y2(f,ni), 30, '.', 'MarkerEdgeColor', clrs(ni,:))
            subplot(122); hold on
            plot(x(ftail, ni), y(ftail,ni), 'Color', clrs(ni,:))
            scatter(x(f, ni), y(f,ni), 30, '.', 'MarkerEdgeColor', clrs(ni,:))
        end
        for ni = 1:length(sleap_out.edge_inds)
            a = sleap_out.edge_inds(1, ni)+1;
            b = sleap_out.edge_inds(2, ni)+1;
            subplot(122); hold on
            plot([x(f, a) x(f, b)], [y(f, a) y(f, b)], 'Color', 'k')
            subplot(121); hold on
            plot([x2(f, a) x2(f, b)], [y2(f, a) y2(f, b)], 'Color', 'k')
        end
        subplot(121); hold on
        axis([0 640 0 480])
        axis
        subplot(122); hold on
        axis([0 640 0 480])
        drawnow
        pause(.05)
    end
end









