function draw_sleap_skeleton(sleappos, edge_inds, inds2draw, plotargs)
e = edge_inds;

if nargin<3
    inds2draw = 1:size(sleappos,1);
end

if nargin<4
    plotargs = {'.-k'};
end
for idx = 1:length(inds2draw)
%      = {'-', 'Color', [.5 .2 .2]};
    for e1 = 1:size(e,2)
        xx = [sleappos(inds2draw(idx), e(1, e1), 1), sleappos(inds2draw(idx), e(2, e1), 1)];
        yy = [sleappos(inds2draw(idx), e(1, e1), 2), sleappos(inds2draw(idx), e(2, e1), 2)];
        
        plot(xx,yy,plotargs{:})
%         drawnow
%         pause(.02)
    end
end