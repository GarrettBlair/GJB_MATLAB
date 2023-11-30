function [p1, p2] = ipos_patch(ipos)

yy = [0, nanmean(ipos,1), 0]; yy(isnan(yy)) = 0; yy(yy<=0) = 0;
xx = [1, 1:size(ipos,2), size(ipos,2)]; 
p1 = patch(xx,yy, 'red', 'EdgeColor', 'none');
yy = [0, nanmean(ipos,1), 0]; yy(isnan(yy)) = 0; yy(yy>=0) = 0;
p2 = patch(xx,yy, 'blue', 'EdgeColor', 'none');

end