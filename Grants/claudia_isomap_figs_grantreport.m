

frxn = frx;
X_data = frxn';
D = squeeze(sqrt(sum(bsxfun(@minus,X_data,reshape(X_data',1,size(X_data,2),size(X_data,1))).^2,2)));
options.dims = 1:3;
[Y, R, E] = Isomap(D, 'k', 50, options); 

frxn = zscore(frx,1,2);
X_data = frxn';
D = squeeze(sqrt(sum(bsxfun(@minus,X_data,reshape(X_data',1,size(X_data,2),size(X_data,1))).^2,2)));
options.dims = 1:3;
[Yn, R, E] = Isomap(D, 'k', 50, options); 
%%
[nsegs, nt] = size(frx);
ns = 1:nt;
% inds = [0 k-1; k k+k2; k2:k+1 nt];
inds1 = ns<k;
inds2 = ns>=k & ns<(k+k2);
inds3 = ns>=(k+k2);
figure(6); clf; hold on;
% plot(Y.coords{2}(1,:), Y.coords{2}(2,:))
a = .5;
scatter(find(inds1), Y.coords{1}(1,inds1), 100, '.', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', a)
scatter(find(inds2), Y.coords{1}(1,inds2), 100, '.', 'MarkerEdgeColor', 'b', 'MarkerEdgeAlpha', a)
scatter(find(inds3), Y.coords{2}(1,inds3), 100, '.', 'MarkerEdgeColor', [.5 .5 .5], 'MarkerEdgeAlpha', a)

figure(7); clf; hold on;
% plot(Y.coords{2}(1,:), Y.coords{2}(2,:))
scatter(Y.coords{2}(1,inds1), Y.coords{2}(2,inds1), 100, '.', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', a)
scatter(Y.coords{2}(1,inds2), Y.coords{2}(2,inds2), 100, '.', 'MarkerEdgeColor', 'b', 'MarkerEdgeAlpha', a)
scatter(Y.coords{2}(1,inds3), Y.coords{2}(2,inds3), 100, '.', 'MarkerEdgeColor', [.5 .5 .5], 'MarkerEdgeAlpha', a)


figure(8); clf; hold on;
% plot3(Y.coords{3}(1,:), Y.coords{3}(2,:), Y.coords{3}(3,:), 'k-')
scatter3(Y.coords{3}(1,inds1), Y.coords{3}(2,inds1), Y.coords{3}(3,inds1), 100, '.', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', a)
scatter3(Y.coords{3}(1,inds2), Y.coords{3}(2,inds2), Y.coords{3}(3,inds2), 100, '.', 'MarkerEdgeColor', 'b', 'MarkerEdgeAlpha', a)
scatter3(Y.coords{3}(1,inds3), Y.coords{3}(2,inds3), Y.coords{3}(3,inds3), 100, '.', 'MarkerEdgeColor', [.5 .5 .5], 'MarkerEdgeAlpha', a)

%%
[nsegs, nt] = size(frx);
ns = 1:nt;
% inds = [0 k-1; k k+k2; k2:k+1 nt];
inds1 = ns<k;
inds2 = ns>=k & ns<(k+k2);
inds3 = ns>=(k+k2);
figure(16); clf; hold on;
% plot(Y.coords{2}(1,:), Y.coords{2}(2,:))
a = .5;
scatter(find(inds1), Yn.coords{1}(1,inds1), 100, '.', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', a)
scatter(find(inds2), Yn.coords{1}(1,inds2), 100, '.', 'MarkerEdgeColor', 'b', 'MarkerEdgeAlpha', a)
scatter(find(inds3), Yn.coords{2}(1,inds3), 100, '.', 'MarkerEdgeColor', [.5 .5 .5], 'MarkerEdgeAlpha', a)

figure(17); clf; hold on;
% plot(Yn.coords{2}(1,:), Yn.coords{2}(2,:))
scatter(Yn.coords{2}(1,inds1), Yn.coords{2}(2,inds1), 100, '.', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', a)
scatter(Yn.coords{2}(1,inds2), Yn.coords{2}(2,inds2), 100, '.', 'MarkerEdgeColor', 'b', 'MarkerEdgeAlpha', a)
scatter(Yn.coords{2}(1,inds3), Yn.coords{2}(2,inds3), 100, '.', 'MarkerEdgeColor', [.5 .5 .5], 'MarkerEdgeAlpha', a)

figure(18); clf; hold on;
% plot3(Yn.coords{3}(1,:), Yn.coords{3}(2,:), Yn.coords{3}(3,:), 'k-')
scatter3(Yn.coords{3}(1,inds1), Yn.coords{3}(2,inds1), Yn.coords{3}(3,inds1), 100, '.', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', a)
scatter3(Yn.coords{3}(1,inds2), Yn.coords{3}(2,inds2), Yn.coords{3}(3,inds2), 100, '.', 'MarkerEdgeColor', 'b', 'MarkerEdgeAlpha', a)
scatter3(Yn.coords{3}(1,inds3), Yn.coords{3}(2,inds3), Yn.coords{3}(3,inds3), 100, '.', 'MarkerEdgeColor', [.5 .5 .5], 'MarkerEdgeAlpha', a)
%%
axis([-30 30 -30 30 -30 30])
for i = 1:300
    %%
    set(gca, 'View', [45+i 30], 'XTick', [], 'YTick', [], 'ZTick', [], 'Box', 'on')
    drawnow
    pause(.01)
end