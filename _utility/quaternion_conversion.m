%%
q0 = o(:,1);
q1 = o(:,2);
q2 = o(:,3);
q3 = o(:,4);

xyz = real([atan2( 2*(q0.*q1 + q2.*q3), 1 - 2*(q1.^2 + q2.^2) ), asin(2*(q0.*q2 - q3.*q1)), atan2( 2*(q0.*q3 + q1.*q2), 1 - 2*(q2.^2 + q3.^2) )]);
figure(99); clf
subplot(131)
plot(time, xyz(:,1), 'k.')
subplot(132)
plot(time, xyz(:,2), 'k.')
subplot(133)
plot(time, xyz(:,3), 'k.')

%%
subplot(121)
hSurface = surf(peaks(20));
direction = [1 1 1];
subplot(122)

hSurface = surf(peaks(20));
for i = 1:5:5000
    rotate(hSurface,[1 0 0],xyz(i,1)');
    rotate(hSurface,[0 1 0],xyz(i,2)')
    rotate(hSurface,[0 0 1],xyz(i,3)')
%     rotate(hSurface,[1 1 1],xyz(i,:)')
    axis([ -20 20 -20 20 -20 20])
    drawnow
end