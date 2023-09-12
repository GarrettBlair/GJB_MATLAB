function [ang_dist, inside_ang_ref] = angular_distance(ang1, ang_ref, bleed_size)
if nargin<3
    bleed_size=0;
end
zonex = [0 length(ang1)+1 length(ang1)+1 0 0];
zoney = [ang_ref-bleed_size ang_ref-bleed_size ang_ref+bleed_size ang_ref+bleed_size ang_ref-bleed_size];
inside_ang_ref = false(size(ang1));
inside_ang_ref(:) = inpolygon(1:length(ang1), ang1, zonex, zoney);
if bleed_size==0
    angdiff = ang1 - ang_ref;
    ang_dist = mod(angdiff + pi, 2*pi) - pi;
    ang_dist = abs(ang_dist);
else
    angdiff = ang1 - ang_ref-bleed_size;
    ang_dist1 = mod(angdiff + pi, 2*pi) - pi;
    ang_dist1 = abs(ang_dist1);
    angdiff = ang1 - ang_ref+bleed_size;
    ang_dist2 = mod(angdiff + pi, 2*pi) - pi;
    ang_dist2 = abs(ang_dist2);
    ang_dist = min(ang_dist1, ang_dist2);
%     ang_dist(ang_dist<=bleed_size) = -.5;
end
plotting = false;
if plotting == true
figure(20); clf; subplot(211)
plot(ang1); hold on; 
rectangle('Position', [0 ang_ref-bleed_size length(ang1) bleed_size*2]); 
plot(find(inside_ang_ref), ang1(inside_ang_ref), 'ro');
subplot(212)
hold on
% plot(ang_dist1); 
% plot(ang_dist2);
plot(ang_dist);
rectangle('Position', [0 0 length(ang1) 0]); 
plot(find(inside_ang_ref), ang_dist(inside_ang_ref), 'ro');
end

