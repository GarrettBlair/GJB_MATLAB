function [ang_dist] = angular_distance(ang1, ang_ref)
angdiff = ang1 - ang_ref;
ang_dist = mod(angdiff + pi, 2*pi) - pi;