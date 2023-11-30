function [sub_str] = bin_session_data(ms, integration_time, frame_string, params)

occupancy_thresh    = params.occupancy_thresh;
visits_thresh       = params.min_samples;

if any(strfind(frame_string, 'polar'))
    binsx = params.yaw_bins;
    binsy = params.rho_bins;
    if contains(frame_string, 'room')
        xo = ms.room.x;
        yo = ms.room.y;
    elseif contains(frame_string, 'arena')
        xo = ms.arena.x;
        yo = ms.arena.y;
    end
    [xo, yo] = cart2pol(xo,yo);
else % euclidean
    binsx = params.pos_bins;
    binsy = params.pos_bins;
    if strcmp(frame_string, 'room')
        xo = ms.room.x;
        yo = ms.room.y;
    elseif strcmp(frame_string, 'arena')
        xo = ms.arena.x;
        yo = ms.arena.y;
    end
end
spks = ms.spks>0;
dt = ms.dt;
to = ms.timestamps./1000;
if isfield(ms, 'goodFrames') && any(~ms.goodFrames)
    spks(:, ~ms.goodFrames) = 0;
end
baddt = (dt./median(dt)) >= 3;
dt(baddt) = median(dt);

[spks_bin, group] = bin_spks_time(spks, integration_time, to, false);

spks_bin = normalize_rows(spks_bin);

x_av = average_spks_time(xo', integration_time, to, false, 'mean');
y_av = average_spks_time(yo', integration_time, to, false, 'mean');
t_av = average_spks_time(to', integration_time, to, false, 'mean');
[dt_av, ~] = bin_spks_time(dt', integration_time, to, false);


[vmap, countmap, xbin, ybin]        = make_occupancymap_2D(x_av, y_av, dt_av, binsx, binsy);

vmap(countmap<visits_thresh) = NaN;
vmap(vmap<occupancy_thresh) = NaN;

sub_str = [];
if any(strfind(frame_string, 'polar'))
    sub_str.type      = 'polar';
    sub_str.theta     = x_av';
    sub_str.rho       = y_av';
    sub_str.thetabin  = xbin;
    sub_str.rhobin    = ybin;
else
    sub_str.type    = 'euclid';
    sub_str.x       = x_av';
    sub_str.y       = y_av';
    sub_str.xbin    = xbin;
    sub_str.ybin    = ybin;
end
sub_str.t                   = t_av';
sub_str.dt                  = dt_av';
sub_str.spks_bin            = spks_bin;
sub_str.vmap                = vmap;
sub_str.countmap            = countmap;
sub_str.spks_bin_group      = group;