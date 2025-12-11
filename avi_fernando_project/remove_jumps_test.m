%%
x = led_x;
y = led_y;
figure(13); clf; 
plot(x,'k.')
% x = movmedian(x, 5);
% y = movmedian(y, 5);
hsmall = gausswin(7); hsmall = hsmall./sum(hsmall(:));

naninds = isnan(x);
naninds = conv(naninds, ones(3,1), 'same') ==1;
x(naninds) = NaN; y(naninds) = NaN;
naninds = isnan(x);
x(naninds)= interp1(find(~naninds), x(~naninds), find(naninds), 'linear');
y(naninds)= interp1(find(~naninds), y(~naninds), find(naninds), 'linear');

x = conv(x, hsmall, 'same');
y = conv(y, hsmall, 'same');

hold on; plot(x, 'r.')
spd = speed_calc_gb(x, y);
yyaxis('right')
plot(spd)


%%
ha = led_dist;
ha(isnan(ha))= interp1(find(~isnan(ha)), ha(~isnan(ha)), find(isnan(ha)), 'nearest');
dd = abs(diff(ha));
dd = [dd(1); dd];

ha = headangle;
ha(isnan(ha))= interp1(find(~isnan(ha)), ha(~isnan(ha)), find(isnan(ha)), 'nearest');
da = abs(diff(unwrap(deg2rad(ha))));
da = [da(1); da];

dh = da.*dd;


dnew = spd;
ha = unwrap(deg2rad(ha));
hnew = ha;
xnew = x; ynew = y;
thresh = 40;
a = 1;
while any(dnew >= thresh)
    departs = find(dnew >= thresh, 1);
    % dists = hnew(departs:end)-hnew(departs-1) <= thresh ;
    dists = speed_calc_gb(x(departs-1:end), y(departs-1:end));

    returns = departs + find(dists, 1);
    if ~isempty(returns)
        disp('in')
        % hnew(departs:returns-1) = NaN;
        xnew(departs:returns-1) = NaN;
        ynew(departs:returns-1) = NaN;
        % dnew = abs(diff(unwrap(deg2rad(hnew))));
        % dnew = [dnew(1); dnew];
        dnew = speed_calc_gb(xnew, ynew);

    else
        warning('t')
        break
    end

end

figure(13); clf; 
subplot(2,1,1)
hold on; plot(spd); plot(dnew)
subplot(2,1,2)
hold on; plot(ha); plot(hnew)