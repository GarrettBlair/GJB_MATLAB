function [right_bias, out_pref] = smart_task_right_bias(room, arena, cut, angsep, plotting)
%     cut = 35;
%     angsep = pi/4;

strs = {'room', 'arena'};
right_bias = NaN(1,2);
out_pref = NaN(1,2);
% if plotting==true; figure; end
for i = 1:2
    x = eval(sprintf('%s.x;', strs{i}));
    y = eval(sprintf('%s.y;', strs{i}));
    [t,r] = cart2pol(x, y);
    out_pref(i) = nansum(r>cut)/nansum(r>=0);

    quad_l = (t>=(pi-angsep)) | (t<=-(pi-angsep));
    quad_r = (t<=angsep & t>=0) | (t>=-angsep & t<=0);
    out_r = quad_r==true & r>cut;
    out_l = quad_l==true & r>cut;
    right_bias(i) = (sum(out_r))/(sum(out_r | out_l));
%     right_bias(i) = (sum(out_r)-sum(out_l))/(sum(out_r)+sum(out_l));
%     right_bias(i) = sum(out_r) / sum(r>cut);
    if plotting==true
        subplot(1,2,i)
        plot(x,y, 'k'); hold on
        title(strs{i})
        scatter(x(out_r), y(out_r), 'r.')
        text(35,35, sprintf('%5d', sum(out_r)), 'Color', 'r')
        scatter(x(out_l), y(out_l), 'c.')
        text(-40,35, sprintf('%5d', sum(out_l)), 'Color', 'c')
    end
end

end