w = 12.5;
h = 53; % total heigt
a = (atan(h/w));

h2=12;
w2 = cos((a))*h2;
t = (w-w2)/cos(a);
t = (h2)/sin(a);

close all

figure(123); clf; hold on
plot([0 w], [0, 0])
plot([0 0], [0, h])
plot([0 w-w2], [h2, h2])
plot([0 w], [h, 0])
plot([w-w2 w], [t, 0])
axis equal

%%
w-w2


figure(124); clf; hold on;
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8.5, 11], 'PaperUnits', 'Inches', 'PaperSize', [8.5, 11])
circle(0,0, w-w2)

circle(0,0, w-w2+t)
axis equal
%%
function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end