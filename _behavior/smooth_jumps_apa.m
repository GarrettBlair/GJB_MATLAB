function [x, y] = smooth_jumps_apa(x, y, ts, pixpercm, jump_thresh, return_thresh, twin_seconds)

nanind = (isnan(x) & isnan(y));
xn = interp1(ts(~nanind), x(~nanind), ts(nanind), 'nearest');
yn = interp1(ts(~nanind), y(~nanind), ts(nanind), 'nearest');
x(nanind) = xn; y(nanind) = yn;

% jump_thresh=10;%*cmperpix; % cm threshold for jumps
% return_thresh=10;%*cmperpix; % cm threshold for jumps

twin= twin_seconds*round(1/nanmedian(abs(diff(ts)))); % time win to look
dxy = sqrt(diff(x).^2 + diff(y).^2)/pixpercm;
bads = find(dxy>jump_thresh)+1;
fprintf('\nSmoothing %d jumps \n', length(bads))
startidx=1;
while_breaker = tic();
figure; plot(dxy)
counter=0;

figure; % 
hold on; plot3(x,y,ts, 'Color', [.4, .4, .4])
for i = 1:length(bads)
    hold on; 
    ii = bads(i)-1:bads(i)+twin;
    ii = ii(ii>0 & ii<length(x));
    plot3(x(ii), y(ii), ts(ii), 'r.'); 
end

while ~isempty(bads) % ind = 1:length(bads)
%                 badx=x(bads(1):bads(1)+twin);
%                 bady=y(bads(1):bads(1)+twin);
    goodx=x(bads(1)-1); % :end);
    goody=y(bads(1)-1);
    
    ii = bads(1)-1:bads(1)+twin;
    ii = ii(ii>0 & ii<length(x));

    dxy_prev = sqrt((goodx - x(ii)).^2 +...
                    (goody - y(ii)).^2)/pixpercm;
    next_good = find(dxy_prev < return_thresh,1)+bads(1);
    
    temp = sqrt(diff(x).^2 + diff(y).^2)/pixpercm;
    dxy = sqrt(diff(x(startidx:end)).^2 + diff(y(startidx:end)).^2)/pixpercm;
%     hold on; 
%     plot(temp+100*counter, 'k')
%     plot(bads, temp(bads)+100*counter, 'r.')
%     plot(bads(1), temp(bads(1))+100*counter, 'ro')
%     plot(bads(1)+twin, temp(bads(1))+100*counter, 'ro')
    if ~isempty(next_good)
        x(bads(1):next_good) = linspace(goodx, x(next_good), next_good-bads(1)+1);
        y(bads(1):next_good) = linspace(goody, y(next_good), next_good-bads(1)+1);
        startidx = next_good;
        clr='g';
    else
        dxy_prev = sqrt((x(bads(1)-1) - x(bads(1):end)).^2 +...
                        (y(bads(1)-1) - y(bads(1):end)).^2)/pixpercm;
        startidx = find(dxy_prev < jump_thresh/pixpercm,1)+bads(1);
        next_good = startidx;
        clr='r';
    end
    plot3(x(bads(1):next_good), y(bads(1):next_good), ts(bads(1):next_good), '*', 'Color', clr)
    temp = sqrt(diff(x).^2 + diff(y).^2)/pixpercm;
    dxy = sqrt(diff(x(startidx:end)).^2 + diff(y(startidx:end)).^2)/pixpercm;
%     hold on; 
%     plot(temp+100*counter, 'k')
%     plot(bads(1), temp(bads(1))+100*counter, 'ro')
%     plot(startidx, dxy(1)+100*counter, 'x', 'Color', clr)
%     plot(startidx:startidx+length(dxy)-1, dxy +100*counter, '--')
    bads = find(dxy>jump_thresh)+startidx;
    counter = counter+1;
    if toc(while_breaker)>300
        warning('Got stuck in smooth jumps!')
        bads=[];
    end
end