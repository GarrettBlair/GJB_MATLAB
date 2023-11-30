%%
clc
freq = 14; % days between meetings
names = { 'Simon' 'Garrett' 'Claudia' 'Darryl' 'Jenessa' 'Jose' 'Gino' 'Andy'};
nnames = length(names);
n_meets = nnames*1; % number of meetings to schedule
% ords = 1:nnames; % randperm(length(names));
ords = randperm(nnames);
year_start = 2023;
month_start = 12;
day_start = 6;
skips = [31 1 2024]; % day month year for skipping
% just print dates
count = 0;
i =0;
while i < n_meets
    i=i+1;
    day = datetime(year_start, month_start, day_start + (i-1)*freq );
    isskipped = false;
    for j = 1:size(skips,1)
        dayskip = datetime(skips(j,3), skips(j,2), skips(j,1));
        if day == dayskip
            isskipped = true;
        end
    end
    if ~isskipped
        ind = ords( mod(count, nnames)+1 );
        count = count+1;
        speaker = names{ind};
    else
        n_meets = n_meets+1;
        speaker = '~skipped~';
    end
    dd = sprintf('%s', day);
    s = strfind(dd, '-');
    fprintf('%s - %s - %s    ->    %s\n', dd(s(1)+1:s(2)-1), dd(1:s(1)-1), dd(s(2)+1:end), speaker)
end
fprintf('\n\n')

% print dates and speakers
for i = 1:n_meets
    day = datetime(year_start, month_start, day_start + (i-1)*freq );
    ind = ords( mod(i, nnames)+1 );
    speaker = names{ind};
    dd = sprintf('%s', day);
    s = strfind(dd, '-');
    fprintf('%s - %s - %s\n', dd(s(1)+1:s(2)-1), dd(1:s(1)-1), dd(s(2)+1:end))
%     fprintf('%s - %s - %s    ->    %s\n', dd(s(1)+1:s(2)-1), dd(1:s(1)-1), dd(s(2)+1:end), speaker)
end
