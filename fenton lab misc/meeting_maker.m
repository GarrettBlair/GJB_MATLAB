%%
clc
freq = 7; % days between meetings
% names = {'Garrett' 'Claudia' 'Darryl' 'Jenessa' 'Jose' 'Diego' 'Andy'};
% names = {'Andy' 'Darryl' 'Diego' 'Jose' 'Claudia' 'Garrett' 'Jenessa', 'EunHye', 'Heng'};
names = {'Andy' 'Darryl' 'Diego' 'Jose' 'Claudia' 'Garrett' 'JiYeon' 'Jenessa', 'Heng', 'Luke'};
nnames = length(names);
n_meets = nnames*2; % number of meetings to schedule
% ords = 1:nnames; % randperm(length(names));
ords = randperm(nnames);
year_start = 2025;
month_start = 1;
day_start = 15;
skips = [1 31; 10 30; 11 27]; % day month year for skipping; MONTH DAY format
% just print dates 
fname = []; % 'C:\Users\gjb326\Documents\MATLAB\MeetingMaker.csv';
if~isempty(fname)
    count = 0;
    i =0;
    f = fopen(fname, 'w');
    fwrite(f, sprintf('%s, %s\n', 'Date', 'Presenter'));
    while i < n_meets
        i=i+1;
        day = datetime(year_start, month_start, day_start + (i-1)*freq );
        isskipped = false;
        for j = 1:size(skips,1)
            dayskip = datetime(day.Year, skips(j,1), skips(j,2));
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
        filestr = sprintf('%s - %s - %s,%s\n', dd(s(1)+1:s(2)-1), dd(1:s(1)-1), dd(s(2)+1:end), speaker);
        fwrite(f, filestr);
    end
    fclose(f);
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

for i = 1:n_meets
    day = datetime(year_start, month_start, day_start + (i-1)*freq );
    ind = ords( mod(i, nnames)+1 );
    speaker = names{ind};
    dd = sprintf('%s', day);
    s = strfind(dd, '-');
    fprintf('%s\n', speaker)
%     fprintf('%s - %s - %s    ->    %s\n', dd(s(1)+1:s(2)-1), dd(1:s(1)-1), dd(s(2)+1:end), speaker)
end
