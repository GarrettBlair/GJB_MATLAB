%%
clc
n_meets = 7;
freq = 14; % days between meetings
names = { 'Simon' 'Garrett' 'Claudia' 'Darryl' 'Jenessa'};
nnames = length(names);
ords = 1:nnames; % randperm(length(names));
for i = 1:n_meets
    day = datetime(2023, 9, 13 + (i-1)*freq );
    ind = ords( mod(i, nnames)+1 );
    speaker = names{ind};
    fprintf('%s    -    %s\n', day, speaker)
end
