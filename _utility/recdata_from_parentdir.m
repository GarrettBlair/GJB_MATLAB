function [sessDate, trialname, trial_type, trial_num] = recdata_from_parentdir(ms_parentDir)
%% GET THE REC TIME
s       = ms_parentDir;
% directory is expected to be similar to
% directorypath/year_month_day/hour_min_sec_TR10/
% first find the year (>=2020)
ystart  = strfind(s, '202');
% forward slashes
fsl     = strfind(s, '/');
fsl     = fsl(fsl > ystart(1));
% fist find the year (>=2020)
us      = strfind(s, '_');
us      = us(us   > ystart(1));

yy = s(ystart:ystart+3);
MM = s(us(1)+1:us(1)+2);
dd = s(us(2)+1:us(2)+2);

us = us (us>fsl(1));
hh = s(fsl(1)+1:fsl(1)+2);

mm = s(us(1)+1:us(1)+2);
ss = s(us(2)+1:us(2)+2);

datevec = sprintf('%s_%s_%s_%s_%s_%s', yy, MM, dd, hh, mm, ss);
ms_timefmt = ['yy_MM_dd_HH_mm_ss'];
sessDate = datetime(datevec, 'InputFormat', ms_timefmt);


trialname = s(us(3)+1:fsl(end)-1);

switch trialname(1:2)
    case 'WT' % WATER MANIP
        trial_type = 'WTR';
        trial_num = str2double( trialname( strfind(trialname, 'WTR')+3:end) );
    case 'TR' % TRAINING
        trial_type = 'TR';
        trial_num = str2double( trialname( strfind(trialname, 'TR')+2:end) );
    case 'CO' % CONFLICT
        trial_type = 'CON';
        trial_num = str2double( trialname( strfind(trialname, 'CON')+3:end) );
    case 'DR' % DARK MANIP
        trial_type = 'DRK';
        trial_num = str2double( trialname( strfind(trialname, 'DRK')+3:end) );
    case 'NE' % NEW TRAIN / CONFLICT
        trial_type = 'CON';
        trial_num = str2double( trialname( strfind(trialname, 'NEW')+3:end) );
    case 'RE' % RETENTION
        trial_type = 'RET';
        trial_num = str2double( trialname( strfind(trialname, 'RET')+3:end) );
    case 'HA' % HABITUATION
        trial_type = 'HAB';
        trial_num = str2double( trialname( strfind(trialname, 'HAB')+3:end) );
    case 'HC' % HABITUATION
        trial_type = 'HC';
        trial_num = str2double( trialname( strfind(trialname, 'HC')+2:end) );
    otherwise
        trial_type = 'IDK';
        trial_num = -1;
end
