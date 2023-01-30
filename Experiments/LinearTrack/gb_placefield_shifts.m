%%


topdir = 'D:\Sample Data\doi_10.5068_D1ZT2S__v1\Blair_et_al._2022_repo\tadblair-main\Blair_et_al_DATA\'; % 'D:\MiniScopeData\Blair_et_al\'; %top level folder where data files are found
shockdir=[topdir 'shocktimes\']; %folder where the goodies are stored
cellmapdir=[topdir 'cellmaps\']; %folder where the goodies are stored
predir=[topdir 'pretrn\'];
postdir=[topdir 'prepost\'];
% datadir=[topdir 'pretrn\'];

% cellmat_shockcols = [...
%     %    3 6 9;...  %Hipp6 3 day
%     5 6 8;...  %Hipp6 ext
%     6 7 10;... %Hipp7
%     3 4 7;...  %Hipp8
%     4 5 8;...  %Hipp9
%     8 9 12;... %Hipp15
%     5 6 9;...  %Hipp31
%     
%     5 6 9;... %Hipp12 (barrier)
%     5 6 9;... %Hipp13 (barrier)
%     6 7 8;... %Hipp18 (barrier) ext
%     7 8 11;... %Hipp35 (barrier)
%     5 6 7;... %Hipp30 (barrier) ext
%     7 8 11;... %Hipp34 (barrier)
%     
%     9 10 13;... %Hipp12 (shock)
%     9 10 13;... %Hipp13 (shock)
%     12 13 14;... %Hipp18 (shock) ext
%     11 12 15;...   %Hipp35 (shock)
%     
%     9 10 13;... %Hipp30 (shock+scop)
%     10 11 14;... %Hipp32 (shock+scop)
%     11 12 15;... %Hipp34 (shock+scop)
%     11 12 15;...   %Hipp36 (shock+scop)
%     
%     22 23 26;... %Hipp12 (shock+scop)
%     22 23 26;... %Hipp13 (shock+scop)
%     22 23 26;... %Hipp15 (shock+scop)
%     22 23 26;...   %Hipp18 (shock+scop)
%     19 20 21;...   %Hipp31 (shock+scop)
%     
%     13 14 15;... %Hipp30 (drug free shock2)
%     19 20 24;... %Hipp32 (drug free shock2)
%     
%     14 15 16;... %Hipp31 (scopolamine alone)
%     5 6 7;... %Hipp32 (scopolamine alone)
%     6 7 8]; %Hipp36 (scopolamine alone)
% sesslabel = [1 1 1 1 1 1 3 3 3 3 3 3 1 1 1 1 0 0 0 0 0 0 0 0 0 2 2 4 4 4];
% % 0==scopo+shock, 1==shock, 2==2nd shock, 3==barrier, 4==scopo alone
% cellmat_Anum = [6 7 8 9 15 31 12 13 18 35 30 34 12 13 18 35 30 32 34 36 12 13 15 18 31 30 32 31 32 36];
% is_shock  = [1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0];

% from fig4A_analysis.m
cellmat_shockcols = [...
    %    3 6 9;...  %Hipp6 3 day
    5 6 8;...  %r=1 Hipp6 ext
    6 7 11;... %r=2 Hipp7
    3 4 7;...  %r=3 Hipp8
    4 5 8;...  %r=4 Hipp9
    8 9 12;... %r=5 Hipp15
    5 6 9;...  %r=6 Hipp31

    5 6 9;... %r=7 Hipp12 (barrier)
    5 6 9;... %r=8 Hipp13 (barrier)
    6 7 9;... %r=9 Hipp18 (barrier) ext
    7 8 11;... %r=10 Hipp35 (barrier)
    5 6 9;... %r=11 Hipp30 (barrier) ext
    7 8 11;... %r=12 Hipp34 (barrier)

    9 10 13;... %r=13 Hipp12 (shock)
    9 10 13;... %r=14 Hipp13 (shock)
    12 13 16;... %r=15 Hipp18 (shock) ext
    11 12 15;...   %r=16 Hipp35 (shock)
    
    9 10 13;... %r=17 Hipp30 (shock+scop)
    10 11 14;... %r=18 Hipp32 (shock+scop)
    11 12 15;... %r=19 Hipp34 (shock+scop)
    11 12 15;...   %r=20 Hipp36 (shock+scop)
    
    22 23 26;... %r=21 Hipp12 (shock+scop)
    22 23 26;... %r=22 Hipp13 (shock+scop)
    22 22 26;... %r=23 Hipp15 (shock+scop)    %%%% 23 -->22
    22 23 26;...   %r=24 Hipp18 (shock+scop)
    19 20 23;...   %r=25 Hipp31 (shock+scop)
    
    13 14 17;... %r=26 Hipp30 (drug free shock2)
    19 20 21;... %r=27 Hipp32 (drug free shock2)
    
    14 15 16;... %r=28 Hipp31 (scopolamine alone)
    5 6 7;... %r=29 Hipp32 (scopolamine alone)
    6 7 8]; %r=30 Hipp36 (scopolamine alone)

sesslabel = [1 1 1 1 1 1 3 3 3 3 3 3 1 1 1 1 0 0 0 0 0 0 0 0 0 2 2 4 4 4];
% 0==scopo+shock, 1==shock, 2==2nd shock, 3==barrier, 4==scopo alone
is_shock  = [1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0];
cellmat_Anum = [6 7 8 9 15 31 12 13 18 35 30 34 12 13 18 35 30 32 34 36 12 13 15 18 31 30 32 31 32 36];

rats_to_analyze=[1:4 6 13:17 20:27]; %dfdf rats
% rats_to_analyze=[1:4 6:27];
cellmat_shockcols = cellmat_shockcols(rats_to_analyze, :);
sesslabel = sesslabel(rats_to_analyze);
is_shock = is_shock(rats_to_analyze);
cellmat_Anum = cellmat_Anum(rats_to_analyze);

% load([shockdir 'shocktimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
% load([shockdir 'shocktimes2']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
% load([shockdir 'scopshocktimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
% load([shockdir 'barriertimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
shocktime_anum = [6 7 8 9 12 13 15 18 30 31 32 34 35 36 37];


%% load matching data
PRE_TRN_comparison = false;
shockcenter = 12;
peakbins = [-23.5:1:23.5];
posbins = [1:23];
% ratebins = [-1.025:.05:1.025];
ratebins = [-.0105:.001:.0105];
shockdistbins = [-11.5:1:11.5];

peak_binshift = NaN(length(cellmat_Anum), length(peakbins)-1);
ratechange_by_bin = NaN(length(cellmat_Anum), length(posbins)-1);
peak_rateshift = NaN(length(cellmat_Anum), length(ratebins)-1);
peak_shockdistshift = NaN(length(cellmat_Anum), length(shockdistbins)-1);

overlap_pre = NaN(1, length(cellmat_Anum)); overlap_post = NaN(1, length(cellmat_Anum));
deltaT = 1; %second integration for spike binning
for aLoop = 1:length(cellmat_Anum)
    %%
    clearvars sessionNums frame*
    ANIMALNUM = cellmat_Anum(aLoop);%6;
    s1 = cellmat_shockcols(aLoop, 1);
    s2 = cellmat_shockcols(aLoop, 2);
    s3 = cellmat_shockcols(aLoop, 3);
    switch sesslabel(aLoop)
        case 0
            load([cellmapdir '\Hipp' num2str(ANIMALNUM) '_scopshock_cmap.mat'])
        case 1
            load([cellmapdir '\Hipp' num2str(ANIMALNUM) '_shock_cmap.mat'])
        case 2
            load([cellmapdir '\Hipp' num2str(ANIMALNUM) '_shock2_cmap.mat'])
        case 3
            load([cellmapdir '\Hipp' num2str(ANIMALNUM) '_barrier_cmap.mat'])
        case 4
            load([cellmapdir '\Hipp' num2str(ANIMALNUM) '_scopo_cmap.mat'])
    end
    % load session data
    f1 = [predir 'Hipp' num2str(ANIMALNUM) '_linear' num2str(s1) '_predata.mat'];
    f2 = [predir 'Hipp' num2str(ANIMALNUM) '_linear' num2str(s2) '_trndata.mat'];
    f3 = [postdir 'Hipp' num2str(ANIMALNUM) '_linear' num2str(s1) '_predata.mat'];
    f4 = [postdir 'Hipp' num2str(ANIMALNUM) '_linear' num2str(s3) '_postdata.mat'];
    if all(isfile({f1, f2, f3, f4}))
        %
    load(f1)
    load(f2)
    load(f3)
    load(f4)
    sess1 = (sessionNums==s1);
    sess2 = (sessionNums==s2);
    sess3 = (sessionNums==s3);
    if PRE_TRN_comparison==true
    frame_pre1 = eval(sprintf('predata'));
    frame_shock = eval(sprintf('trndata'));
    match = cmap(:, [sess1|sess2]);
    else % evaluate PRE-POST effect
    frame_pre1 = eval(sprintf('trndata'));
    frame_shock = eval(sprintf('postdata'));
    match = cmap(:, [sess2|sess3]);
    end
shared = sum(match>0,2)==2;
segs1 = match(shared, 1);
segs2 = match(shared, 2);
a = [frame_pre1.dcurve_LR(segs1,:) frame_pre1.dcurve_RL(segs1,:)];
b = [frame_shock.dcurve_LR(segs2,:) frame_shock.dcurve_RL(segs2,:)];
%  smoothing place fields to find peak?
s = ones(5,1); s = s./sum(s); %s = ones(sum(shared), 1)*s';
a1 = conv2(1, s, frame_pre1.dcurve_LR(segs1,:), 'same');
a2 = conv2(1, s, frame_pre1.dcurve_RL(segs1,:), 'same');
a = [a1 a2];
b1 = conv2(1, s, frame_shock.dcurve_LR(segs2,:), 'same');
b2 = conv2(1, s, frame_shock.dcurve_RL(segs2,:), 'same');
b = [b1 b2];
% a = normalize_rows(a); b = normalize_rows(b);
[peakval1, maxbin1] = max(a,[],2);
[peakval2, maxbin2] = max(b,[],2);
[nsegs, nbins] = size(frame_pre1.dcurve_LR(segs1,:));
peakval2 = peakval2*0;
for j = 1:nsegs
peakval2(j) = b(j, maxbin1(j));
end
[~, ord] = sort(maxbin1);
a = a(ord,:);
b = b(ord,:);

maxbin1(maxbin1>nbins) = maxbin1(maxbin1>nbins)-nbins;
maxbin2(maxbin2>nbins) = maxbin2(maxbin2>nbins)-nbins;
for j = 1:nbins
    thisbin = j==maxbin1;
    ratechange_by_bin(aLoop, j) = nanmean((peakval1(thisbin)-peakval2(thisbin)));
end
shockdist1 = abs(maxbin1-shockcenter); shockdist2 = abs(maxbin2 - shockcenter);
% maxbin1 = mod(maxbin1, nbins);
[peak_binshift(aLoop,:), ~, bin_x] = histcounts((maxbin1-maxbin2), peakbins, 'Normalization', 'count');
prc = (peakval1-peakval2); %(peakval1-peakval2)./max([peakval1], [], 2);
[peak_rateshift(aLoop,:)] = histcounts(prc, ratebins, 'Normalization', 'count');
[peak_shockdistshift(aLoop,:)] = histcounts((shockdist1-shockdist2), shockdistbins, 'Normalization', 'count');
    
    else
        if ~isfile(f1); fprintf('%s\n', f1); end
        if ~isfile(f2); fprintf('%s\n', f2); end
        if ~isfile(f3); fprintf('%s\n', f3); end
    end
end
%
figure(93); clf
%
% peak_binshift = normalize_rows(peak_binshift);
% peak_rateshift = normalize_rows(peak_rateshift);
% peak_shockdistshift = normalize_rows(peak_shockdistshift);
% % % peak_shockdistshift = ratechange_by_bin;
for i = 1:length(cellmat_Anum)
subplot(141); hold on
plot(peakbins(1:end-1)+mean(abs(diff(peakbins))), peak_binshift(i,:))
subplot(142); hold on
plot(ratebins(1:end-1)+mean(abs(diff(ratebins))), peak_rateshift(i,:))
subplot(143); hold on
plot(posbins+mean(abs(diff(posbins))), ratechange_by_bin(i,:))
subplot(144); hold on
plot(shockdistbins(1:end-1)+mean(abs(diff(shockdistbins))), peak_shockdistshift(i,:))


end
%
figure(95); clf
scop = sesslabel==0; % scoposhock
shk = sesslabel==1; % shock
binsize_cm = 250/23;

a1 = nanmean(peak_binshift(scop,:),1);
b1 = nanmean(peak_rateshift(scop,:),1);
c1 = nanmean(ratechange_by_bin(scop,:),1);
d1 = nanmean(peak_shockdistshift(scop,:),1);
c1std = nanstd(ratechange_by_bin(shk,:),1);
a1 = a1./nansum(a1); b1 = b1./nansum(b1); d1 = d1./nansum(d1);

subplot(141); hold on
x = peakbins(1:end-1)+mean(abs(diff(peakbins)))/2;
plot(binsize_cm*x, (a1), 'm')
subplot(142); hold on
x = ratebins(1:end-1)+mean(abs(diff(ratebins)))/2;
plot(x, (b1), 'm')
subplot(143); hold on
x = posbins+mean(abs(diff(posbins)))/2;
plot(binsize_cm*x, c1, 'm')
shadedErrorBar(binsize_cm*x, c1, c1std, 'lineprops', '-m')
subplot(144); hold on
x = shockdistbins(1:end-1)+mean(abs(diff(shockdistbins)))/2;
plot(binsize_cm*x, (d1), 'm')


a2 = nanmean(peak_binshift(shk,:),1);
b2 = nanmean(peak_rateshift(shk,:),1);
c2 = nanmean(ratechange_by_bin(shk,:),1);
c2std = nanstd(ratechange_by_bin(shk,:),1);
d2 = nanmean(peak_shockdistshift(shk,:),1);
a2 = a2./nansum(a2); b2 = b2./nansum(b2); d2 = d2./nansum(d2);

subplot(141); hold on
x = peakbins(1:end-1)+mean(abs(diff(peakbins)))/2;
plot(binsize_cm*x, (a2), 'k')
subplot(142); hold on
x = ratebins(1:end-1)+mean(abs(diff(ratebins)))/2;
plot(x, (b2), 'k')
subplot(143); hold on
x = posbins+mean(abs(diff(posbins)))/2;
plot(binsize_cm*x, c2, 'k')
shadedErrorBar(binsize_cm*x, c2, c2std, 'lineprops', '-k')
subplot(144); hold on
x = shockdistbins(1:end-1)+mean(abs(diff(shockdistbins)))/2;
plot(binsize_cm*x, (d2), 'k')
















