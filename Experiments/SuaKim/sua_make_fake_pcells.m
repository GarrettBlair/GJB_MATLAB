clear;

temp = load("D:\GarrettBlair\APA\HPCACC34990\processed_files\2025_04_09_H12_32_43_IL5_@placecells_HPC_miniscope1.mat");
t = temp.ms.timestamps;
x = temp.ms.arena.x;
y = temp.ms.arena.y;

%%
load('C:\Users\gjb326\Documents\mentorship\Sua Kim\data\example_position_34990_IL5.mat');
binsx = [-45:10:45];
binsy = binsx;
%%
rng('default')
speed = sqrt((x(1:end-1) - x(2:end)).^2 + (y(1:end-1) - y(2:end)).^2);
speed = [speed(1); speed];
spd = conv(speed, ones(30,1)/30, 'same');
speed(16:end-15) = spd(16:end-15);
dt = abs(diff(t));
dt = [dt(1); dt];

[occ_map, binx, biny] = binned_statistic2d(x, y, dt,  binsx, binsy, 'sum');
[spd_map, ~, ~] = binned_statistic2d(x, y, speed, binsx, binsy, 'median');

ux = unique(binx);
uy = unique(biny);

nsegs = 100; nt = length(x);

[xp] = randi(max([ux;uy]), [nsegs,2]);
randdrive = randn(nsegs, nt)/3 +1.5;
spks = zeros(nsegs, nt);
spks2 = zeros(nsegs, nt);
refractory = [0*ones(30,1); 1; -1000*ones(30,1)];
for ii = 1:nsegs
    prefdist = sqrt((binx-xp(ii,1)).^2 + (biny-xp(ii,2)).^2);
    prefdist = prefdist - min(prefdist);
    prefdist = 1 - prefdist/max(prefdist);
%     spk = ((speed+1).*prefdist)' > randdrive(ii,:);
    spk = (prefdist)' > randdrive(ii,:);
%     figure; hold on
%     plot(x,y);
%     plot(binsx(binx(s)), binsy(biny(s)), 'ro')
%     plot(binsx(xp(ii,1)), binsy(xp(ii,2)), 'mx')
%     spks2(ii,:) = double(spk>0);
    spk = conv(spk, refractory, 'same');
    spks(ii,:) = spk>0;
    
end
smap = zeros(length(binsy)-1, length(binsx)-1, nsegs);
pmap = zeros(length(binsy)-1, length(binsx)-1, nsegs);
% figure(1); clf
% for ii = 1:nsegs
%     [sm, ~, ~] = binned_statistic2d(x, y, spks(ii,:), binsx, binsy, 'sum');
%     sm(isnan(sm))=0;
%     smap(:,:,ii) = sm;
%     pmap(:,:,ii) = squeeze(smap(:,:,ii))./occ_map;
%     if ii<=50
%         subplot(5,10,ii)
%     imagesc(squeeze(pmap(:,:,ii)));
%     axis image off
%     end
% end

%%
reward_time_ms = temp.ms.room.entranceTimes;
save('C:\Users\gjb326\Documents\mentorship\Sua Kim\data\example_position_34990_IL5.mat', 'x', 'y', 't', 'spks', 'reward_time_ms');

