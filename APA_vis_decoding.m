ddir = 'Z:\f\fentonlab\RAWDATA\CaImage\GarrettBlair\ImagingData\';
h = load([ddir 'APA_HPCACC\HPCACC24500\processed_files\2023_06_20_H13_11_53_TR15_@placecells_HPC_miniscope1.mat']);
a = load([ddir 'APA_HPCACC\HPCACC24500\processed_files\2023_06_20_H13_11_53_TR15_@placecells_ACC_miniscope2.mat']);


%%
clearvars -except a h ddir

h_spks = h.ms.spks;
a_spks = a.ms.spks;
t = average_spks_time(h.ms.timestamps'./1000, .5, h.ms.timestamps./1000, false, 'mean');
roomx = average_spks_time(h.ms.room.x', .5, h.ms.timestamps./1000, false, 'mean');
roomy = average_spks_time(h.ms.room.y', .5, h.ms.timestamps./1000, false, 'mean');
% roomspd = average_spks_time(h.ms.room.speed', .5, h.ms.timestamps./1000, false, 'mean');
arenax = average_spks_time(h.ms.arena.x', .5, h.ms.timestamps./1000, false, 'mean');
arenay = average_spks_time(h.ms.arena.y', .5, h.ms.timestamps./1000, false, 'mean');
arenaspd = average_spks_time(h.ms.arena.speed', .5, h.ms.timestamps./1000, false, 'median');
hpc_spks = average_spks_time(h_spks, .5, h.ms.timestamps./1000, false, 'sum');
acc_spks = average_spks_time(a_spks, .5, a.ms.timestamps./1000, false, 'sum');

vars2save = {'hpc_spks', 'acc_spks', 't', 'roomx', 'roomy', 'arenax', 'arenay', 'arenaspd'};

save('D:\Sample Data\conjoint_ipos_ca_imaging\HPCACC24500_TR15.mat', vars2save{:})
%%
hgood = (h.ms.room.pcell_stats.infoProb)<=1;%.05;
agood = (a.ms.room.pcell_stats.infoProb)<=1;%.05;
pfh = reshape(h.ms.room.pfields_smooth(hgood,:,:), [sum(hgood), 20*20]);
pfa = reshape(a.ms.room.pfields_smooth(agood,:,:), [sum(agood), 20*20]);
hpc_spks_std = std_scale_rows(hpc_spks);
acc_spks_std = std_scale_rows(acc_spks);
pp = [pfh;pfa]'*[hpc_spks_std(hgood,:); acc_spks_std(agood,:)];
nf = size(hpc_spks,2);
pp2 = reshape(normalize_cols(pp), [20,20,nf]);
%%
figure(99);clf
binx = normalize_matrix(roomx)*40;
biny = normalize_matrix(roomy)*40;
for i=600:nf
    %%
    p = squeeze(pp2(:,:,i));
    p(p<.9) = .25;
    p = imresize(p, 2,'Method', 'bilinear');
    imagesc(p, [0 1]); hold on
    colormap magma
%     plot(binx(i), biny(i), 'ro');
    scatter(binx(i), biny(i), 60, 'ro', 'MarkerFaceColor', 'r');
    axis tight
    hold off
    drawnow
    pause(.15)
end