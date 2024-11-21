load('D:\Sample Data\2023_06_20_H12_34_14_TR15_@placecells_HPC_miniscope1.mat')

% binraw = average_spks_time(craw, .2, time_ms./1000, false, 'median');
% binspk = average_spks_time(spks, .2, time_ms./1000, false, 'median');
% 
% figure; 
% stacked_traces(normalize_rows(binraw), .8)
% close all
% data=[]; data.y=binraw; RunParam=[]; RunParam.N_LATENT = 3;
% [EstParam, Result] = cilds(data, RunParam);
%%
d = average_spks_time(shock_zone_dist', .2, time_ms./1000, false, 'median');

figure; stacked_traces(normalize_rows(Result.z));

hold on; plot(normalize_rows(d)+1, 'r')