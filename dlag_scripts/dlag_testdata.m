%%
% load('C:\Users\gjb326\Documents\MATLAB\Darryl_simdata.mat');
h = load('D:\GarrettBlair\APA\HPCACC24500\simple_files\2023_06_20_H13_11_53_TR15_@placecells_HPC_miniscope1.mat');
a = load('D:\GarrettBlair\APA\HPCACC24500\simple_files\2023_06_20_H13_11_53_TR15_@placecells_ACC_miniscope2.mat');
%
[hns, hnt] = size(h.spks);
[aans, ant] = size(a.spks);
ns = min([hns, aans]);

binWidth = 500;

rng('default')
hspks = h.spks(randperm(hns),:)>0;
aspks = a.spks(randperm(aans),:)>0;
hspks = hspks(1:ns,:);
aspks = aspks(1:ns,:);

ht = h.time_ms - min(h.time_ms(1), a.time_ms(1));
at = a.time_ms - min(h.time_ms(1), a.time_ms(1));

[hspks, ~] = bin_spks_time(hspks, binWidth, ht,  false);
[aspks, ~] = bin_spks_time(aspks, binWidth, at,  false);
min_samples = min(size(hspks,2), size(aspks,2));
%% Split data into entries towards shock zone
spks = cat(1, hspks(:, 1:min_samples), aspks(:, 1:min_samples));
seqTrue = [];
i = 1;
seqTrue(i).trialId = i;
seqTrue(i).T = size(spks,2);
% seqTrue(i).xsm = [];
spks(isnan(spks)) = 0;
seqTrue(i).y = spks;

seqTrue = cutTrials(seqTrue);
seqTrue.y = rand(size(spks));
% seqTrue(i).inds = inds(validinds);
%% Split data into entries towards shock zone

if false
    ht = h.time_ms'./1000;
    at = a.time_ms'./1000;

    win = 300;
    ne = length(h.shock_approach_ind);

    clearvars seqTrue;
    % seqTrue = struct('T', win*2+1, 'trialId', NaN, 'y', haspks, 'xsm', []);
    seqTrue = struct();
    bin_ds = 10;
    dt = (median(abs(diff(ht))));
    bin_size = 1000*dt*bin_ds; % milliseconds
    trial_length = floor(  (win*2+1)/bin_ds  );
    for i = 1:ne
        inds = h.shock_approach_ind(i)-win:h.shock_approach_ind(i)+win;
        validinds = inds>0 & inds<hnt;
        hs = hspks(:, inds(validinds));
    %     [hsb, gt] = bin_spks_time(hs, timeres, ht2, true);


        inds = a.shock_approach_ind(i)-win:a.shock_approach_ind(i)+win;
        validinds = inds>0 & inds<ant;
        as = aspks(:, inds(validinds));
    %     [asb, gt] = bin_spks_time(as, timeres, at2, true);


        spks = cat(1, hs, as);
        spks = bin_spks(spks, bin_ds, false);
        if size(spks,2) == trial_length
            seqTrue(i).trialId = i;
            seqTrue(i).T = size(spks,2);
            seqTrue(i).xsm = [];
            seqTrue(i).y = spks;
            seqTrue(i).inds = inds(validinds);
        end
    end
end
%%
% hspks = bin_spks_time(hspks(1:ns,:), 1, ht, false);
% aspks = bin_spks_time(aspks(1:ns,:), 1, at, false);
% haspks = cat(1,hspks,aspks);
% halabel = cat(1,hspks*0,aspks*0+1);


% seqTrue = struct('T', size(bs,2), 'trialId', 1, 'y', cat(1, bs, bs2), 'xsm', []);
% seqTrue = struct('T', size(haspks,2), 'trialId', 1, 'y', haspks, 'xsm', []);

% seqTrue(1).T = size(haspks,2);
% seqTrue(1).trialId = 1;
% seqTrue(1).y = haspks;
% seqTrue(1).xsm = [];









