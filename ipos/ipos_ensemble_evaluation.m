clear; close all;
fr = 30;
bin_timeres = .25;
samplesPerBin = fr*bin_timeres;
segs_range = 2.^[0:11];
sample_range = round(samplesPerBin*50*2.^[1:9] + 2*samplesPerBin) ;
nsamples_binned = (sample_range./(fr*bin_timeres))-2;
nsamples_time = nsamples_binned*bin_timeres;
swap_range = round(100*[0, .01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99, 1]);
seg_mod_range = round(100*[0, .01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99, 1]);

isomap_segrange = 2.^[0:8];
randLoops = 1%:10;

stand_prop = swap_range(end-3);
stand_ncells = segs_range(end-3);%nsegs(end-2);
stand_nsamples = sample_range(end-3);
stand_binsamples = nsamples_binned(end-3);
stand_prop_cells = seg_mod_range(end-3);

max_len = max([size(seg_mod_range,2), size(swap_range,2), size(sample_range,2), size(segs_range,2)]);
% topdir = 'D:\Sample Data\ensemble_prob_median';
topdir = 'D:\Sample Data\ensemble_prob_linear';
% model_fitting_test;
% topdir = 'D:\Sample Data\ensemble_prob';
%
props = swap_range;
nrand = length(randLoops);

nantemplate = NaN(max_len, nrand);
temp = [];
temp.ensemMod       = nantemplate;
temp.iposMod        = nantemplate;
temp.spatialDistMod = nantemplate;
temp.IsoMapDist     = nantemplate;
temp.spkDist        = nantemplate;
temp.xval           = NaN(max_len, 1);
ncells      = temp;
nsamples    = temp;
prop_cells  = temp;
prop_samples= temp;
stand_comp  = temp;
isomap_dims = temp;
isomap_dims.iso_ensemMod        = nantemplate;
isomap_dims.iso_iposMod         = nantemplate;
isomap_dims.iso_spatialDistMod  = nantemplate;
temp = [];
%%
nrand = 1;
plotting = false;
ncells = eval_data_conjointipos(ncells, [topdir '\ncells\'], segs_range, nrand, plotting);
nsamples = eval_data_conjointipos(nsamples, [topdir '\nsamples\'], sample_range, nrand, plotting);
prop_cells = eval_data_conjointipos(prop_cells, [topdir '\prop_cells\'], seg_mod_range, nrand, plotting);
prop_samples = eval_data_conjointipos(prop_samples, [topdir '\prop_samples\'], swap_range, nrand, plotting);

isomap_dims = eval_data_conjointipos(isomap_dims, [topdir '\isomap_ncells\'], isomap_segrange, nrand, plotting);
%%
stand_range = {'256_0', '256_1' '256_2'};
stand_comp = eval_data_conjointipos(stand_comp, [topdir '\standard\'], stand_range, 1:3, plotting);



%%

for i = 1:0%length(stand_range)
    fn = sprintf('%s%d_%d.mat', ddir, 256, stand_range(i)-1);
    if isfile(fn)
        %%
        disp(fn)
        temp = load(fn);
        [spks_val, sd_val, ipos_val, iso_val, e_val] = i_e_powermod(temp);
        e = temp.ensem_av1;
        e = e-min(e); 
        e = e./max(e)+ .001;
        thresh_e = median(e)+2*std(e);
%         edges_of_bins  = [-.0015:.001:.0015;
        edges_of_bins  = [.001:.001:1];
        [m] = gb_conjoint_ensemble_model_fitting(e, thresh_e, edges_of_bins, true);
        m.product
%         e_val
        if i<11
            stand_comp.spkDist(i,1)       = spks_val;
            stand_comp.spatialDistMod(i,1)= sd_val;
            stand_comp.iposMod(i,1)       = ipos_val;
            stand_comp.IsoMapDist(i,1)    = iso_val;
            stand_comp.ensemMod(i,1)      = e_val;
        else
            stand_comp.spkDist(i-10,2)       = spks_val;
            stand_comp.spatialDistMod(i-10,2)= sd_val;
            stand_comp.iposMod(i-10,2)       = ipos_val;
            stand_comp.IsoMapDist(i-10,2)    = iso_val;
            stand_comp.ensemMod(i-10,2)      = e_val;
        end
        
    end
end

%%

vars2plot = {'ncells', 'nsamples', 'prop_samples', 'prop_cells'};
dep2plot = {'ensemMod', 'iposMod', 'spatialDistMod', 'spkDist' , 'IsoMapDist'};
clrs2plot = [.2 .2 .2; plasma(6)];
% clrs2plot = [.2 .2 .2; .6 .2 .2; .2 .2 .8; .8 .2 .8];
stand_clr = [.8 1 .8];
medge = 'none';
msize = 10;
figure(102); clf
set(gcf, 'Name', 'Parameter comp', 'Position', [480   301   828   677])

subplot(2,2, 1); cla; hold on
xs = log2(segs_range);
s = rectangle('Position', [log2(stand_ncells)-.2 -1.1 .4 2.2]); s.FaceColor = stand_clr; s.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XTickLabel', segs_range, 'XTickLabelRotation', -90)
axis([-1  length(xs)+1 -1.1 1.1]); axis square
xlabel('Number of cells')

subplot(2,2, 2); cla; hold on
s = rectangle('Position', [log2(stand_nsamples)-.2 -1.1 .4 2.2]); s.FaceColor = stand_clr; s.EdgeColor = 'none';
set(gca, 'XTick', log2(sample_range), 'XTickLabel', nsamples_time, 'XTickLabelRotation', -90)
axis([log2(sample_range(1))/1.1  log2(sample_range(end))*1.1 -1.1 1.1]); axis square
p99 = plot([-1 -1], [0 1], '-', 'Color', stand_clr, 'LineWidth', 10);
xlabel('Total time (sec)')

subplot(2,2, 3); hold on
s = rectangle('Position', [(stand_prop)-2.5 -1.1 5 2.2]); s.FaceColor = stand_clr; s.EdgeColor = 'none';
gt = [0 20:20:80 100];
set(gca, 'XTick', gt)
axis([-10 110 -1.1 1.1]); axis square
xlabel('Percent samples in C_1')

subplot(2,2, 4); cla; hold on
s = rectangle('Position', [(stand_prop_cells)-2.5 -1.1 5 2.2]); s.FaceColor = stand_clr; s.EdgeColor = 'none';
gt = [0 20:20:80 100];
set(gca, 'XTick', gt)
axis([-10 110 -1.1 1.1]); axis square
xlabel('Percent cells modulated')

for varLoop = 1:length(vars2plot)
    subplot(2,2, varLoop); hold on
        x = eval(sprintf('%s.xval', vars2plot{varLoop}));
    if varLoop==1 || varLoop == 2
        x = log2(x);
    else
    end
    for depLoop = 1:length(dep2plot)
        y = eval(sprintf('%s.%s', vars2plot{varLoop}, dep2plot{depLoop}));
        c = clrs2plot(depLoop,:);
        eval(sprintf('p%d = plot(x, nanmean(y, 2), ''Color'', shift_colormap(c, 2), ''LineWidth'', 3);', depLoop));
        for i = 1:length(x)
            xs = x(i)*ones(length(y(i,:)),1);
            scatter(xs, y(i,:), msize, 'MarkerFaceColor', c, 'MarkerEdgeColor', medge, 'MarkerFaceAlpha', .8);
        end
    end
    if varLoop==4
    legend([p1 p2 p3 p4 p5], 'Conj. Ipos', 'Avg. Ipos', 'Spatial d', 'Spk dist',  'IsoMap dist',...
        'Location', 'northwest', 'EdgeColor', 'none')
    end
end
subplot(2,2, 2); 
legend(p99, {'Standard val'}, 'Location', 'northwest', 'EdgeColor', 'none')
subplot(2,2, 1); 
ylabel(sprintf('Sesnsitivity to C_2\n (\\SigmaC_2-\\SigmaC_1)/(\\SigmaC_2+\\SigmaC_1)'));
subplot(2,2, 3); 
ylabel(sprintf('Sesnsitivity to C_2\n (\\SigmaC_2-\\SigmaC_1)/(\\SigmaC_2+\\SigmaC_1)'));

%% Showing the standard data with mixed and chunked switch signal
stan1 = load('D:\Sample Data\ensemble_prob\standard\256_0.mat');
stan2 = load('D:\Sample Data\ensemble_prob\standard\256_11.mat');
% stan1.isoDist = squareform(pdist(stan1.isoMap));
% stan2.isoDist = squareform(pdist(stan2.isoMap));
st = {'stan1' 'stan2'};
[ns, nt] = size(stan1.bin_spks);
figure(101); clf; colormap viridis
set(gcf, 'Name', 'Example spikes, grouped vs interleved', 'Position', [335         228        1288         750])
set(gcf, 'Color', 'w')

nsubs = 6;
clr = plasma(nsubs);

inds = [nt-1000:nt];
inds = [2000:3000];
nt = length(inds);
for i = 1:2
    %
    temp_struct         = eval(sprintf('%s;', st{i}));
    spk_dist            = normalize_matrix(squareform(pdist(temp_struct.bin_spks')));
%     spk_corr            = corr(temp_struct.bin_spks(:,inds), 'Type','Kendall');
%     time_corr           = nanmean(spk_corr,1);
    ensemble_val        = normalize_matrix(temp_struct.ensem1);
    ipos1               = normalize_matrix(nanmean(abs(temp_struct.ipos1), 1));
    spatialDistinct     = normalize_matrix(2.^(temp_struct.sd1));
    isoDist             = normalize_matrix(squareform(pdist(temp_struct.isoMap)));
    %
    sig = temp_struct.switch_binned*1.25 - .15;
    sig = sig(inds);
    spk_dist            = spk_dist(inds);
    ensemble_val        = ensemble_val(inds);
    ipos1               = ipos1(inds);
    spatialDistinct     = spatialDistinct(inds);
    isoDist             = isoDist(inds);
    
    
    subplot_tight(2,2,0+i); cla; hold on
    plot(temp_struct.switch_binned(inds)*120 -10, 'r')
    imagesc(temp_struct.bin_spks(1:100,inds), [0 4])
    axis([ 0 nt -20.5 120.5])
    set(gca, 'YTick', [0:20:100], 'XTick', [0:250:1000], 'XTickLabel', inds(1:250:1001))
    if i==1; ylabel('Example cell #'); end
    if i==2; yyaxis('right'); set(gca, 'Ycolor', 'r'); ylabel('Context', 'Color', 'r'); yticks([]); end

    jj = nsubs*2 + i;
    q = 1;
    subplot_tight(nsubs*2, 2, jj); cla; hold on; jj = jj+2;
    plot(sig,'r')
    plot(ensemble_val,'k')
    ylabel(sprintf('Conjoint\nIpos'))
    axis([ 0 nt -.25 1.1 ])
    set(gca, 'YTick', [], 'XTick', [0:250:1000], 'XTickLabel', '')

    subplot_tight(nsubs*2, 2, jj); cla; hold on; jj = jj+2;
    plot(sig,'r')
    plot(ipos1,'-', 'Color', clr(q,:)); q = q+1;
    ylabel(sprintf('average\nIpos'))
    axis([ 0 nt -.25 1.1 ])
    set(gca, 'YTick', [], 'XTick', [0:250:1000], 'XTickLabel', '')
    
    subplot_tight(nsubs*2, 2, jj); cla; hold on; jj = jj+2;
    plot(sig,'r')
    plot(spatialDistinct,'-', 'Color', clr(q,:)); q = q+1;
    ylabel(sprintf('Spatial\ndistinctive'))
    axis([ 0 nt -.25 1.1 ])
    set(gca, 'YTick', [], 'XTick', [0:250:1000], 'XTickLabel', '')
    
    subplot_tight(nsubs*2, 2, jj); cla; hold on; jj = jj+2;
    plot(sig,'r')
    plot(spk_dist,'-', 'Color', clr(q,:)); q = q+1;
    ylabel(sprintf('Spk\ndist'))
    axis([ 0 nt -.25 1.1 ])
    set(gca, 'YTick', [], 'XTick', [0:250:1000], 'XTickLabel', '')
    
    subplot_tight(nsubs*2, 2, jj); cla; hold on; jj = jj+2;
    plot(sig,'r')
    plot(isoDist,'-', 'Color', clr(q,:)); q = q+1;
    ylabel(sprintf('IsoMap\ndist'))
    axis([ 0 nt -.25 1.1 ])
    set(gca, 'YTick', [], 'XTick', [0:250:1000], 'XTickLabel', inds(1:250:1001))
%     ylabel('Sample')
    xlabel('Sample')
end

%%
figure(104); clf; hold on
set(gcf, 'Name', 'Grouped vs interleaved switch', 'Position', [981   235   531   397])
set(gcf, 'Color', 'w')
i = 0;
b1 = bar(-99, 1, 'FaceColor', [.3 .3 .3], 'FaceAlpha', .7, 'EdgeColor', 'k', ...
    'BarWidth', .6, 'LineStyle', '-', 'LineWidth', 2);
b2 = bar(-98, 1, 'FaceColor', [.3 .3 .3]./3, 'FaceAlpha', .7, 'EdgeColor', 'k', ...
    'BarWidth', .6, 'LineStyle', 'none', 'LineWidth', 2);

gb_quickbar(i-.4, stand_comp.ensemMod(:,1), [.3 .3 .3], '-')
gb_quickbar(i+.4, stand_comp.ensemMod(:,2), [.3 .3 .3]./3, 'none')
i = 2;
gb_quickbar(i-.4, stand_comp.iposMod(:,1), clr(1,:), '-')
gb_quickbar(i+.4, stand_comp.iposMod(:,2), clr(1,:)./1.5, 'none')
i = 4;
gb_quickbar(i-.4, stand_comp.spatialDistMod(:,1), clr(2,:), '-')
gb_quickbar(i+.4, stand_comp.spatialDistMod(:,2), clr(2,:)./1.5, 'none')
i = 6;
gb_quickbar(i-.4, stand_comp.spkDist(:,1), clr(3,:), '-')
gb_quickbar(i+.4, stand_comp.spkDist(:,2), clr(3,:)./1.5, 'none')
i = 8;
gb_quickbar(i-.4, stand_comp.IsoMapDist(:,1), clr(4,:), '-')
gb_quickbar(i+.4, stand_comp.IsoMapDist(:,2), clr(4,:)./1.5, 'none')
plot([-1 9], [0 0], 'k', 'LineWidth', 3)
axis([ -1 9 -1.1 1.1]) 
legend([b1 b2], {'Grouped', 'Interleaved'}, 'EdgeColor', 'none')
set(gca, 'XTick', 0:2:8, 'XTickLabel', {'Conj. Ipos', 'Avg. Ipos', 'Spatial d', 'Spk dist',  'IsoMap dist'})
%%


figure(103); clf

hold on
xs = log2(isomap_segrange);
s = rectangle('Position', [log2(stand_ncells)-.2 -1.1 .4 2.2]); s.FaceColor = stand_clr; s.EdgeColor = 'none';
set(gca, 'XTick', xs, 'XTickLabel', isomap_segrange, 'XTickLabelRotation', -90)
axis([-1 xs(end)+1 -1.1 1.1])
xlabel('Number of IsoMap dims')

vars2plot = {'isomap_dims'};
dep2plot = {'ensemMod', 'iso_ensemMod'};
for varLoop = 1:length(vars2plot)
        x = eval(sprintf('%s.xval', vars2plot{varLoop}));
        x = log2(x);
    for depLoop = 1:length(dep2plot)
        y = eval(sprintf('%s.%s', vars2plot{varLoop}, dep2plot{depLoop}));
        c = clrs2plot(depLoop,:);
        eval(sprintf('p%d = plot(x, nanmean(y, 2), ''Color'', shift_colormap(c, 2), ''LineWidth'', 3)', depLoop))
        for i = 1:length(x)
            xs = x(i)*ones(length(y(i,:)),1);
            scatter(xs, y(i,:), msize, 'MarkerFaceColor', c, 'MarkerEdgeColor', medge, 'MarkerFaceAlpha', .8);
        end
    end
    legend([p1 p2], '256 cells', 'IsoMap', 'Location', 'northwest')
end


%%
function struct_out = eval_data_conjointipos(struct_in, ddir, loop_range, randLoops, plotting)
%     ddir = [topdir '\ncells\'];
    struct_out = struct_in;
    for i = 1:length(loop_range)
        for j = randLoops
            %
            if iscell(loop_range(i)) % only used in standard case
                fn = sprintf('%s%s.mat', ddir, loop_range{i});
%                 j = i;
                struct_out.xval(i)       = NaN;
            else
                fn = sprintf('%s%d_%d.mat', ddir, loop_range(i), j-1);
                struct_out.xval(i)       = loop_range(i);
            end
            if isfile(fn)
                disp(fn)
                temp = load(fn);
                [spks_val, sd_val, ipos_val, iso_val, e_val] = i_e_powermod(temp);

                struct_out.spkDist(i,j)       = spks_val;
                struct_out.spatialDistMod(i,j)= sd_val;
                struct_out.iposMod(i,j)       = ipos_val;
                struct_out.IsoMapDist(i,j)    = iso_val;
                struct_out.ensemMod(i,j)      = e_val;  
                if plotting % j==Inf
                    figure(i+2000); clf;
                    subplot(1,3,1);
                    hold on;
                    inds = temp.switch_binned==0;
                    scatter(temp.isoMap(inds,1), temp.isoMap(inds,2), 'k.')
                    inds = temp.switch_binned==1;
                    scatter(temp.isoMap(inds,1), temp.isoMap(inds,2), 'r.')


                    subplot(1,3,2);
                    hold on;
                    inds = temp.switch_binned==0;
                    scatter(temp.isoMap(inds,3), temp.isoMap(inds,2), 'k.')
                    inds = temp.switch_binned==1;
                    scatter(temp.isoMap(inds,3), temp.isoMap(inds,2), 'r.')

                    subplot(1,3,3);
                    hold on;
                    inds = temp.switch_binned==0;
                    scatter(temp.isoMap(inds,1), temp.isoMap(inds,3), 'k.')
                    inds = temp.switch_binned==1;
                    scatter(temp.isoMap(inds,1), temp.isoMap(inds,3), 'r.')
                end
            else
                warning('Could not find file: %s', fn)
            end
        end
    end
end

function [spks_val, sd_val, ipos_val, iso_val, e_val] = i_e_powermod(temp_struct)
inds2 = temp_struct.switch_binned==1; % temp_struct.switch_binned==0 is first ref frame

% ensemble_val = abs(temp_struct.ensem1);
ensemble_val = abs(temp_struct.ensem_av1);
ipos1 = nanmean(abs(temp_struct.ipos1), 1);
spatialDistinct = 2.^(temp_struct.sd1);% - abs(temp_struct.sd2); % sd1 and sd2 are essentially the same
isoDist = squareform(pdist(temp_struct.isoMap));
% calculate selectivity for each (a-b)/(a+b)
if isfield(temp_struct, 'bin_spks')
    spk_dist = squareform(pdist(temp_struct.bin_spks'));
    spks_val    = (nansum(spk_dist(inds2)) - nansum(spk_dist(~inds2))) /...
        ( nansum(spk_dist(inds2)) + nansum(spk_dist(~inds2)) );
    
else
    spks_val = NaN;
end
% sd_val      = (nansum(spatialDistinct(inds2)) - nansum(spatialDistinct(~inds2))) /...
%     ( nansum(spatialDistinct(inds2)) + nansum(spatialDistinct(~inds2)) );
% ipos_val    = (nansum(ipos1(inds2)) - nansum(ipos1(~inds2))) /...
%     ( nansum(ipos1(inds2)) + nansum(ipos1(~inds2)) );
% iso_val    = (nansum(isoDist(inds2)) - nansum(isoDist(~inds2))) /...
%     ( nansum(isoDist(inds2)) + nansum(isoDist(~inds2)) );
% e_val       = (nansum(ensemble_val(inds2)) - nansum(ensemble_val(~inds2))) /...
%     ( nansum(ensemble_val(inds2)) + nansum(ensemble_val(~inds2)) );

sd_val      = (nanmean(spatialDistinct(inds2)) - nanmean(spatialDistinct(~inds2))) /...
    ( nanmean(spatialDistinct(inds2)) + nanmean(spatialDistinct(~inds2)) );
ipos_val    = (nanmean(ipos1(inds2)) - nanmean(ipos1(~inds2))) /...
    ( nanmean(ipos1(inds2)) + nanmean(ipos1(~inds2)) );
iso_val    = (nanmean(isoDist(inds2)) - nanmean(isoDist(~inds2))) /...
    ( nanmean(isoDist(inds2)) + nanmean(isoDist(~inds2)) );
e_val       = (nanmean(ensemble_val(inds2)) - nanmean(ensemble_val(~inds2))) /...
    ( nanmean(ensemble_val(inds2)) + nanmean(ensemble_val(~inds2)) );

end
            
function gb_quickbar(x, y, clr, lineedge)
xs = gb_rand_jitter(y,3);
bar(x, nanmean(y), 'FaceColor', clr, 'FaceAlpha', .7, 'EdgeColor', 'k', ...
    'BarWidth', .6, 'LineStyle', lineedge, 'LineWidth', 2)
scatter(x+y*0+xs, y, 50, 'o', 'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5)
errorbar(x, nanmean(y), std(nanmean(y)), 'Color', 'k')
end