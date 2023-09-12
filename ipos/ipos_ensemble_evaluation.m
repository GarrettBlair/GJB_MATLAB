fr = 30;
bin_timeres = .25;
samplesPerBin = fr*bin_timeres;
segs_range = 2.^[0:11];
sample_range = round(samplesPerBin*50*2.^[0:7] + 2*samplesPerBin) ;
swap_range = round(100*[0, .01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99, 1]);
seg_mod_range = round(100*[0, .01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99, 1]);
randLoops = 1:10;
topdir = 'D:\Sample Data\ensemble_ipos_cosine';
%%
ddir = [topdir '\prop_mod\'];
props = swap_range;
props_min = abs(swap_range - 50) + 50;
nrand = length(randLoops);
props_e_mod = NaN(length(props),nrand);
props_i_mod = NaN(length(props),nrand);
props_sd_mod = NaN(length(props),nrand);
props_iso_mod = NaN(length(props),nrand);
props_e_signal = NaN(length(props),nrand);
props_i_signal = NaN(length(props),nrand);
for i = 1:length(props)
%                 figure(i); clf;
    for j = randLoops
    %%
        fn = sprintf('%sipos_prop%d_%d.mat', ddir, props(i), j-1);
        if isfile(fn)
            temp = load(fn);
            [i_val, e_val, sd_val, iso_val, ns] = i_e_powermod(temp);
            props_e_mod(i,j) = e_val; %nansum(abs(e(s==smin)))/nansum(abs(e));
            props_i_mod(i,j) = i_val;%nansum(abs(ipos(s==smin)))/nansum(abs(ipos));
            props_sd_mod(i,j) = sd_val;%nansum(abs(ipos(s==smin)))/nansum(abs(ipos));
            props_iso_mod(i,j) = iso_val;%nansum(abs(ipos(s==smin)))/nansum(abs(ipos));
            e1 = temp.ensem1;
            e2 = temp.ensem2;
            sd1 = temp.sd1;
            sd2 = temp.sd2;
            e = (e1-e2);
            i1 = nanmean(abs(temp.ipos1), 1);
            i2 = nanmean(abs(temp.ipos2), 1);
            ipos = nanmean(abs(temp.ipos1) - abs(temp.ipos2), 1);
            s = temp.switch_binned;
            [so, ord] = sort(s, 'ascend');

            if j==1 && i==6
            sd = 2.^(temp.sd2);
            figure(100+i); clf;
            subplot(4,2,3)
            plot(e, 'b')
            subplot(4,2,5)
            plot(ipos, 'k')
            subplot(4,2,7)
            plot(sd, 'r')
            subplot(4,2,1)
            imagesc(temp.bin_spks1, [0 2])
            subplot(122)
            scatter(ipos,sd)
            [h,p] = corr(ipos',sd');
            
            axis equal
            drawnow
            end
            
            if j==1
                d = temp.isoMap;
                d = normalize_cols(d);
                figure(i+1000); clf; 
                subplot(4,1,1)
                plot(temp.ensem1-temp.ensem2, 'b')
                title(props(i))
                subplot(4,1,2:4)
                hold on
                plot(d(:,1), 'k')
                plot(d(:,2)+1, 'k')
                plot(d(:,3)+2, 'k')
            end
        end
    end
end
drawnow

%%
ddir = [topdir '\shuffling\'];
isodims = [1 3 5 10 30 60];
nrand = 1;%length(randLoops);
shuff_e_mod = NaN(length(isodims),nrand);
shuff_i_mod = NaN(length(isodims),nrand);
figure(5001); clf; 
for i = 1:length(isodims)
%                 figure(i); clf;
    for j = randLoops
    %%
        fn = sprintf('%sshuffiso%d_%d.mat', ddir, isodims(i), j-1);
        if isfile(fn)
            temp = load(fn);
%             [i_val, e_val, sd_val, iso_val, ns] = i_e_powermod(temp);
%             props_e_mod(i,j) = e_val; %nansum(abs(e(s==smin)))/nansum(abs(e));
%             props_i_mod(i,j) = i_val;%nansum(abs(ipos(s==smin)))/nansum(abs(ipos));
%             props_sd_mod(i,j) = sd_val;%nansum(abs(ipos(s==smin)))/nansum(abs(ipos));
            
            if j==3
                e = normalize_rows(temp.ensem1);
                s = normalize_rows(temp.switch_binned);
                d = temp.isoMap;
                d = normalize_cols(d);
%                 figure(i+1000); clf; 
                figure(5001); %clf; 
                subplot(2,1,1); hold on 
                sep = i*2;
                plot(s+sep, 'k')
                plot(e+sep)
                title(isodims(i))
                subplot(2,1,2)
                hold on
                size(d)
%                 stacked_traces(d', .8)
                plot(sep + s, 'k')
                plot(sep + d(:,1)')
            end
        end
    end
end
drawnow

%%
nsegs = segs_range; % [1 5 10 25 50 100 250 500];
nsegs_e_mod = NaN(length(nsegs),nrand);
nsegs_i_mod = NaN(length(nsegs),nrand);
nsegs_sd_mod = NaN(length(nsegs),nrand);
nsegs_iso_mod = NaN(length(nsegs),nrand);
for i = 1:length(nsegs)
    for j = randLoops
        %%
        ddir = [topdir '\segs_mod\'];
        fn = sprintf('%sipos_nsegs%d_%d.mat', ddir, nsegs(i), j-1);
        if isfile(fn)
            temp = load(fn);
            [i_val, e_val, sd_val, iso_val, ns] = i_e_powermod(temp);
%             e1 = temp.ensem1;
%             e2 = temp.ensem2;
%             e = (e1-e2);
%             i1 = nanmean(abs(temp.ipos1), 1);
%             i2 = nanmean(abs(temp.ipos2), 1);
%             ipos = nanmean(abs(temp.ipos1) - abs(temp.ipos2), 1);
%             s = temp.switch_binned;
%             nsegs_e_mod(i,j) = nansum(abs(e(s==1)))/nansum(abs(e));
%             nsegs_i_mod(i,j) = nansum(abs(ipos(s==1)))/nansum(abs(ipos));
            nsegs_e_mod(i,j) = e_val;
            nsegs_i_mod(i,j) = i_val;
            nsegs_sd_mod(i,j) = sd_val;
            nsegs_iso_mod(i,j) = iso_val;
            
            if j==1
                d = temp.isoMap;
                d = normalize_cols(d);
                figure(i+2000); clf; 
                subplot(4,1,1)
                plot(temp.ensem1-temp.ensem2, 'b')
                title(nsegs(i))
                subplot(4,1,2:4)
                hold on
                stacked_traces(d')
%                 plot(d(:,1), 'k')
%                 plot(d(:,2)+1, 'k')
%                 plot(d(:,3)+2, 'k')
            end
        end
    end
end
drawnow
%%
ddir = [topdir '\nsample_mod\'];
nsamples = sample_range; % [50, 250, 500, 1000, 2500, 5000];
randLoops_ns = 1:124;
nsamples_e_mod = NaN(length(nsamples),length(randLoops_ns));
nsamples_i_mod = NaN(length(nsamples),length(randLoops_ns));
nsamples_sd_mod = NaN(length(nsamples),length(randLoops_ns));
nsamples_iso_mod = NaN(length(nsamples),length(randLoops_ns));
nsamples_binned = NaN(length(nsamples), 1);
for i = 1:length(nsamples)
    for j = randLoops_ns
        fn = sprintf('%sipos_nsample%d_%d.mat', ddir, nsamples(i), j-1);
        if isfile(fn)
            temp = load(fn);
            [i_val, e_val, sd_val, iso_val, ns] = i_e_powermod(temp);
            nsamples_e_mod(i,j) = e_val;% nansum(abs(e(s==1)))/nansum(abs(e));
            nsamples_i_mod(i,j) = i_val;%nansum(abs(ipos(s==1)))/nansum(abs(ipos));
            nsamples_sd_mod(i,j) = sd_val;%nansum(abs(ipos(s==1)))/nansum(abs(ipos));
            nsamples_iso_mod(i,j) = iso_val;%nansum(abs(ipos(s==1)))/nansum(abs(ipos));
            if isnan(nsamples_binned(i,1))
                nsamples_binned(i,1) = ns;%length(e);
            end
            if j==1
                d = temp.isoMap;
                d = normalize_cols(d);
                figure(i+3000); clf; 
                subplot(4,1,1)
                plot(temp.ensem1-temp.ensem2, 'b')
                title(nsamples(i))
                subplot(4,1,2:4)
                hold on
                plot(d(:,1), 'k')
                plot(d(:,2)+1, 'k')
                plot(d(:,3)+2, 'k')
            end
        end
    end
end
drawnow
%%
% nsamples = seg_mod_range; % [50, 250, 500, 1000, 2500, 5000];
randLoops_ns = 1:10;
seg_prop_e_mod = NaN(length(seg_mod_range),length(randLoops_ns));
seg_prop_i_mod = NaN(length(seg_mod_range),length(randLoops_ns));
seg_prop_sd_mod = NaN(length(seg_mod_range),length(randLoops_ns));
seg_prop_iso_mod = NaN(length(seg_mod_range),length(randLoops_ns));
for i = 1:length(seg_mod_range)
    for j = randLoops_ns
        %%
        ddir = [topdir '\seg_prop_mod\'];
        fn = sprintf('%sipos_segprop%d_%d.mat', ddir, seg_mod_range(i), j-1);
        if isfile(fn)
            temp = load(fn);
            e1 = temp.ensem1;
            e2 = temp.ensem2;
            e = (e1-e2);
%             plot(e2+.1*i)
            i1 = nanmean(abs(temp.ipos1), 1);
            i2 = nanmean(abs(temp.ipos2), 1);
            ipos = nanmean(abs(temp.ipos1) - abs(temp.ipos2), 1);
            sd = 2.^(temp.sd1);
            if j==1 && i==6
            figure(400+i); clf;
            subplot(4,2,3)
            plot(e, 'b')
            subplot(4,2,5)
            plot(ipos, 'k')
            subplot(4,2,7)
            plot(sd, 'r')
            subplot(4,2,1)
            imagesc(temp.bin_spks1, [0 1])
            subplot(122)
            scatter(ipos,sd)
            [c,p] = nancorr(ipos',sd');
            axis equal
            end
            if j==1
                d = temp.isoMap;
                d = normalize_matrix(d);
                figure(i+4000); clf; 
                subplot(4,1,1)
                plot(temp.ensem1-temp.ensem2, 'b')
                title(seg_mod_range(i))
                subplot(4,1,2:4)
                hold on
                stacked_traces(d')
            end
%             input('')
            [i_val, e_val, sd_val, iso_val, ~] = i_e_powermod(temp);
            seg_prop_e_mod(i,j) = e_val;% nansum(abs(e(s==1)))/nansum(abs(e));
            seg_prop_i_mod(i,j) = i_val;%nansum(abs(ipos(s==1)))/nansum(abs(ipos));
            seg_prop_sd_mod(i,j) = sd_val;%nansum(abs(ipos(s==1)))/nansum(abs(ipos));
            seg_prop_iso_mod(i,j) = iso_val;%nansum(abs(ipos(s==1)))/nansum(abs(ipos));
        end
    end
end

%%
stand_prop = 70;
stand_ncells = 512;%nsegs(end-2);
stand_nsamples = nsamples(end-1);
stand_prop_cells = 100;
seg_prop = seg_mod_range;

vars2plot = {'props', 'nsegs', 'nsamples', 'seg_prop'};
% dep2plot = {'i_mod', 'sd_mod', 'e_mod'};
dep2plot = {'i_mod', 'sd_mod', 'e_mod', 'iso_mod'};
clrs2plot = [.2 .2 .2; .6 .2 .2; .2 .2 .8; .8 .2 .8];
stand_clr = [.8 1 .8];
% iclr = [.2 .2 .2];
% eclr = [.2 .2 .8];
medge = 'none';
msize = 10;
figure(102); clf
subplot(2,2, 1); hold on
s = rectangle('Position', [(stand_prop)-2.5 -.05 5 2]); s.FaceColor = stand_clr; s.EdgeColor = 'none';
plot([0 100], [1 0], 'k:', 'LineWidth', 1)
gt = [props(1:4:end)];
set(gca, 'XTick', (gt), 'XTickLabel', gt)
axis([-10 110 -.1 1.1])
xlabel('Percent samples modulated')

subplot(2,2, 2); cla; hold on
xs = log2(nsegs);
s = rectangle('Position', [log2(stand_ncells)-.2 -.05 .4 2]); s.FaceColor = stand_clr; s.EdgeColor = 'none';
plot([-10 110], 1-[stand_prop/100 stand_prop/100], 'k:', 'LineWidth', 1)
gt = [nsegs(1:4:end) nsegs(end)];
set(gca, 'XTick', log2(gt), 'XTickLabel', gt)
axis([-1 log2(nsegs(end)*2) -.1 1.1])
xlabel('Number of cells')

subplot(2,2, 3); cla; hold on
xs = log2(nsamples);
s = rectangle('Position', [log2(stand_nsamples)-.2 -.05 .4 2]); s.FaceColor = stand_clr; s.EdgeColor = 'none';
plot([-10 110], 1-[stand_prop/100 stand_prop/100], 'k:', 'LineWidth', 1)
gt = [nsamples(2:2:end-2) nsamples(end)];
gtl = bin_timeres*[nsamples_binned(2:2:end-2, 1); nsamples_binned(end, 1)];
set(gca, 'XTick', log2(gt), 'XTickLabel', gtl)
axis([log2(nsamples(1)/2) log2(nsamples(end)*2) -.1 1.1])
xlabel('Total time (sec)')

subplot(2,2, 4); cla; hold on
s = rectangle('Position', [(stand_prop_cells)-2.5 -.05 5 2]); s.FaceColor = stand_clr; s.EdgeColor = 'none';
plot([-10 110], 1-[stand_prop/100 stand_prop/100], 'k:', 'LineWidth', 1)
gt = [0:20:100];
set(gca, 'XTick', (gt), 'XTickLabel', gt)
axis([-10 110 -.1 1.1])
xlabel('Percent cells modulated')

for varLoop = 1:length(vars2plot)
    subplot(2,2, varLoop); hold on
        x = eval(sprintf('%s', vars2plot{varLoop}));
    if varLoop==2 || varLoop == 3
        x = log2(x);
    else
    end
    for depLoop = 1:length(dep2plot)
        y = eval(sprintf('%s_%s', vars2plot{varLoop}, dep2plot{depLoop}));
        c = clrs2plot(depLoop,:);
        plot(x, nanmean(y, 2), 'Color', shift_colormap(c, 2), 'LineWidth', 3)
        for i = 1:length(x)
            xs = x(i)*ones(length(y(i,:)),1);
            scatter(xs, y(i,:), msize, 'MarkerFaceColor', c, 'MarkerEdgeColor', medge, 'MarkerFaceAlpha', .8);
        end
    end
end
%%
function [i_val, e_val, sd_val, iso_val, ns] = i_e_powermod(temp_struct)
sd = 2.^(temp_struct.sd1);% - abs(temp_struct.sd2); % sd1 and sd2 are essentially the same
e1 = temp_struct.ensem1;
e2 = temp_struct.ensem2;
e = (e1-e2);
i1 = nanmean(abs(temp_struct.ipos1), 1);
i2 = nanmean(abs(temp_struct.ipos2), 1);
ipos = nanmean(abs(temp_struct.ipos1) - abs(temp_struct.ipos2), 1);
s = temp_struct.switch_binned;
iso = normalize_cols(temp_struct.isoMap);
iso = iso(:,1)';

e = abs(e);
ipos = abs(ipos);
e_val = nansum(e(s==1))/nansum(e);
sd_val = nansum(sd(s==1))/nansum(sd);
i_val = nansum(ipos(s==1))/nansum(ipos);
iso_val = nansum(iso(s==1))/nansum(iso);
%             i_val = nansum(abs(ipos(s==1)))/nansum(abs(ipos));
ns = length(e);
end
            
