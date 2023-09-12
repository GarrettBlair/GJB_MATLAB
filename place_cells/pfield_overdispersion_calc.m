%% function [od] = pfield_overdispersion_calc(spks, placefields, x, y, t)
spike_matrix = h.ms.spks; 
t = h.ms.timestamps./1000;
dt = abs(diff(t)); dt = [dt(1); dt];
placefields = h.ms.room.pfields_smooth;
[nsegs, nbins,~] = size(placefields);
nsamples = length(t);
bins = h.ms.params.pos_bins;
x = h.ms.room.x;
y = h.ms.room.y;
spd = h.ms.arena.speed_smooth;

linear_pfields = reshape(placefields, [nsegs, nbins^2]);
linear_pfields(isnan(linear_pfields)) = 0;
linear_pfields = normalize_rows(linear_pfields);
linear_pfields(linear_pfields<.2) = 0;
pf = reshape(linear_pfields, [nsegs, nbins, nbins]);


[~, ~, ~, xbin, ybin] = histcounts2(x,y,bins,bins);
bads = isnan(xbin) | isnan(ybin) | (xbin==0) | (ybin==0);
if any(bads)
    xbin(bads) = interp1(find(~bads), xbin(~bads), find(bads), 'nearest', 'extrap');
    ybin(bads) = interp1(find(~bads), ybin(~bads), find(bads), 'nearest', 'extrap');
end
lin_ind = sub2ind([nbins nbins], ybin, xbin);
[yy,xx] = ind2sub([nbins nbins], lin_ind);
ex_spks = NaN(nsegs, nsamples);
for i = 1:nsamples
    ex_spks(:, i) = linear_pfields(:, lin_ind(i));
end
ex_spks(isnan(ex_spks)) = 0;
spike_matrix(isnan(spike_matrix)) = 0;

% ex_spks = ex_spks(:, spd);
spks = spike_matrix;

sum_t   = bin_spks_time(dt',    5, t', false);
av_spd = average_spks_time(spd',    5, t, false, 'mean');
av_exp = average_spks_time(ex_spks, 5, t, false, 'mean');
% av_obs = average_spks_time(spks,    5, t, false, 'mean');
% av_exp = bin_spks_time(ex_spks, 5, t', false);
av_obs = bin_spks_time(spks,    5, t', false);
% t_mat = ones(size(av_exp))*sum_t;
% t_mat = ones(nsegs,1)*sum_t;
av_exp = av_exp./t_mat;
av_obs = av_obs./t_mat;
av_exp = av_exp(:,av_spd>=5);
av_obs = av_obs(:,av_spd>=5);
%%
n_ds = size(av_obs,2);
rate_thresh = h.ms.room.pcell_stats.peakRate;
% rr = max(av_obs,[],2);
ratemat = rate_thresh*ones(1,n_ds);
% sum(av_exp>=me, 2)./n_ds;

z = NaN(nsegs,n_ds);
goods = (av_exp>=ratemat);
for i = 1:n_ds
%     goods = av_exp(i,:)>=me(i,:);
    o = av_obs(:,i);
    e = av_exp(:,i);
    z(goods(:,i), i) = (o(goods(:,i))-e(goods(:,i)))./(e(goods(:,i)).^.5);
end



%%
for i = 1:15%nsegs
    figure(89); clf
    histogram(av_obs(i,:), [0:.05:1], 'Normalization', 'probability');
    drawnow
    pause(.5)
end














