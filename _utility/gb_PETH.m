function [peth_mean, peth_std, peth_pre, num_events, peth_p, peth_diff, P] = gb_PETH(sig, event, pre, post)

if islogical(sig)
%     bin_e = 1:8:length(sig);
%     bin_width = mean(diff(bin_e))/2;
%     bin_center = bin_e(1:end-1) + bin_width;
%     
%     ts = find(sig);
%     [binLP, ~, ~] = histcounts(ts, bin_e);
%     ts_hs = interp1(bin_center, 1*(binLP), 1:length(sig));
%     ts_hs(isnan(ts_hs)) = 0;
%     sig = ts_hs;
else
    sig = sig/max(sig);
    
end

if size(sig,1) > size(sig,2)
    sig = sig';
end
sig_pad = [NaN(1, pre+1) sig NaN(1, post+1)];

if size(event,1) > size(event,2)
    event = event';
end
    eventdiff = [0 diff(event)];

% event_pad = [NaN(1, pre+1) event NaN(1, post+1)];
eventdiff_pad = [NaN(1, pre+1) eventdiff NaN(1, post+1)];

% et = find(event_pad>0);
et  = find(eventdiff_pad>0);

P = NaN(pre+post+1, length(et));
for i = 1:length(et)
    pre_ind = et(i)-pre;
    post_ind = et(i)+post+1;
    P(:, i) = sig_pad(pre_ind:post_ind-1);
end
% P_pre = P(pre-29:pre-8, :);
% P_post = P(pre+8:pre+29, :);
P_pre = P(1:pre, :);
P_post = P(pre+1:end, :);
pre_dist  = sum(P_pre(:,  sum([P_pre; P_post], 1)>.000001), 1);
post_dist = sum(P_post(:, sum([P_pre; P_post], 1)>.000001), 1);
% [~, peth_p, ~, stats] = ttest(mean(P(pre-14:pre, :),2), mean(P(pre+1:pre+15,:),2));
% [~, peth_p, ~, stats] = ttest( pre_dist,  post_dist);
% peth_t = stats.tstat;
if ~isempty(pre_dist) || ~isempty(pre_dist)
    peth_p = signrank(pre_dist, post_dist);
    peth_diff = [sum(pre_dist) - sum(post_dist)];
else 
    peth_p = 0;
    peth_diff = [sum(P_pre(:)) - sum(P_post(:))];
end
peth_mean = nanmean(P,2);
peth_pre = nanmean(P(1:pre, :), 1);
num_events = size(P, 2);
peth_std = nanstd(P, [], 2);
