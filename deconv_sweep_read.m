function [smat, smat_weighted, good_idx, good_idx_weighted] = deconv_sweep_read(fname, smin)
%%
load(fname);
i = 1;
if smin(i)==0
    varName = sprintf('smin_None');
    smin_weight = 1;
elseif smin(i)<0
    varName = sprintf('smin_neg%d', abs(smin(i)));
    smin_weight = abs(smin(i));
else
    varName = sprintf('smin_%d', smin(i));
    smin_weight = abs(smin(i));
end
template_s = eval(varName);

smat              = false(length(smin), size(template_s,1), size(template_s,2));
smat_weighted     = zeros(length(smin), size(template_s,1), size(template_s,2));
good_idx          = false(length(smin), size(template_s,1));
good_idx_weighted = zeros(length(smin), size(template_s,1));

for i = 1:length(smin)
    if smin(i)==0
        varName = sprintf('smin_None');
        idxGoodName = sprintf('smin_None_idx_components');
        smin_weight = 1;
    elseif smin(i)<0
        varName = sprintf('smin_neg%d', abs(smin(i)));
        idxGoodName = sprintf('smin_neg%d_idx_components', abs(smin(i)));
        smin_weight = abs(smin(i));
    else
        varName = sprintf('smin_%d', smin(i));
        idxGoodName = sprintf('smin_%d_idx_components', abs(smin(i)));
        smin_weight = abs(smin(i));
    end
    s       = eval(varName);
    idx_g   = eval(idxGoodName);
    good_idx(i, idx_g+1) = true;
    good_idx_weighted(i, idx_g+1) = smin_weight;
    smat(i,:,:) = s;% = cat(3, smat, s);
    smat_weighted(i,:,:) = s*smin_weight;% = cat(3, smat_weighted, s*smin_weight);
end
S1 = squeeze(sum(smat,1));
% S1 = S1(ms.neuron.idx_components,:);
S2 = squeeze(sum(smat_weighted,1));
% S2 = S2(ms.neuron.idx_components,:);
% figure(1); clf
% stacked_traces(S1, 1, {'r-', 'LineWidth', 1.5});
% stacked_traces(S2, 1, {'b-', 'LineWidth', 1});
