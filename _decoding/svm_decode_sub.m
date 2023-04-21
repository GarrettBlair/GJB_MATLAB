function [str_in] = svm_decode_sub(str_in, predict_var, predict_bins, spks, x_fold_training, distance_method)
% spks = str_in.spks;
%%
[nsegs, nsamples] = size(spks);
[~, ranord] = sort(rand(nsegs, 1));
spks_rand = spks(ranord,:);
valid_shifts = [round(-.9*nsamples):round(-.1*nsamples), round(.1*nsamples):round(.9*nsamples)];
shift_ind = randi([1, length(valid_shifts)]);
for i = 1:nsegs
%     shiftval = randi([-1*floor(nsamples/2) 1*floor(nsamples/2)]);
    spks_rand(i,:) = circshift(spks_rand(i,:), valid_shifts(shift_ind));
end

str_in.pred_real        = NaN(nsamples, 1);
str_in.pred_rand        = NaN(nsamples, 1);
str_in.var_binned       = NaN(nsamples, 1);
str_in.test_inds        = NaN(nsamples, 1);
test_inds               = NaN(nsamples, 1);

decode_2d = iscell(predict_var) == true;
if decode_2d == false 
[str_in.counts, ~,  str_in.var_binned(:)] = histcounts(predict_var, predict_bins);
    str_in.bin_center = predict_bins(2:end) - abs(diff(predict_bins))/2;

else
    [angcounts, ~, ~,  angbin, radbin] = histcounts2(predict_var{1}, predict_var{2}, predict_bins{1}, predict_bins{2});
    str_in.var_binned = sub2ind(size(angcounts), angbin, radbin);
    str_in.counts = angcounts(:)';
    str_in.bin_center = [];
    str_in.bin_center{1} = predict_bins{1}(2:end) - abs(diff(predict_bins{1}))/2;
    str_in.bin_center{2} = predict_bins{2}(2:end) - abs(diff(predict_bins{2}))/2;
end
%%%% Training and testing indices based by equalizing the count histograms
for i = 1:length(str_in.counts)
    ind = find(str_in.var_binned==i);
    [~, ranord] = sort(rand(length(ind),1));
    ind = ind(ranord);
    xi = floor(linspace(0, length(ind), x_fold_training+1));
    for j = 1:x_fold_training
        test_inds(ind(xi(j)+1:xi(j+1))) = j;
    end
end
str_in.test_inds = test_inds;
%%%% Construct the weights based on inverse frequency
sample_weights = 1 - (str_in.counts./sum(str_in.counts));
wmat = ones(nsamples, 1) * sample_weights ; 
str_in.weights = NaN(nsamples, 1);
for i = 1:length(str_in.counts)
    ind = str_in.var_binned==i;
    str_in.weights(ind) = wmat(ind, i);
end
bad_pred = isnan(sum(spks,1));
bad_rand_pred = isnan(sum(spks_rand,1));
test_inds2 = test_inds;
test_inds2(bad_pred) = NaN;
test_inds3 = test_inds;
test_inds3(bad_rand_pred) = NaN;
%%
% fprintf('Decoding using xfold=%d;  ', x_fold_training)
for i = 1:x_fold_training
    % using SVM to decode
    istesting  = find(test_inds2==i);
    istraining = find(test_inds2~=i & ~isnan(test_inds2));
    train_s = spks(:, istraining)';
    train_x = str_in.var_binned(istraining);
    train_w = str_in.weights(istraining);
    test_s = spks(:, istesting)';
    test_x = str_in.var_binned(istesting)';

    fprintf(' ..')
    Mdl = fitcecoc(train_s, train_x, 'Weights', train_w);
    str_in.pred_real(istesting) = predict(Mdl, test_s);
    
    % random shuffle comparison
    istesting_rand  = find(test_inds3==i);
    istraining_rand = find(test_inds3~=i & ~isnan(test_inds3));
    train_s = spks_rand(:, istraining_rand)';
    train_x = str_in.var_binned(istraining_rand);
    train_w = str_in.weights(istraining_rand);
    test_s = spks_rand(:, istesting_rand)';


    Mdl = fitcecoc(train_s, train_x, 'Weights', train_w);
    str_in.pred_rand(istesting_rand) = predict(Mdl, test_s);
    fprintf('%d%%', round(100*(i/x_fold_training)))
end
fprintf(' Done!\n')
valid_bin = find(~bad_pred);
valid_bin_rand = find(~bad_rand_pred);
if decode_2d == true
    % 2D eucidean distance
    [angbin_real, radbin_real] = ind2sub(size(angcounts), str_in.pred_real);
    %     [ax, ay] = pol2cart( str_in.bin_center{1}(angbin), str_in.bin_center{2}(radbin));
    [ax, ay] = pol2cart( predict_var{1}, predict_var{2});
    
    px = ax*NaN; py = px;
    pxr = ax; pyr = px;
    
    [px(valid_bin), py(valid_bin)] = ...
        pol2cart( str_in.bin_center{1}(angbin_real(valid_bin)), str_in.bin_center{2}(radbin_real(valid_bin)));
    str_in.pred_err = sqrt((ax-px).^2 + (ay-py).^2);
    str_in.pred_x = px;
    str_in.pred_y = py;
    
    [angbin_rand, radbin_rand] = ind2sub(size(angcounts), str_in.pred_rand);
    [px(valid_bin_rand), py(valid_bin_rand)] = ...
        pol2cart( str_in.bin_center{1}(angbin_rand(valid_bin_rand)), str_in.bin_center{2}(radbin_rand(valid_bin_rand)));
%     [px, py] = pol2cart( str_in.bin_center{1}(angbin_rand), str_in.bin_center{2}(radbin_rand));
    str_in.rand_err = sqrt((ax-px).^2 + (ay-py).^2);
    str_in.rand_x = px;
    str_in.rand_y = py;
else
    a = predict_var*NaN; ap = a; apr = a;
    ap(valid_bin) = str_in.bin_center(str_in.pred_real(valid_bin));
    apr(valid_bin_rand) = str_in.bin_center(str_in.pred_rand(valid_bin_rand));
    switch distance_method
        case 'polar_dist_real'
            a(:) = predict_var;
            angdiff = a - ap;
%             str_in.pred_err = abs(mod(angdiff + pi, 2*pi) - pi);
            str_in.pred_err = (mod(angdiff + pi, 2*pi) - pi);
            
            % random difference
            angdiff = a - apr;
            str_in.rand_err = (mod(angdiff + pi, 2*pi) - pi);
        case 'polar_dist_binned'
            a(:) = str_in.bin_center(str_in.var_binned);
            % predicted difference
            angdiff = a - ap;
            str_in.pred_err = (mod(angdiff + pi, 2*pi) - pi);
            
            % random difference
            angdiff = a - apr;
            str_in.rand_err = (mod(angdiff + pi, 2*pi) - pi);
        case 'eucl_dist_real'
            a(:) = predict_var;
            % predicted difference
            str_in.pred_err = sqrt((a - ap).^2);
            
            % random difference
            str_in.rand_err = sqrt((a - apr).^2);
        case 'eucl_dist_binned'
            a(:) = str_in.bin_center(str_in.var_binned);
            % predicted difference
            ap(:) = str_in.bin_center(str_in.pred_real);
            str_in.pred_err = sqrt((a - ap).^2);
            
            % random difference
            str_in.rand_err = sqrt((a - apr).^2);
    end
end