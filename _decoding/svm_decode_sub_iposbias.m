function [str_in] = svm_decode_sub_iposbias(str_in, ipos_mean, predict_var, predict_bins, spks, x_fold_training, distance_method)
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

pos_bias = ipos_mean>0;

str_in.pred_same        = NaN(nsamples, 1);
str_in.pred_diff        = NaN(nsamples, 1);
str_in.pred_same_rand   = NaN(nsamples, 1);
str_in.pred_diff_rand   = NaN(nsamples, 1);
str_in.pred_err_same    = NaN(nsamples, 1);
str_in.pred_err_diff    = NaN(nsamples, 1);
str_in.rand_err_same    = NaN(nsamples, 1);
str_in.rand_err_diff    = NaN(nsamples, 1);
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
%%%% Training and testing based on the room/arena ipos signal
%
test_inds               = NaN(nsamples, 1);
pos_bias = ipos_mean>0;
pos_inds = find(ipos_mean>0);
neg_inds = find(ipos_mean<=0);

labels = [-1*x_fold_training:-1 1:x_fold_training];
% x_fold_training = length(labels);

test_inds(pos_inds) = mod(1:length(pos_inds), x_fold_training);
test_inds(test_inds==0) = x_fold_training;
test_inds(neg_inds) = -1*mod(1:length(neg_inds), x_fold_training);
test_inds(test_inds==0) = -1*x_fold_training;

% for i = 1:length(str_in.counts)
%     ind = find(str_in.var_binned==i);
%     [~, ranord] = sort(rand(length(ind),1));
%     ind = ind(ranord);
%     xi = floor(linspace(0, length(ind), x_fold_training+1));
%     for j = 1:x_fold_training
%         test_inds(ind(xi(j)+1:xi(j+1))) = j;
%     end
% end
% str_in.test_inds = test_inds;
%%%% Construct the weights based on inverse frequency of occupancy
sample_weights = 1-(str_in.counts./sum(str_in.counts));
wmat = ones(nsamples, 1) * sample_weights ; 
str_in.weights = NaN(nsamples, 1);
for i = 1:length(str_in.counts)
    ind = str_in.var_binned==i;
    str_in.weights(ind) = wmat(ind, i);
end
bad_ipos = isnan(ipos_mean);
bad_pred = isnan(sum(spks,1)) | bad_ipos;
% bad_rand_pred = isnan(sum(spks_rand,1));
bad_rand_pred = isnan(sum(spks_rand,1)) | bad_ipos;
test_inds2 = test_inds;
test_inds2(bad_pred) = NaN;
% test_inds3 = test_inds;
% test_inds3(bad_rand_pred) = NaN;
%
% fprintf('Decoding using xfold=%d;  ', x_fold_training)
% labels = [-2 -1 1 2];
%%
for i = 1:length(labels)
    % using SVM to decode
    istesting = test_inds2 == labels(i); % get the training inds    
    test_spks = spks(:, istesting)';
    
    train_same = test_inds2/0 == labels(i)/0 ; % get the same sign for 'same'
    train_same(test_inds2 == labels(i)) = 0;
    train_diff = test_inds2/0 == -labels(i)/0 ; % get the opposite sign for 'diff'
    train_diff(test_inds2 == -labels(i)) = 0;
    
        train_s = spks(:, train_same)';
        train_x = str_in.var_binned(train_same);
        train_w = str_in.weights(train_same);
    Mdl_same = fitcecoc(train_s, train_x, 'Weights', train_w, 'Coding', 'onevsall');
        train_s = spks(:, train_diff)';
        train_x = str_in.var_binned(train_diff);
        train_w = str_in.weights(train_diff);
    Mdl_diff = fitcecoc(train_s, train_x, 'Weights', train_w, 'Coding', 'onevsall');
    
    test_spks_rand = spks_rand(:, istesting)';
    
    fprintf(' ..')
    % onevsall seems faster and better than onevsone
    str_in.pred_same(istesting) = predict(Mdl_same, test_spks);
    str_in.pred_diff(istesting) = predict(Mdl_diff, test_spks);
    
    str_in.pred_same_rand(istesting) = predict(Mdl_same, test_spks_rand);
    str_in.pred_diff_rand(istesting) = predict(Mdl_diff, test_spks_rand);
    
    % random shuffle comparison
    fprintf('%d%%', round(100*(i/length(labels))))
end
fprintf(' Done!\n')
%%
valid_bin = find(~bad_pred);
valid_bin_rand = find(~bad_rand_pred);
if decode_2d == true
    % 2D eucidean distance
    [angbin_same, radbin_same] = ind2sub(size(angcounts), str_in.pred_same);
    [angbin_diff, radbin_diff] = ind2sub(size(angcounts), str_in.pred_diff);
    [angbin_rand_same, radbin_rand_same] = ind2sub(size(angcounts), str_in.pred_same_rand);
    [angbin_rand_diff, radbin_rand_diff] = ind2sub(size(angcounts), str_in.pred_diff_rand);
    [predx, predy] = pol2cart( predict_var{1}, predict_var{2});
    
    posx = predx*NaN; posy = posx;
    posx_r = posx; posy_r = posx;
    
    [posx(valid_bin), posy(valid_bin)] = ...
        pol2cart( str_in.bin_center{1}(angbin_same(valid_bin)), str_in.bin_center{2}(radbin_same(valid_bin)));
    str_in.pred_err_same = sqrt((predx-posx).^2 + (predy-posy).^2);
    [posx(valid_bin), posy(valid_bin)] = ...
        pol2cart( str_in.bin_center{1}(angbin_diff(valid_bin)), str_in.bin_center{2}(radbin_diff(valid_bin)));
    str_in.pred_err_diff = sqrt((predx-posx).^2 + (predy-posy).^2);
    
    
    [posx_r(valid_bin_rand), posy_r(valid_bin_rand)] = ...
        pol2cart( str_in.bin_center{1}(angbin_rand_same(valid_bin_rand)), str_in.bin_center{2}(radbin_rand_same(valid_bin_rand)));
    str_in.rand_err_same = sqrt((predx-posx_r).^2 + (predy-posy_r).^2);
    [posx_r(valid_bin_rand), posy_r(valid_bin_rand)] = ...
        pol2cart( str_in.bin_center{1}(angbin_rand_diff(valid_bin_rand)), str_in.bin_center{2}(radbin_rand_diff(valid_bin_rand)));
    str_in.rand_err_diff = sqrt((predx-posx_r).^2 + (predy-posy_r).^2);
else
    pos = predict_var*NaN; 
    pred_pos_same = pos; pred_pos_diff = pos; 
    rand_pos_same = pos; rand_pos_diff = pos; 
    pred_pos_same(valid_bin) = str_in.bin_center(str_in.pred_same(valid_bin));
    pred_pos_diff(valid_bin) = str_in.bin_center(str_in.pred_diff(valid_bin));
    rand_pos_same(valid_bin_rand) = str_in.bin_center(str_in.pred_same_rand(valid_bin_rand));
    rand_pos_diff(valid_bin_rand) = str_in.bin_center(str_in.pred_diff_rand(valid_bin_rand));
    switch distance_method
        case 'polar_dist_real'
            pos(:) = predict_var;
            angdiff = pos - pred_pos_same;
            str_in.pred_err_same = (mod(angdiff + pi, 2*pi) - pi);
            angdiff = pos - pred_pos_diff;
            str_in.pred_err_diff = (mod(angdiff + pi, 2*pi) - pi);
            
            % random difference
            angdiff = pos - rand_pos_same;
            str_in.rand_err_same = (mod(angdiff + pi, 2*pi) - pi);
            angdiff = pos - rand_pos_diff;
            str_in.rand_err_diff = (mod(angdiff + pi, 2*pi) - pi);
% %         case 'polar_dist_binned'
% %             pos(:) = str_in.bin_center(str_in.var_binned);
% %             % predicted difference
% %             angdiff = pos - pred_pos;
% %             str_in.pred_err = (mod(angdiff + pi, 2*pi) - pi);
% %             
% %             % random difference
% %             angdiff = pos - pred_rand;
% %             str_in.rand_err = (mod(angdiff + pi, 2*pi) - pi);
% %         case 'eucl_dist_real'
% %             pos(:) = predict_var;
% %             % predicted difference
% %             str_in.pred_err = sqrt((pos - pred_pos).^2);
% %             
% %             % random difference
% %             str_in.rand_err = sqrt((pos - pred_rand).^2);
% %         case 'eucl_dist_binned'
% %             pos(:) = str_in.bin_center(str_in.var_binned);
% %             % predicted difference
% %             pred_pos(:) = str_in.bin_center(str_in.pred_real);
% %             str_in.pred_err = sqrt((pos - pred_pos).^2);
% %             
% %             % random difference
% %             str_in.rand_err = sqrt((pos - pred_rand).^2);
    end
end