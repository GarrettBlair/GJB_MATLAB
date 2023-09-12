function [str_in] = svm_decode_sub(str_in, predict_var, predict_bins, spks, x_fold_training, distance_method)
% spks = str_in.spks;
%%
n_rand     = str_in.num_random_shuffle;
train_flag = str_in.train_flag;
parfor_progbar = str_in.parfor_progbar;

[nsegs, nsamples] = size(spks);
valid_shifts = [round(-.9*nsamples):round(-.1*nsamples), round(.1*nsamples):round(.9*nsamples)];
shift_ind = randi([1, length(valid_shifts)]);
spks_rand = spks;
% % for i = 1:nsegs
% %     spks_rand(i,:) = circshift(spks_rand(i,:), valid_shifts(shift_ind));
% % end
% spks_rand = spks(ranord,:);
% [~, rand_ord] = sort(rand(nsegs,1));
% spks_rand = spks_rand(rand_ord,:);
str_in.pred_real        = NaN(nsamples, 1);
str_in.pred_rand        = NaN(nsamples, n_rand);
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
% str_in.real_post        = NaN(nsamples, length(unique(str_in.var_binned)));
% str_in.rand_post        = NaN(nsamples, length(unique(str_in.var_binned)));
%%%% Training and testing indices based by equalizing the count histograms
% for i = 1:length(str_in.counts)
%     ind = find(str_in.var_binned==i);
%     [~, ranord] = sort(rand(length(ind),1));
%     ind = ind(ranord);
%     xi = floor(linspace(0, length(ind), x_fold_training+1));
%     for j = 1:x_fold_training
%         test_inds(ind(xi(j)+1:xi(j+1))) = j;
%     end
% end
%%%%%%%%%%%%%%%%%
odd_ind = [];
split_equal = 1;
for i = 1:length(str_in.counts)
    ind = find(str_in.var_binned==i);
    [~, randord] = sort(rand(length(ind),1));
    ind = ind(randord);
    if mod(length(ind), 2) == 1
        % this is to keep the last split value from always geting one extra
        % when there's an odd number
        odd_ind = ind(end);
        ind = ind(1:end-1);
    else
        odd_ind = [];
    end
    if ~isempty(ind)
        xi = floor(linspace(0, length(ind), x_fold_training+1));
        for j = 1:x_fold_training
            test_inds(ind(xi(j)+1:xi(j+1))) = j;
        end
    end
    if ~isempty(odd_ind)
        test_inds(odd_ind) = split_equal;
        split_equal = mod(split_equal, x_fold_training)+1;
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
no_spks = nansum(spks,1)==0;
no_spks_rand = nansum(spks_rand,1)==0;
%%
% fprintf('Decoding using xfold=%d;  ', x_fold_training)
% warning('off', 'stats:fitSVMPosterior:PerfectSeparation')
% warning('off', 'stats:ClassificationECOC:localFitECOC:CannotFitPosterior')
% % Mdl = fitcecoc(spks', str_in.var_binned, 'Weights', str_in.weights, 'FitPosterior', true);
% % [label1, loss1, score1, pred1] = resubPredict(Mdl);
% % 
% % Mdl = fitcecoc(spks_rand', str_in.var_binned, 'Weights', str_in.weights, 'FitPosterior', true);
% % [label2, loss2, score2, pred2] = resubPredict(Mdl);
svmtemplate = templateSVM('Standardize', false);
[~, rand_spks_ord] = sort(rand(nsegs, n_rand));
ppm = 0;
for i = 1:x_fold_training
    % using SVM to decode
    istesting       = (test_inds==i & no_spks'==false);
    istraining      = (test_inds~=i & no_spks'==false & train_flag'==true);
    istesting_rand  = (test_inds==i & no_spks_rand'==false);
    train_s         = spks(:, istraining)';
    train_x         = str_in.var_binned(istraining);
    train_w         = str_in.weights(istraining);
    test_s          = spks(:, istesting)';
    test_x          = str_in.var_binned(istesting)';
    
    fprintf(' ..')
    Mdl = fitcecoc(train_s, train_x, 'Learners', svmtemplate, 'Weights', train_w, 'FitPosterior', false);
%     [templabel,negloss, tempscore, temppred] = resubPredict(Mdl);
%     [str_in.pred_real(istesting),~,~,str_in.real_post(istesting,:)] = predict(Mdl, test_s);
%     [str_in.pred_rand(istesting),~,~,str_in.rand_post(istesting_rand,:)] = predict(Mdl, test_s_rand);
    [str_in.pred_real(istesting)] = predict(Mdl, test_s);
    % random shuffle comparison
    test_s_rand     = spks_rand(:, istesting_rand)';
    
    rand_pred = NaN(sum(istesting_rand), n_rand);
    
    if parfor_progbar == true
        ppm = ParforProgressbar(n_rand,'parpool', {'local'}, 'showWorkerProgress',true,...
        'progressBarUpdatePeriod',5,'title','SVM Decoding randLoop');
    end
    parfor randLoop = 1:n_rand
        [rand_pred(:,randLoop)] = predict(Mdl, test_s_rand(:, rand_spks_ord(:, randLoop)));
        if parfor_progbar == true
            pause(.001)
            ppm.increment();
        end
    end
    for randLoop = 1:n_rand
        str_in.pred_rand(istesting_rand, randLoop) = rand_pred(:,randLoop);
    end
    fprintf('%d%%', round(100*(i/x_fold_training)))
end
%%
fprintf(' Done!\n')
valid_bin = ~no_spks;% find(~bad_pred);
valid_bin_rand = ~no_spks_rand;% find(~bad_rand_pred);
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
    
    for randLoop = 1:n_rand
        r = str_in.pred_rand(:, randLoop);
        [angbin_rand, radbin_rand] = ind2sub(size(angcounts), r);
        [px(valid_bin_rand), py(valid_bin_rand)] = ...
            pol2cart( str_in.bin_center{1}(angbin_rand(valid_bin_rand)), str_in.bin_center{2}(radbin_rand(valid_bin_rand)));
        dx = sqrt((ax-px).^2 + (ay-py).^2);
        str_in.rand_err_median(randLoop) = nanmedian(dx);
        str_in.rand_err_q25(randLoop) = quantile(dx, .25);
        str_in.rand_err_q75(randLoop) = quantile(dx, .75);
        if randLoop ==1
            str_in.rand_dist_ex = dx;
            str_in.rand_x = px;
            str_in.rand_y = py;
        end
    end
else
    a = predict_var*NaN; ap = a; apr = a;
    ap(valid_bin) = str_in.bin_center(str_in.pred_real(valid_bin));
%     apr(valid_bin_rand) = str_in.bin_center(str_in.pred_rand(valid_bin_rand));
    switch distance_method
        case 'polar_dist_real'
            a(:) = predict_var;
            angdiff = a - ap;
            str_in.pred_err = (mod(angdiff + pi, 2*pi) - pi);
            % random difference
            for randLoop = 1:n_rand
                apr(valid_bin_rand) = str_in.bin_center(str_in.pred_rand(valid_bin_rand, randLoop));
                angdiff = a - apr;
                str_in.rand_err_median(randLoop) = (mod(angdiff + pi, 2*pi) - pi);
            end
        case 'polar_dist_binned'
            a(:) = str_in.bin_center(str_in.var_binned);
            % predicted difference
            angdiff = a - ap;
            str_in.pred_err = (mod(angdiff + pi, 2*pi) - pi);
            % random difference
            for randLoop = 1:n_rand
                apr(valid_bin_rand) = str_in.bin_center(str_in.pred_rand(valid_bin_rand, randLoop));
                angdiff = a - apr;
                str_in.rand_err_median(randLoop) = (mod(angdiff + pi, 2*pi) - pi);
            end
        case 'eucl_dist_real'
            a(:) = predict_var;
            % predicted difference
            str_in.pred_err = sqrt((a - ap).^2);
            % random difference
            for randLoop = 1:n_rand
                apr(valid_bin_rand) = str_in.bin_center(str_in.pred_rand(valid_bin_rand, randLoop));
                str_in.rand_err_median(randLoop) = sqrt((a - apr).^2);
            end
        case 'eucl_dist_binned'
            a(:) = str_in.bin_center(str_in.var_binned);
            % predicted difference
            ap(:) = str_in.bin_center(str_in.pred_real);
            str_in.pred_err = sqrt((a - ap).^2);
            % random difference
            for randLoop = 1:n_rand
                apr(valid_bin_rand) = str_in.bin_center(str_in.pred_rand(valid_bin_rand, randLoop));
                str_in.rand_err_median(randLoop) = sqrt((a - apr).^2);
            end
    end
end

str_in.pred_rand = str_in.pred_rand(:, 1); % to reduce size of save file; just keep one example

    
end