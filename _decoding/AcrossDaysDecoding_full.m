function [full, random ] = AcrossDaysDecoding_full(trainday, testday, segsday1, segsday2, rand_flag)
train_set = load(trainday);

params = train_set.params;
params.decoding_bin_edges = -40:80/12:40;
params.decoding_nbins = length(params.decoding_bin_edges) - 1;
conv_kern = ones(1,7); % temporal convolution window to spread out spike info, ~1 sec
train_spk = train_set.spike_mat(segsday1,:);
conv_train_spk = nan(size(train_spk));
for seg = 1:size(train_spk,1)
    conv_train_spk(seg, :) = conv(train_spk(seg,:), conv_kern, 'same');
    conv_train_spk(seg, conv_train_spk(seg,:)<1) = 0;
    conv_train_spk(seg, conv_train_spk(seg,:)>=1) = 1;
end
xtr = train_set.ms.x;
ytr = train_set.ms.y;
spdt = train_set.ms.speed;

test_set = load(testday);
test_spk = test_set.spike_mat(segsday2,:);

conv_test_spk = nan(size(test_spk));
for seg = 1:size(test_spk,1)
    conv_test_spk(seg, :) = conv(test_spk(seg,:), conv_kern, 'same');
    conv_test_spk(seg, conv_test_spk(seg,:)<1) = 0;
    conv_test_spk(seg, conv_test_spk(seg,:)>=1) = 1;
end

if strcmp(trainday, testday)
    numparts = 10; % ten fold validation
else
    numparts = 1; % ten fold validation
end
train_partition_ID = Partition_session(numparts, spdt, params.spd_thresh, 1*30/4);

  
xte = test_set.ms.x;
yte = test_set.ms.y;
spdte = test_set.ms.speed;
test_partition_ID = Partition_session(numparts, spdte, params.spd_thresh, 1*30/4);


predictions = NaN(length(test_partition_ID), params.decoding_nbins^2);
Z = NaN(params.decoding_nbins, params.decoding_nbins);

if rand_flag
    conv_rand_spk = conv_train_spk;
    for seg = 1:size(test_spk,1)
        conv_rand_spk(seg, :) = Shuffle_ISI(conv_rand_spk(seg, :));
    end
    rand_predictions = NaN(length(train_partition_ID), params.decoding_nbins^2);
end

for part = 1:numparts
    if numparts > 1
        trainspikes = conv_train_spk(:, [train_partition_ID ~= part & train_partition_ID > 0]);
        testspikes = conv_test_spk(:, test_partition_ID == part);
        x = xtr(train_partition_ID ~= part & train_partition_ID > 0);
        y = ytr(train_partition_ID ~= part & train_partition_ID > 0);
    else
        trainspikes = conv_train_spk(:, train_partition_ID == part);
        testspikes = conv_test_spk(:, test_partition_ID == part);
        x = xtr(train_partition_ID > 0);
        y = ytr(train_partition_ID > 0);
    end
    
    vmap_train = make_vmap_2D(x, y, params.decoding_bin_edges, params.decoding_nbins, 1);
    [~, ~, ~, xtrainbin, ytrainbin] = histcounts2(x, y, params.decoding_bin_edges, params.decoding_bin_edges);
    binID = sub2ind(size(Z), xtrainbin, ytrainbin);

    [trained_model] = train_SVM_decoder(vmap_train, trainspikes, binID);
    [pred] = test_SVM_decoder(trained_model, testspikes);
    predictions(test_partition_ID == part, :) = pred;
    if rand_flag    
        [rand_trained_model] = train_SVM_decoder(vmap_train, conv_rand_spk, binID);
        [rand_pred] = test_SVM_decoder(rand_trained_model, testspikes);
        rand_predictions(test_partition_ID == part, :) = rand_pred;
    end  
end    

%%
xte(test_partition_ID == 0) = NaN;
yte(test_partition_ID == 0) = NaN;

vmap_test = make_vmap_2D(xte, yte, params.decoding_bin_edges, params.decoding_nbins, 1);
[~, ~, ~, xtestbin, ytestbin] = histcounts2(xte, yte, params.decoding_bin_edges, params.decoding_bin_edges);
xtestbin(isnan(xte)) = NaN;
ytestbin(isnan(yte)) = NaN;
predicted_x = NaN(size(conv_test_spk,2), 1);
predicted_y = NaN(size(conv_test_spk,2), 1);
% ubin = unique(binID);
% ubin = ubin(~isnan(ubin));
for t = 1:size(conv_test_spk, 1)
    temp = vmap_test*NaN;
    temp(:) = predictions(t, :);
    [px, py] = find(max(temp(:)) == temp);
    predicted_x(t) = mean(px);
    predicted_y(t) = mean(py);
end
full.predx = predicted_x;
full.predy = predicted_y;
full.distances = sqrt((predicted_x - xtestbin).^2 + (predicted_y - ytestbin).^2);
full.dist_med = nanmedian(full.distances);
full.dist_std = std(full.distances(~isnan(full.distances)));
% distances when convolved spike occurred
full.has_spike = logical(sum(conv_test_spk,1));
full.good_distances = full.distances(full.has_spike);
full.good_dist_med = nanmedian(full.good_distances);
full.good_dist_std = std(full.good_distances(~isnan(full.good_distances)));
% distances when no convolved spike occurred
full.bad_distances = full.distances(~full.has_spike);
full.bad_dist_med = nanmedian(full.bad_distances);
full.bad_dist_std = std(full.bad_distances(~isnan(full.bad_distances)));

if rand_flag    
    rand_predicted_x = NaN(size(conv_test_spk,2), 1);
    rand_predicted_y = NaN(size(conv_test_spk,2), 1);
    for t = 1:size(conv_test_spk, 1)
        rtemp = vmap_test*NaN;
        rtemp(:) = rand_predictions(t, :);
        [rx, ry] = find(max(rtemp(:)) == rtemp);
        rand_predicted_x(t) = mean(rx);
        rand_predicted_y(t) = mean(ry);
    end
    random.predx = rand_predicted_x;
    random.predy = rand_predicted_y;

    random.distances = sqrt((rand_predicted_x - xtestbin).^2 + (rand_predicted_y - ytestbin).^2);
    random.dist_med = nanmedian(random.distances);
    random.dist_std = std(random.distances(~isnan(random.distances)));
    % distances when convolved spike occurred
    random.has_spike = logical(sum(conv_test_spk,1));
    random.good_distances = random.distances(random.has_spike);
    random.good_dist_med = nanmedian(random.good_distances);
    random.good_dist_std = std(random.good_distances(~isnan(random.good_distances)));
    % distances when no convolved spike occurred
    random.bad_distances = random.distances(~random.has_spike);
    random.bad_dist_med = nanmedian(random.bad_distances);
    random.bad_dist_std = std(random.bad_distances(~isnan(random.bad_distances)));
else
    random = NaN;
end
