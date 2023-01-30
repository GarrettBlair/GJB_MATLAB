function [mdl] = train_SVM_decoder(weight_vector, response_matrix, binID)
% Make a conditional probability map for every time step (columns of
% response_matrix).
%
% Inputs: 
%   pfield   - conditional probability of firing given location (place field)
%   vmap     - stimulus probability; occupancy probability
%   seg_prob - response probability; firing probability
%   response_matrix     - response matrix (neuronID x logical firing trace)
%
% Outputs:
%   prob_map - conditional probability map
%                                       
% based on the method in Stefanini et al 2020
% -Garrett Blair 9/6/2018

%%%%%%%%%%%%%%% USE Mdl = fitcecoc(response_matrix(:, train_inds)', bin_ang(train_inds)', 'Weights', wvec(train_inds));
%%%%%%%%%%%%%%% USE pred = predict(Mdl, response_matrix(:, test_inds)');

% vmapdims = length(vmap)^2;
ubin = unique(binID);
ubin = ubin(~isnan(ubin));
[nr, nc] = size(weight_vector);
vmapdims = nr*nc;
mdl = cell(vmapdims);

for i = 1:length(ubin)-1
    for j = i+1:length(ubin)
        c1 = ubin(i);
        act1 = response_matrix(:, c1 == binID);
        ans1 = ones(size(act1,2), 1);
        % positive response == row bin
%         weight1 = 1./( vmap(c1)*ones(size(act1,2), 1) );
        weight1 = weight_vector(c1)*ones(size(act1,2), 1);
        c2 = ubin(j);
        act2 = response_matrix(:, c2 == binID);
        ans2 = ones(size(act2,2), 1)*-1;
        % negative response == column bin
%         weight2 = 1./( vmap(c2)*ones(size(act2,2), 1) );
        weight2 = weight_vector(c2)*ones(size(act2,2), 1);
        
        mdl{ubin(i), ubin(j)} = fitcsvm([act1 act2]', [ans1; ans2], 'Weights', [weight1; weight2]);
    end
end
