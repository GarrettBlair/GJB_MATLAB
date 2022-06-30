function [comp] = test_SVM_decoder(mdl, spks)
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
comp = zeros(size(spks, 2), size(mdl, 1));
parfor t = 1:size(spks, 2)
    s_vec = spks(:,t);
    c = comp(t,:);
    for i = 1:size(mdl, 1) - 1
        for j = i+1:size(mdl, 2)
            if ~isempty(mdl{i,j})
                [pred, ~, ~] = predict(mdl{i,j}, s_vec');
                switch pred
                    case 1 % positive response == row bin
                        c(i) = c(i) + 1;
%                         c(j) = c(j) - 1;
                    case -1 % negative response == column bin
%                         c(i) = c(i) - 1;
                        c(j) = c(j) + 1;
                end
            end
        end
    end
    comp(t,:) = c;
end
