function [prob_map] = Bayesian_decoder(pfield, vmap, seg_prob, response_matrix, x, y)
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
%
% -Garrett Blair 9/6/2018

figure(3); clf

prob_map = NaN(size(vmap,1), size(vmap,2), size(response_matrix,2));
flag = .05;
plotting = true;
stopping = false;
for sample = 4:size(response_matrix, 2) % Step through every sample
    p = NaN(size(vmap,1), size(vmap,2), size(response_matrix,1));
    ind = 1;
    for segment = 1:size(response_matrix, 1)
        switch response_matrix(segment, sample)
            case 1
                cpr = squeeze(pfield(:, :, segment));
                pr  = seg_prob(segment);
                if plotting && ind<16
                    subplot(4,4,ind); cla
                    imagesc(cpr)
                    hold on;
                    plot(x(sample), y(sample),'g*')
                    hold off
                    axis square off
                    colorbar
                    title(sprintf('pr = %2.4f', pr))
                    ind = ind+1; 
                    stopping = true;
                end
            case 0
                cpr = 1 - squeeze(pfield(:, :, segment));
                pr  = 1 - seg_prob(segment);
        end
        p(:,:,segment) = (cpr.*vmap)./pr;
    end
    prob = nansum((p),3);
    if plotting
    prob_map(:,:,sample) = prob/nansum(prob(:));
    subplot(4,4,3*3)
    [centx, centy] = find(max(prob(:)) == prob);
    centx(sample) = mean(centx);
    centy(sample) = mean(centy);
    imagesc(prob_map(:,:,sample))
    axis square off
    hold on
    plot(centy(sample), centx(sample), 'r*')
    plot(x(sample), y(sample), 'g*')
    drawnow
    pause(.1)
    if stopping
        pause(.5)
        subplot(3,3,1:8)
        cla 
        stopping = false;
    end
%     stopping = false;
    end
    if sample/(size(response_matrix, 2)*flag)>1
        fprintf('..%d%%', round(flag*100))
        flag = flag + .05;
    end
end

fprintf('...Done!\n')
