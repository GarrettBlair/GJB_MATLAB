function [models] = gb_conjoint_ensemble_model_fitting(e, thresh_e, edges_of_bins, plotting)
% clear
% load('D:\Sample Data\ensemble_prob_linear\standard\256_0.mat')
% load('D:\Sample Data\ensemble_prob_linear\standard\256_0.mat')
% load('D:\Sample Data\ensemble_prob_linear\standard\256_3.mat')
% load("D:\Sample Data\ensemble_prob_median\prop_samples\100_0.mat")
%
if plotting==true
figure; clf;
hold on; 
end
% bin_size = .001;
% 
% % e = 2.^ensem1';
% e = normalize_matrix(ensem1)+bin_size;
% thresh_e = median(e)+2*std(e);
% 
% edges_of_bins=[bin_size:bin_size:1.2];
options = statset('MaxIter',10000, 'MaxFunEvals',20000, 'Display', 'off');

cents = edges_of_bins(1:end-1) + median(abs(diff(edges_of_bins)));
centers_of_bins{1}=cents;

microns_per_pixel = 1;
data = e;

number_of_bins=length(centers_of_bins{1});
bin_center=centers_of_bins{1};
data_span = (bin_center(end)-bin_center(1));

[data_dist,~]=histcounts(data,edges_of_bins);
data_dist=data_dist./sum(data_dist)*(number_of_bins/(microns_per_pixel*(bin_center(2)-bin_center(1))+microns_per_pixel*(bin_center(end)-bin_center(1))));

% finding initial parameters for lsqcurvefit function:
maximal_distance_to_fit=thresh_e; % this distance is used only to find the initial paramters
data_to_fit=data(data<maximal_distance_to_fit);
parmhat=lognfit(data_to_fit);
optimal_delta=length(data_to_fit)/length(data);
p_0=optimal_delta;
c_0=6;
a_0=1;
% b_0=(data_dist(end)-data_dist(round(number_of_bins/2)))/((bin_center(end)-bin_center(round(number_of_bins/2))));
b_0=(data_to_fit(end)-data_to_fit(round(number_of_bins/2)))/((bin_center(end)-bin_center(round(number_of_bins/2))));
b_0=b_0/(1-p_0);
initial_parameters=[p_0 parmhat a_0 c_0 b_0];
F = @(x,xdata)...
    x(1)*(1./(xdata.*x(3).*sqrt(2*pi)).*exp(-(log(xdata)-x(2)).^2./(2*x(3)^2))...
    + (1-x(1)).*x(6).*xdata./(1+exp(-x(4).*(xdata-x(5)))));

lb = [0 -Inf 0   0   0   0];
ub = [1 Inf  Inf Inf Inf Inf];
% options = statset('MaxIter',1000, 'MaxFunEvals',2000, 'Display', 'off');

% finding the parameters that best fit the data:
% centroid_distances_model_parameters=lsqcurvefit(F,initial_parameters,bin_center,data_dist,lb,ub,options);
data_sub = data_dist(bin_center<(thresh_e));
bin_sub = bin_center(bin_center<(thresh_e));
centroid_distances_model_parameters=lsqcurvefit(F,initial_parameters,bin_sub,data_sub,lb,ub,options);

% calculating the distribution for same cells:
model_lower=lognpdf(bin_center,centroid_distances_model_parameters(2),centroid_distances_model_parameters(3));
model_lower=model_lower./sum(model_lower)*(number_of_bins/((bin_center(2)-bin_center(1))+microns_per_pixel*(bin_center(end)-bin_center(1))));


maximal_distance_to_fit = thresh_e; % this distance is used only to find the initial paramters
data_to_fit=data(data>=maximal_distance_to_fit);
parmhat=lognfit(data_to_fit);
optimal_delta=1-length(data_to_fit)/length(data);
p_0=optimal_delta;
c_0=6;
a_0=1;
b_0=(data_dist(end)-data_dist(round(number_of_bins/2)))/((bin_center(end)-bin_center(round(number_of_bins/2))));
b_0=b_0/(p_0);
initial_parameters=[p_0 parmhat a_0 c_0 b_0];
F = @(x,xdata)...
    x(1)*(1./(xdata.*x(3).*sqrt(2*pi)).*exp(-(log(xdata)-x(2)).^2./(2*x(3)^2))...
    + (1-x(1)).*x(6).*xdata./(1+exp(-x(4).*(xdata-x(5)))));

lb = [0 -Inf 0   0   0   0];
ub = [1 Inf  Inf Inf Inf Inf];
% options = statset('MaxIter',1000, 'MaxFunEvals',2000, 'Display', 'off');

% finding the parameters that best fit the data:
% centroid_distances_model_parameters=lsqcurvefit(F,initial_parameters,bin_center,data_dist,lb,ub,options);
data_sub = data_dist(bin_center>=(thresh_e));
bin_sub = bin_center(bin_center>=(thresh_e));
centroid_distances_model_parameters=lsqcurvefit(F,initial_parameters,bin_sub,data_sub,lb,ub,options);

% calculating the distribution for same cells:
model_upper=lognpdf(bin_center,centroid_distances_model_parameters(2),centroid_distances_model_parameters(3));
model_upper=model_upper./sum(model_upper)*(number_of_bins/((bin_center(2)-bin_center(1))+microns_per_pixel*(bin_center(end)-bin_center(1))));

model_intersect = model_upper.*model_lower;%./(sum(model_resid.*model_out));
model_product = model_upper*model_lower';

indep_sum = model_upper*model_upper' + model_lower*model_lower';

if plotting==true
subplot(3,1,1); hold on
plot(data); axis tight
subplot(3,1,2:3); hold on
plot(bin_center, data_dist, 'k')
plot(bin_center, model_lower, 'b')
plot(bin_center, model_upper, 'r-')
yyaxis('right')
plot(bin_center, model_intersect, 'm:')
plot([thresh_e thresh_e], max([model_lower model_upper model_intersect]), 'k')

% subplot(3,1,3); hold on
% plot(bin_center, cumsum(model_lower), 'b')
% plot(bin_center, cumsum(model_upper), 'r')
% 
% yyaxis('right')
% plot(bin_center, cumsum(model_intersect), 'm:')
title(sprintf('overlap = %1.3f', model_product/indep_sum))
end

models = [];
models.data      = data;
models.split_val = thresh_e;
models.bin_center= cents;
models.lower     = model_lower;
models.upper     = model_upper;
models.intersect = model_intersect;
models.product   = model_product;
models.indep_sum = indep_sum;










