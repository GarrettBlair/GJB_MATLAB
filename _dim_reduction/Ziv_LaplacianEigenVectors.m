function [ReducedDataArray, legitimacy_vec, data_out] = Ziv_LaplacianEigenVectors(neuronal_activity_mat, plotting, p_neighbors_vec)
%% Adapted by GJB from [Rubin,... Ziv et al. (2020) NatComm] methods paper
%% This script allows performing analysis on hippocampal data. 

ActivityThreshold=1; % orig=1
% ActivityThreshold=3; % orig=1
% ActivityThreshold=ceil(size(neuronal_activity_mat,1)*.1); % GJB
% % fprintf('Activity threshold: %d, NumSegs = %d\n', ActivityThreshold, size(neuronal_activity_mat,1))
% downsample_flag=0;%%GJB
downsample_flag=0;%% orig=2, conv method; 0=no ds; 1=slicing ds
% downsample_flag=1;%% GJB
ActvitySmoothingSize=0;%% orig=2
% ActvitySmoothingSize=0;%%GJB

if nargin<3
    p_neighbors_vec=[0.075/15 0.075];%% orig
end
% p_neighbors_vec=[0.075/15 0.075];%% GJB
% p_neighbors_vec=[0.075 0.075/15];%% GJB
% p_neighbors_vec=[4 15] %GJB 19AQ
% p_neighbors_vec=[.25/15 .25] %GJB

ChosenDirectionArray={'Positive' 'Negative'};

global cmap_position
cmap_position=jet(24);

global SubTypeColorMapArrayByNumer
SubTypeColorMapArrayByNumer{1}=[0 0 0];
SubTypeColorMapArrayByNumer{2}=[1 0 0; 0 0 1];
SubTypeColorMapArrayByNumer{3}=[1 0 0; 0 1 0; 0 0 1];
SubTypeColorMapArrayByNumer{4}=[1 0 0; 0 1 0; 0 0 1; 1 0 1];

%%GJB for checking random data
% random_act = NaN*neuronal_activity_mat;
% for j = 1:size(neuronal_activity_mat, 1)
%     a = neuronal_activity_mat(j,:);
%     [~, ord] = sort(rand(size(neuronal_activity_mat,2), 1));
%     random_act(j,:) = a(ord);
% end
%%GJB for checking random data
position_per_frame = zeros(size(neuronal_activity_mat, 2), 1);
[all_trials_activity_smoothed_mat, position_active_frames, activity_mat_active_frames, legitimacy_vec]=...
    SmoothingAndThresholdingData(position_per_frame,neuronal_activity_mat,downsample_flag,ActvitySmoothingSize,ActivityThreshold);
%% Dimentionality Reduction
rng('default')
activity_mat_active_frames = neuronal_activity_mat;
ReducedDataArray=DimentionalityReduction_Ver1(activity_mat_active_frames',p_neighbors_vec);
data_out = ReducedDataArray{3};
% data_out = NaN(size(neuronal_activity_mat,2), size(ReducedDataArray{3},2));
% data_out(legitimacy_vec, :) = ReducedDataArray{3};
if plotting
    
v2=ReducedDataArray{3};
%%
%
max_v2_2=max(v2(:,2));
min_v2_2=min(v2(:,2));
normed_v2_2=24*(v2(:,2)-min_v2_2)/(max_v2_2-min_v2_2);
% if corr(normed_v2_2,position_active_frames')<0
%     normed_v2_2=25-normed_v2_2;
% end

figure, 
plot3(v2(:,2),v2(:,3),v2(:,4),'k.','markersize',5)
hold on
view([20,20])
xlabel('Comp. 1')
ylabel('Comp. 2')
zlabel('Comp. 3')
% Panel1A_lim=0.025;
% xlim([-Panel1A_lim Panel1A_lim])
% ylim([-Panel1A_lim Panel1A_lim])
% zlim([-Panel1A_lim Panel1A_lim])
set(gca,'XTick',[])
set(gca,'XTickLabel',[]);
set(gca,'YTick',[])
set(gca,'YTickLabel',[]);
set(gca,'ZTick',[])
set(gca,'ZTickLabel',[]);
axis square
%%
% nsurr = 10;
% for ii = 101:length(v2(:,2))-100
% cla
% plot3(v2(:,2),v2(:,3),v2(:,4),'k.','markersize',5)
% plot3(v2(ii-nsurr:ii,2),v2(ii-nsurr:ii,3),v2(ii-nsurr:ii,4),'r-','markersize',5)
% plot3(v2(ii:ii+nsurr,2),v2(ii:ii+nsurr,3),v2(ii:ii+nsurr,4),'g-','markersize',5)
% plot3(v2(ii,2),v2(ii,3),v2(ii,4),'m.','markersize',30)
% drawnow
% end
end