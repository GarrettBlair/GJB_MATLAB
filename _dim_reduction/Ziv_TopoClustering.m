function [topo_out] = Ziv_TopoClustering(Reduced_data, position_active_frames)
%% Create barcodes plot and data
plotting = false;
data=Reduced_data;
cmap_position = jet(length(unique(position_active_frames)));
tic
cluster_size=450; %minimal number of data points per cluster
% cluster_size=200; % GJB minimal number of data points per cluster

distances = squareform(pdist(data));
distances(distances == 0) = 1; %
r_vec = 0.000015:0.00001:0.003; 

% all_bins_mat=nan(length(distances),length(r_vec));
% values_array={};%%% GJB preallocating
% indices_array={};%%% GJB preallocating
% occurrences_array={};%%% GJB preallocating
values_array=cell(length(r_vec),1);
indices_array=cell(length(r_vec),1);
occurrences_array=cell(length(r_vec),1);
% num_of_clusters=cell(length(r_vec),1);
for i = 2:length(r_vec)
    
    current_graph = distances < r_vec(i);
    
    bins = conncomp(graph(current_graph));
    all_bins_mat(:,i)=bins;
    [values, indices] = unique(bins);
    occurrences = hist(bins, values);
    values_array{i}=values;
    indices_array{i}=indices;
    occurrences_array{i}=occurrences;
    
    values = values(occurrences > cluster_size);
    indices = indices(occurrences > cluster_size); %% GJB not used
   
%     num_of_clusters(i)=length(values);
end

% Create fixed indexing arrays by keeping the cluster number of the larger
% cluster when merging

% new_ind_array={}; %%% GJB preallocating
new_ind_array=cell(length(r_vec),1);
new_ind_array{2} = indices_array{2};
for i = 3:length(r_vec)
    current_indices = indices_array{i};
    current_bins = all_bins_mat(:, i);
    prev_indices = indices_array{i - 1};
    prev_bins = all_bins_mat(:, i - 1);
    prev_occurr = occurrences_array{i - 1};
    
    prev_new_ind = new_ind_array{i - 1};
    
    new_ind_vec = nan(1, length(current_indices));
    for run_cluster=1:length(current_indices) %%% GJB par
        sons_of_run_clusters = prev_new_ind(current_bins(prev_indices) == prev_bins(prev_indices(run_cluster)));
        occurr_of_sons = prev_occurr(current_bins(prev_indices) == prev_bins(prev_new_ind(run_cluster)));
        [~, most_freq_son_ind] = max(occurr_of_sons);

        new_ind_vec(run_cluster) = prev_new_ind(prev_new_ind == sons_of_run_clusters(most_freq_son_ind));
            
    end
    
    new_ind_array{i} = new_ind_vec; 
end

% Plot fixed indices barcodes

if plotting %%%%%%%%%%%%%%%
figure;
hold on;

for i = 2:length(r_vec)
    indices = new_ind_array{i};
    occurrences = occurrences_array{i};
    
    indices = indices(occurrences > cluster_size);
    
    plot(r_vec(i) * ones(length(indices), 1), indices, 'k*');
end

xlabel('Radius')
ylabel('H_0')
end %%%%%%%%%%%%%%%

TopoClusteringMat=[];
for i = 2:length(r_vec)
    indices = new_ind_array{i};
    occurrences = occurrences_array{i}; 
    indices = indices(occurrences > cluster_size);
    
    for run=1:length(indices)       
    TopoClusteringMat(i,indices(run))=1;
    end
end
toc
%%
active_clusters=find(sum(TopoClusteringMat,1));
figure;
hold on;
for run = 1:length(active_clusters)
    current_cluster_ind=active_clusters(run);
    current_cluster_r_vec=r_vec(~~TopoClusteringMat(:,current_cluster_ind));
    plot(min(current_cluster_r_vec),max(current_cluster_r_vec), 'bo');  
end
plot([0 max(r_vec)],[0 max(r_vec)],'k-')
xlim([0 max(r_vec)])
ylim([0 max(r_vec)])
axis square
xlabel('Birth')
ylabel('Death')
%%
% Now change to cluster_size=20
%%%%%%%%%%%%%%% GJB
cluster_size = 20;
%%%%%%%%%%%%
%%

segments_array=nan(length(active_clusters),2);
FatherClusterVec=nan(length(active_clusters),1);
for run = 1:length(active_clusters)
    current_cluster_ind=active_clusters(run);
    current_cluster_r_vec=r_vec(~~TopoClusteringMat(:,current_cluster_ind));
%     plot(current_cluster_r_vec,current_cluster_ind*ones(1,length(current_cluster_r_vec)), 'b-');
    last_frame_cluster_existence=max(find(~~TopoClusteringMat(:,current_cluster_ind)));
    
    IdentityClusterBeforeClusterDie=new_ind_array{last_frame_cluster_existence};
    
    if last_frame_cluster_existence~=length(new_ind_array)
    IdentityClusterAfterClusterDie=new_ind_array{last_frame_cluster_existence+1};
    PointsPerClusterBeforeClusterDie=occurrences_array{last_frame_cluster_existence};
    PointsPerClusterAfterClusterDie=occurrences_array{last_frame_cluster_existence+1};
    DiffPointsPerCluster=PointsPerClusterAfterClusterDie-PointsPerClusterBeforeClusterDie(ismember(IdentityClusterBeforeClusterDie,IdentityClusterAfterClusterDie));
    
    FatherClusterVec(run)=IdentityClusterAfterClusterDie(DiffPointsPerCluster==max(DiffPointsPerCluster));    
    segments_array(run,:)=[min(current_cluster_r_vec) max(current_cluster_r_vec)];
    else
        FatherClusterVec(run)=current_cluster_ind;  %cluster that survive is its own fater  
    segments_array(run,:)=[min(current_cluster_r_vec) max(current_cluster_r_vec)];
    end
end


FatherClusterVec2=nan(1,length(active_clusters));
for run=1:length(FatherClusterVec)
    FatherClusterVec2(run)=find(active_clusters==FatherClusterVec(run));
end
SonPerCluster={};
for run=1:length(FatherClusterVec)
    SonPerCluster{run}=find(FatherClusterVec2==(run));
end

loc=zeros(1,length(active_clusters));
[~,AncestorInd]=max(segments_array(:,2));
CurrInd=AncestorInd;
CurrLoc=1;
loc(CurrInd)=CurrLoc;
while 1
   %['CurrInd = ' num2str(CurrInd)]
    SonsList=SonPerCluster{CurrInd};
    if length(SonsList)>0 & ~(length(SonsList)==1 & SonsList==AncestorInd)
        [~,Ind]=min(segments_array(SonsList,2));
        CurrInd=SonsList(Ind);
        CurrLoc=CurrLoc+1;
        loc(CurrInd)=CurrLoc;
    else
        if CurrInd==AncestorInd
            break
        else
            PrevInd=CurrInd;
            CurrInd=FatherClusterVec2(CurrInd);
            CurrentSonsList=SonPerCluster{CurrInd};
            SonPerCluster{CurrInd}=CurrentSonsList(CurrentSonsList~=PrevInd);
        end
    end 
end
if plotting %%%%%%%%%%%%%%%
figure
hold on
for run = 1:length(active_clusters)  
    plot(segments_array(run,:),loc(run)*[1 1],'k-')
    plot([1 1]*segments_array(run,2),[loc(run) loc(FatherClusterVec2(run))],'k-')
end  
xlabel('Radius')
ylabel('H_0')
end %%%%%%%%%%%%%%%
%% Plot clusters by r
R = 0.0006350; %0.2; %0.05; %
% R = 0.2; %0.05; %
% R = 0.05; %

if plotting %%%%%%%%%%%%%%%
figure;
hold on;

scatter(data(:, 1), data(:, 2), 5, '.');
% scatter(data(:, 2), data(:, 3), 5, '.');%% GJB
end %%%%%%%%%%%%%%%


[~,current_ind]=min(abs(r_vec - R));

indices = new_ind_array{current_ind};
occurrences = occurrences_array{current_ind};
    
indices = indices(occurrences > cluster_size);
NUMBER_OF_CLUSTERS = length(indices);
color_map=hsv(NUMBER_OF_CLUSTERS);

ind = ones(size(data, 1), 1) * (NUMBER_OF_CLUSTERS + 1);

non_chosen_points=ones(1,size(data,1));
for i = 1:length(indices)  
    chosen_points=~~(all_bins_mat(:,current_ind)==all_bins_mat(indices(i),current_ind));
%     scatter(data(chosen_points, 2), data(chosen_points, 3), 20, color_map(i,:), '.');
if plotting %%%%%%%%%%%%%%%
    scatter(data(chosen_points, 1), data(chosen_points, 2), 20, color_map(i,:), '.');
end %%%%%%%%%%%%%%%
    non_chosen_points(chosen_points)=0;
    ind(chosen_points) = i;
end
if plotting %%%%%%%%%%%%%%%
title(num2str(R));
end %%%%%%%%%%%%%%%

%%

if plotting %%%%%%%%%%%%%%%
figure
NumOfSubPlotsCol=ceil(sqrt(NUMBER_OF_CLUSTERS+1));
NumOfSubPlotsRow=ceil((NUMBER_OF_CLUSTERS+1)/NumOfSubPlotsCol);
for run=1:NUMBER_OF_CLUSTERS+1
    subplot(NumOfSubPlotsRow,NumOfSubPlotsCol,run)
    plot(data(:, 1), data(:, 2),'k.','MarkerSize',1)
    hold on
    plot(data(ind==run, 1), data(ind==run, 2),'ro','markerfacecolor', 'r','markersize',1);
    xticks([])
    yticks([])
end
end %%%%%%%%%%%%%%%
%%
position_of_chosen_data=position_active_frames;
num_of_survived_clusters=NUMBER_OF_CLUSTERS;
last_ind_per_index=current_ind*ones(1,NUMBER_OF_CLUSTERS);
current_ind2=current_ind;
disp('A')

%%%%%%%%%%%%%% GJB if statement
while num_of_survived_clusters>2
    num_of_survived_clusters
    current_ind2=current_ind2+1;
    if current_ind2 <= length(new_ind_array)
    current_ind_vec=new_ind_array{current_ind2};
    survival_of_clusters=ismember(indices,current_ind_vec);
    last_ind_per_index=last_ind_per_index+survival_of_clusters;
    num_of_survived_clusters=sum(survival_of_clusters);
    else
        break
    end
end
%%%%%%%%%%%%%% GJB
disp('B')
if num_of_survived_clusters<2
    last_ind_per_index=last_ind_per_index-survival_of_clusters;
end

[index_of_cluster_ordered,order_of_ind]=sort(last_ind_per_index,'descend')
sorted_ind_vec=indices(order_of_ind);

ind_augmented=(NUMBER_OF_CLUSTERS+1)*ones(size(ind));
for run=1:length(indices)
    relevat_index_for_the_cluster=index_of_cluster_ordered(run);
    current_cluster_vector=(all_bins_mat(:,relevat_index_for_the_cluster)==all_bins_mat(sorted_ind_vec(run),relevat_index_for_the_cluster));
    ind_augmented(current_cluster_vector)=order_of_ind(run);
end

if plotting %%%%%%%%%%%%%%%
figure
for run=1:NUMBER_OF_CLUSTERS
    hold on
    scatter3(data(ind_augmented==run, 1), data(ind_augmented==run, 2), data(ind_augmented==run, 3), 20,  color_map(run,:), '.');
end

figure
NumOfSubPlotsCol=ceil(sqrt(NUMBER_OF_CLUSTERS+1));
NumOfSubPlotsRow=ceil((NUMBER_OF_CLUSTERS+1)/NumOfSubPlotsCol);
for run=1:NUMBER_OF_CLUSTERS+1
    subplot(NumOfSubPlotsRow,NumOfSubPlotsCol,run)
hist(position_of_chosen_data(ind_augmented==run),[1:24])
end

figure
NumOfSubPlotsCol=ceil(sqrt(NUMBER_OF_CLUSTERS+1));
NumOfSubPlotsRow=ceil((NUMBER_OF_CLUSTERS+1)/NumOfSubPlotsCol);
for run=1:NUMBER_OF_CLUSTERS+1
    subplot(NumOfSubPlotsRow,NumOfSubPlotsCol,run)
    plot(data(:, 1), data(:, 2),'k.','MarkerSize',1)
    hold on
    %scatter3(data(ind==run, 1), data(ind==run, 2), data(ind==run, 3), 20, cmap_position(position_of_chosen_data(ind==run),:), '.');
    scatter(data(ind_augmented==run, 1), data(ind_augmented==run, 2), 20, cmap_position(position_of_chosen_data(ind_augmented==run),:), '.');
    %scatter3(data(ind==run, 1), data(ind==run, 2), data(ind==run, 3), 20, cmap_velocity(velocity_ind(ind==run),:), '.');
end
end %%%%%%%%%%%%%%%

%%

ind_augmented2=ind;
ind_augmented2(ind==NUMBER_OF_CLUSTERS+1)=nan;
for run_R=current_ind:length(r_vec)
    current_R_threshold=r_vec(run_R);
    labeled_ind=ind_augmented2(~isnan(ind_augmented2));
    unlabeled_ind=ind_augmented2(isnan(ind_augmented2));
%     distance_labeled_to_unlabeled=d2(isnan(ind_augmented2),~isnan(ind_augmented2));
    distance_labeled_to_unlabeled=distances(isnan(ind_augmented2),~isnan(ind_augmented2));
    [min_value_labeling entry_of_min_value_labeling]=min(distance_labeled_to_unlabeled,[],2);
    
    should_be_labeled=(min_value_labeling<current_R_threshold);
    unlabeled_ind(should_be_labeled)=labeled_ind(entry_of_min_value_labeling(should_be_labeled));
    ind_augmented2(isnan(ind_augmented2))=unlabeled_ind;
end

ind_augmented2(isnan(ind_augmented2))=NUMBER_OF_CLUSTERS+1;

if plotting %%%%%%%%%%%%%%%
figure
for run=1:NUMBER_OF_CLUSTERS
    hold on
    scatter3(data(ind_augmented2==run, 1), data(ind_augmented2==run, 2), data(ind_augmented2==run, 3), 20,  color_map(run,:), '.');
end

figure
NumOfSubPlotsCol=ceil(sqrt(NUMBER_OF_CLUSTERS+1));
NumOfSubPlotsRow=ceil((NUMBER_OF_CLUSTERS+1)/NumOfSubPlotsCol);
for run=1:NUMBER_OF_CLUSTERS+1
    subplot(NumOfSubPlotsRow,NumOfSubPlotsCol,run)
hist(position_of_chosen_data(ind_augmented2==run),[1:24])
end

figure
NumOfSubPlotsCol=ceil(sqrt(NUMBER_OF_CLUSTERS+1));
NumOfSubPlotsRow=ceil((NUMBER_OF_CLUSTERS+1)/NumOfSubPlotsCol);
for run=1:NUMBER_OF_CLUSTERS+1
    subplot(NumOfSubPlotsRow,NumOfSubPlotsCol,run)
    plot(data(:, 1), data(:, 2),'k.','MarkerSize',1)
    hold on
    %scatter3(data(ind==run, 1), data(ind==run, 2), data(ind==run, 3), 20, cmap_position(position_of_chosen_data(ind==run),:), '.');
    scatter(data(ind_augmented2==run, 1), data(ind_augmented2==run, 2), 20, cmap_position(position_of_chosen_data(ind_augmented2==run),:), '.');
    %scatter3(data(ind==run, 1), data(ind==run, 2), data(ind==run, 3), 20, cmap_velocity(velocity_ind(ind==run),:), '.');
end

figure
NumOfSubPlotsCol=ceil(sqrt(NUMBER_OF_CLUSTERS+1));
NumOfSubPlotsRow=ceil((NUMBER_OF_CLUSTERS+1)/NumOfSubPlotsCol);
for run=1:NUMBER_OF_CLUSTERS+1
    subplot(NumOfSubPlotsRow,NumOfSubPlotsCol,run)
    plot(data(:, 1), data(:, 2),'k.','MarkerSize',1)
    hold on
    %scatter3(data(ind==run, 1), data(ind==run, 2), data(ind==run, 3), 20, cmap_position(position_of_chosen_data(ind==run),:), '.');
    %scatter(data(ind_augmented2==run, 1), data(ind_augmented2==run, 2), 20, cmap_position(position_of_chosen_data(ind_augmented2==run),:), '.');
    %scatter3(data(ind==run, 1), data(ind==run, 2), data(ind==run, 3), 20, cmap_velocity(velocity_ind(ind==run),:), '.');
    scatter(data(ind_augmented2==run, 1), data(ind_augmented2==run, 2), 20,cmap_velocity(velocity_ind(ind_augmented2==run),:), '.');
end
end %%%%%%%%%%%%%%%
%%
%%%%%% GJB commented out if state
% % % % % % % if 1 %labeling unlabled data points in one shut -based on label of nearest labled point
    unlabeled_data=data(ind_augmented2==NUMBER_OF_CLUSTERS+1,:);
    guess_of_umlabeled=entry_of_min_value_labeling(~should_be_labeled);
if plotting %%%%%%%%%%%%%%%
    figure;
    plot3(data(:, 1), data(:, 2), data(:, 3),'k.','MarkerSize',1);
    hold on; 
    scatter3(unlabeled_data(:, 1), unlabeled_data(:, 2), unlabeled_data(:, 3),5,color_map(labeled_ind(guess_of_umlabeled),:))
end %%%%%%%%%%%%%%%
    
    ind_augmented3=ind_augmented2;
    ind_augmented3(ind_augmented2==NUMBER_OF_CLUSTERS+1)=labeled_ind(guess_of_umlabeled);
    %TopoClustering.global.ind_augmented3=ind_augmented3;
    
if plotting %%%%%%%%%%%%%%%
    figure
    for run=1:NUMBER_OF_CLUSTERS
        hold on
        scatter3(data(ind_augmented3==run, 1), data(ind_augmented3==run, 2), data(ind_augmented3==run, 3), 20,  color_map(run,:), '.');
    end
    
    figure
    NumOfSubPlotsCol=ceil(sqrt(NUMBER_OF_CLUSTERS+1));
    NumOfSubPlotsRow=ceil((NUMBER_OF_CLUSTERS+1)/NumOfSubPlotsCol);
    for run=1:NUMBER_OF_CLUSTERS+1
        subplot(NumOfSubPlotsRow,NumOfSubPlotsCol,run)
        hist(position_of_chosen_data(ind_augmented3==run),[1:24])
    end
    
    figure
    NumOfSubPlotsCol=ceil(sqrt(NUMBER_OF_CLUSTERS+1));
    NumOfSubPlotsRow=ceil((NUMBER_OF_CLUSTERS+1)/NumOfSubPlotsCol);
    for run=1:NUMBER_OF_CLUSTERS+1
        subplot(NumOfSubPlotsRow,NumOfSubPlotsCol,run)
        plot(data(:, 1), data(:, 2),'k.','MarkerSize',1)
        hold on
        %scatter3(data(ind==run, 1), data(ind==run, 2), data(ind==run, 3), 20, cmap_position(position_of_chosen_data(ind==run),:), '.');
        scatter(data(ind_augmented3==run, 1), data(ind_augmented3==run, 2), 20, cmap_position(position_of_chosen_data(ind_augmented3==run),:), '.');
        %scatter3(data(ind==run, 1), data(ind==run, 2), data(ind==run, 3), 20, cmap_velocity(velocity_ind(ind==run),:), '.');
    end
end %%%%%%%%%%%%%%%
    
    
%     figure
%     NumOfSubPlotsCol=ceil(sqrt(NUMBER_OF_CLUSTERS+1));
%     NumOfSubPlotsRow=ceil((NUMBER_OF_CLUSTERS+1)/NumOfSubPlotsCol);
%     for run=1:NUMBER_OF_CLUSTERS+1
%         subplot(NumOfSubPlotsRow,NumOfSubPlotsCol,run)
%         plot(data(:, 1), data(:, 2),'k.','MarkerSize',1)
%         hold on
%         %scatter3(data(ind==run, 1), data(ind==run, 2), data(ind==run, 3), 20, cmap_position(position_of_chosen_data(ind==run),:), '.');
%         %scatter(data(ind_augmented3==run, 1), data(ind_augmented3==run, 2), 20, cmap_position(position_of_chosen_data(ind_augmented3==run),:), '.');
%         %scatter3(data(ind==run, 1), data(ind==run, 2), data(ind==run, 3), 20, cmap_velocity(velocity_ind(ind==run),:), '.');
%         scatter(data(ind_augmented3==run, 1), data(ind_augmented3==run, 2), 20,cmap_velocity(velocity_ind(ind_augmented3==run),:), '.');
%     end
% % % % % % % % end

%%%%%%%% GJB ADDED
topo_out.global.ind = ind;
topo_out.global.ind_augmented = ind_augmented;
topo_out.global.ind_augmented2 = ind_augmented2;
topo_out.global.ind_augmented3 = ind_augmented3;