function [neighbors_spatial_correlations,neighbors_centroid_distances]=GB_neighbor_dist_corrs(maximal_distance,spatial_footprints,centroid_locations)
% This function performs an initial cell registration across sessios
% based on a chosen spatial correlation thresohld. 

% Inputs:
% 1. maximal_distance
% 2. spatial_correlation_threshold
% 3. spatial_footprints
% 4. centroid_locations

% Outputs:
% 1. cell_to_index_map - list of registered cells and their index in each session
% 2. registered_cells_spatial_correlations
% 3. non_registered_cells_spatial_correlations

number_of_sessions=size(spatial_footprints,2);

% initializing the registration with the cells from session #1:
initial_number_of_cells=size(spatial_footprints{1},1);
registered_spatial_footprints=spatial_footprints{1};
registered_centroid_locations=centroid_locations{1};

% allocating space:
neighbors_spatial_correlations=zeros(1,number_of_sessions^2*initial_number_of_cells);
neighbors_centroid_distances=zeros(1,number_of_sessions^2*initial_number_of_cells);

neighbor_count=0;
num_candidates=0;
disp('Registering cells:');
disp('Initializing list with the cells from session #1');
display_progress_bar('Terminating previous progress bars',true)    
for n=2:number_of_sessions % registering the rest of the sessions
    display_progress_bar(['Registering cells in session #' num2str(n) ' - '],false)
    new_spatial_footprints=spatial_footprints{n};
    new_centroids=centroid_locations{n};    
    for k=1:size(new_spatial_footprints,1) % for each cell
        display_progress_bar(100*(k)/size(new_spatial_footprints,1),false)
        new_spatial_footprint=squeeze(new_spatial_footprints(k,:,:));
        centroid=repmat(new_centroids(k,:),size(registered_centroid_locations,1),1);
        distance_vec=sqrt(sum((centroid-registered_centroid_locations).^2,2));
        spatial_footprints_to_check=find(distance_vec<maximal_distance);
        if ~isempty(spatial_footprints_to_check) % finding the best candidate for each cell
            corr_vec=zeros(1,length(spatial_footprints_to_check));
            num_candidates=num_candidates+length(spatial_footprints_to_check);
            for m=1:length(spatial_footprints_to_check)
                suspected_spatial_footprint=squeeze(registered_spatial_footprints(spatial_footprints_to_check(m),:,:));
                if sum(sum(suspected_spatial_footprint))==0 || sum(sum(new_spatial_footprint))==0
                    corr_vec(m)=0;
                else
                    neighbor_count=neighbor_count+1;                  
                    corr_vec(m)=corr2(suspected_spatial_footprint,new_spatial_footprint);
                    neighbors_spatial_correlations(neighbor_count)=corr_vec(m);
                    neighbors_centroid_distances(neighbor_count)=distance_vec(spatial_footprints_to_check(m));
                end
            end
        end
    end
    display_progress_bar(' done',false)
end
goods = neighbors_spatial_correlations>0 & neighbors_centroid_distances>0;
[a,b,c] = histcounts2(neighbors_spatial_correlations(goods), neighbors_centroid_distances(goods), 50);
% registered_cells_spatial_correlations(assigned_count+1:end)=[];
% non_registered_cells_spatial_correlations(non_assigned_count+1:end)=[];
% non_registered_cells_spatial_correlations(non_registered_cells_spatial_correlations<0.01)=[];


end

