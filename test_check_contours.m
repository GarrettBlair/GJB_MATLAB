%%
animals = {'Hipp18240'};%{'Hipp16942'};
numAnimals = length(animals);
analysis_version = 'v1.0';
dir_list_fname = '_directory_list.csv';
experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_water/';



vars2plot = {'ms.neuron.minFrame', 'ms.neuron.maxFrame', 'ms.neuron.pnr_im', 'contour_im'};
for animalLoop = 1:numAnimals
    %%
    animal_name = animals{animalLoop};
    dir_file            = sprintf('%s%s/%s%s', experiment_folder, animal_name, animal_name, dir_list_fname);
%     processedDir        = sprintf('%s%s/processed_files/', experiment_folder, animal_name);
%     contourDir        = sprintf('%s%s/matching_contours/', experiment_folder, animal_name);
    
    AnimalDir = setup_imaging_Sessionfiles(animal_name, dir_file, experiment_folder);%processedDir, contourDir);
    numsess = size(AnimalDir.processedFile,1);
    nr = length(vars2plot);
    nc = numsess;
    figure(animalLoop); clf
    for i = 1:nc
        for j = 1:nr
            fname = AnimalDir.processedFile{i};
            fname_c = AnimalDir.contourFile{i};
            temp = load(fname);
            temp_c = load(fname_c);
            
            if strcmp(vars2plot{j}, 'contour_im')
                x = squeeze(sum(temp_c.contours, 1));
            else
                eval(sprintf('x = temp.%s;', vars2plot{j}))
            end
            ind = i + (j-1)*nc;
            subplot_tight(nr, nc, ind, [.01 .01])
            imagesc(x)
            axis image off
            drawnow
        end
    end
    
end