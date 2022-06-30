%% sample script for manifold testing
dfile = 'C:\Users\gblai\MATLAB\18AQ_PR\4_11_2019\H9_M25_S55\CNMFE_NE.mat';

load(dfile)

%%
s = normalize_matrix(neuron.S);
[ncells, nsamples] = size(s);
g = gausswin(15,1); g = g./sum(g);

conv_s = s*0;
for i = 1:ncells
    conv_s(i,:) = conv(s(i,:), g, 'same');
    
end

global X
X = conv_s;
%remove neurons with no spikes
%s = sum(X,2);
%X = X(s>0,:);

options.dims = 1:25; %1:3 %range of dimensions to calculate typ 2:20 for dim calculations, 3 for plot
options.landmarks = 1:size(X,2); %recommended value
options.display = 1;
nNeighbours = 5;
[Y, R, E] = IsomapII('dfun', 'k', nNeighbours, options);
