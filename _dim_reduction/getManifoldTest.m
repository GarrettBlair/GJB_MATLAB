%get manifold embedding
% Computes Isomap embedding using the algorithm of Tenenbaum, de Silva, and Langford (2000). 
clear; close all; clc

% addpath('IsoMap/');

apFS = 30000; %sample rate

timestampsNP = 1:100*apFS; %timestamps for spikes (10 seconds)
        
%create gaussian with 100ms SD
tg = -0.5*apFS:1:0.5*apFS; %500ms
sd = 2; %sdev
mn = 0; %mean
yg = gaussmf(tg,[0.1*apFS mn]); %*100ms
%plot(tg,yg)

%create matrix with rows = cells, columns = time
Ncells = 20;
Nsamples = length(timestampsNP);

spks = zeros(Ncells,Nsamples);

%populate each cell with 100 random spikes
Nspikes = 1000;
for i = 1:Ncells
    r = randsample(Nsamples,Nspikes);
    spks(i,r) = 1;
end

%create firing rate profile - convolve with gaussian

frx = zeros(size(spks));

for i = 1:Ncells

    x = spks(i,:);
    xc = conv(x,yg); %convolution with gaussian

    %fix offset
    N = round((length(yg)-1)/2);
    xc = xc(N:length(xc)-N); 
    if length(xc) > length(timestampsNP)
        xc = xc(1:length(timestampsNP));
    end

    %sqrt to stabilize variance
    xc = xc.^0.5;

    frx(i,:) = xc;
end
        
    
%computing full embedding taking advantage of sparsity, allowing larger datasets

%downsample 1000x from 30kHz to 30Hz
downSampleRatio = 1000;

%takes more memory - use for server
%frxD = downsample(frx',downSampleRatio)';
        
%takes less memory - use for laptop
N = ceil(size(frx,2)/downSampleRatio);
frxD = zeros(size(frx,1),N);
for i = 1:size(frx,1)
   x = frx(i,:);
   x = downsample(x,downSampleRatio);
   frxD(i,:) = x;
end
    
%downsample timestamps
txD = timestampsNP(1:downSampleRatio:end);
    

global X %globally defined data

X = frxD;

%remove neurons with no spikes
%s = sum(X,2);
%X = X(s>0,:);

options.dims = 3; %1:3 %range of dimensions to calculate typ 2:20 for dim calculations, 3 for plot
options.landmarks = 1:size(X,2); %recommended value
nNeighbours = 5;
[Y, R, E] = IsomapII('dfun', 'k', nNeighbours, options); 

%3D plot
figure; hold on;
dimID = 1; %index of dim to plot from options.dims
plot3(Y.coords{dimID}(1,:),Y.coords{dimID}(2,:),Y.coords{dimID}(3,:))

