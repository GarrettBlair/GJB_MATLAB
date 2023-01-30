X_data = importdata('Mnist_data.txt');
%X_data = (X_data(1:10,:));
D = squeeze(sqrt(sum(bsxfun(@minus,X_data,reshape(X_data',1,size(X_data,2),size(X_data,1))).^2,2)));
%D = L2_distance1(X_data', X_data', 1);
% Both upper and lower D will be performed same


%x = load('swiss_roll_data');
%D = L2_distance(X_data(:,1:1000), X_data(:,1:1000), 1);
options.dims = 1:2;
[Y, R, E] = Isomap(D, 'k', 7, options); 

%%
fname =   'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_water\Hipp18240\processed_files/2022_09_24_H16_46_19_TR20_@placecells.mat';
load(fname, 'ms')
X_data = double(ms.neuron.S_matw);
D = squeeze(sqrt(sum(bsxfun(@minus,X_data,reshape(X_data',1,size(X_data,2),size(X_data,1))).^2,2)));
options.dims = 1:3;
[Y, R, E] = Isomap(D, 'k', 7, options); 
