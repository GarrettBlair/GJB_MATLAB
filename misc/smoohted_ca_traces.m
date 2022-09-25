load('C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\2022_06_18___16_44_33_ms_placecells_data.mat')
%
cells2use = [1 7 8 53 98 117 140 141 149 155 168 169];
ind = 5000:10000;
yra = ms.neuron.YrA(cells2use,ind);
c = ms.neuron.C(cells2use,ind);
% n = n(cells2use, ind);
t = ms.timestamps(ind)./1000; t = t-t(1);
% stacked_traces(n, 3)
%%
s = 3;
figure(2); clf; hold on
set(gcf, 'Position', [100, 100, 200*3.66, 200*2.16])
yc = yra*0;
n = normalize_rows(yra+c);
for i = 1:size(n,1)
yc(i,:) = movmedian(yra(i,:), [5 6]);%, [2 3]);    
end
nc = normalize_rows(yc+c);
ncc=nc*0;
for i = 1:size(n,1)
ncc(i,:) = movmedian(nc(i,:), [2 3]);%, [2 3]);    
end
for i = 1:size(n,1)
%     plot(t, s*(n(i,:))+i-1, 'r')
    plot(t, s*(ncc(i,:))+i-1, 'k', 'LineWidth', 1)
end
axis tight
set(gca, 'XTick', [30:30:600])