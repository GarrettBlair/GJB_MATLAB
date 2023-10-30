clear
%%
topdir = 'D:\Sample Data\Ketamine\conjoint_ipos\';
%
ni = 14;
nj = 4;
figure(78); clf
m=NaN(ni,nj);
        figure(179); clf
for i = 1:ni
        figure(78+i); clf
        figure(179); 
        subplot_tight(ni, 1, i); hold on
    for j = 1:nj
        fn = sprintf('%sfr_i%d_j%d.mat', topdir, i, j);
        fn2 = sprintf('%sstandard\\standard_i%d_j%d.mat', topdir, i, j);
        disp(fn)
        clearvars temp* e
        temp = load(fn);
        temp2 = load(fn2);
        nsm(i,j) = size(temp.ensembl,2);
        temp.ensembl(temp.ensembl>.1) = .1;
        e = temp.ensembl./max(temp.ensembl);
        
%         subplot_tight(nj,ni, ni*(j-1)+i)
        figure(78+i); 
        subplot_tight(nj, 1, j)
        hold on
        plot(e)
%         subplot_tight(nj, 2, j*2)
        imagesc(temp.bin_spk)
        plot(e*size(temp.bin_spk,1))
        axis tight
        figure(179); 
%         histogram(e, [0.001:.0001:.02], 'Normalization', 'probability', 'FaceAlpha', .2)
        hc = histcounts(e, [0.000:.0001:.02], 'Normalization', 'probability');
        plot([0.000:.0001:.02-.00005], cumsum(hc))
        xlim([-.001 .021])
    end
end

