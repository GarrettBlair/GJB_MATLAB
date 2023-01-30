%%% Analysis Garrett

load('G:\My Drive\Colabs\Garrett\GJB_Hipp8_ex.mat')
opts.threshold.method = 'circularshift';
opts.threshold.number_of_permutations = 200;
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;
opts.threshold.permutations_percentile = 95;

Patterns = assembly_patterns(spk,opts);

Activity = assembly_activity(Patterns,spk);

%% Plotting assembly activity

time = t/1000;

Xlims = [3.5e2 7.5e2];
Xlims2 = [4.55e2 4.65e2];
cells1 = [265 272 283 312 328 329 357];
cells2 = [200 211 213 218 219 226 229 259];
Ylims = [180 380];

h(1)=sub(5,3,5,1:2);
g=plot(time,Activity([15 24],:)','linewidth',2);
hold on
plot([Xlims2;Xlims2],[0 0;600 600],'r--','linewidth',2)
hold off
ylabel({'Assembly','Activation Strength'})
xlabel('Time (s)')
legend(g,{'Assembly 15','Assembly 24'},'box','off','location','northwest')
set(gca,'fontsize',12)
h(2)=sub(5,3,1:4,1:2);
imagesc(time,[],spk)
hold on
plot([Xlims2;Xlims2],[0 0;900 900],'r--','linewidth',2)
hold off
ylabel('Cell #')
title('Raster Plot and Assembly Activity')
set(gca,'xticklabel','','fontsize',12)
linkaxes(h,'x')
xlim(Xlims)

h(1)=sub(5,3,5,3);
plot(time,Activity([15 24],:)','linewidth',2)
set(gca,'fontsize',12)
ylim([0 600])
h(2)=sub(5,3,1:4,3);
imagesc(time,[],spk)
hold on
k=plot(Xlims2',[cells1;cells1],'linewidth',0.2,'color',[0 0 1]);
kk=plot(Xlims2',[cells2;cells2],'linewidth',0.2,'color',[1 0 0]);
hold off
linkaxes(h,'x')
set(gca,'xticklabel','','fontsize',12)
ylim(Ylims)
xlim(Xlims2)
legend([k(1) kk(1)],{'Cells from Assembly 15','Cells from Assembly 24'},'box','off',...
    'location','northeastoutside')


%% Calculating Assemblies in subset of data
part = floor(length(spk)/3);
spk1 = spk(:,1:part);
spk2 = spk(:,2*part+1:end);

Patterns1 = assembly_patterns(spk1,opts);
Patterns2 = assembly_patterns(spk2,opts);

[Corr12]=AssemblySI_fast(Patterns1,Patterns2);


%% Plotting patterns
Caxis = [-.42 .42];

sub(1,4,1,1);
imagesc(Patterns)
title({'Assembly Patterns','Calculated using all Data'})
ylabel('Cell #')
caxis(Caxis)
set(gca,'fontsize',12)
hold on
plot([14 14 16 16],[0 738 738 0],'color',[.5 .5 .5])
plot([23 23 25 25],[0 738 738 0],'color',[.5 .5 .5])
hold off
sub(1,4,1,2);
imagesc(Patterns1)
title({'Assembly Patterns','Calculated using First 1/3'})
xlabel('Assembly #')
caxis(Caxis)
set(gca,'fontsize',12)

sub(1,4,1,3);
imagesc(Patterns2)
colorbar
title({'Assembly Patterns','Calculated using Last 1/3'})
text(34,size(spk,1)/2,'Assembly Weights','rotation',270,'horizontalalignment',...
    'center','fontsize',12)
caxis(Caxis)
set(gca,'fontsize',12)

sub(3,5,2,5)
imagesc(Corr12)
caxis([0 1])
ylabel('Assemblies of First 1/3')
xlabel('Assemblies of Last 1/3')
text(37,15,'Similarity Index','rotation',270,'horizontalalignment','center','fontsize',12)
title({'Are assemblies from the First 1/3','the same as the Last 1/30?'})
colorbar
set(gca,'fontsize',12)



