%%
ddir = 'D:\Sample Data\ensemble_ipos\prop_mod\';
props = [0 1 5 10 20 40:10:90 95 99 100];
props_e_mod = NaN(length(props),1);
props_i_mod = NaN(length(props),1);
for i = 1:length(props)
    fn = sprintf('%sipos_prop%d.mat', ddir, props(i));
    temp = load(fn);
    e1 = temp.ensem1;
    e2 = temp.ensem2;
    e = abs(e1-e2);
    i1 = nanmean(temp.ipos1, 1);
    i2 = nanmean(temp.ipos2, 1);
    ipos = abs(abs(i1) - abs(i2));
    s = temp.switch_binned;
    props_e_mod(i) = nansum(e(s==1))/nansum(e);
    props_i_mod(i) = nansum(ipos(s==1))/nansum(ipos);

%     figure(i); clf;
%     subplot(211)
%     plot(e1 + 1);
%     hold on
%     plot(e2 - 1)
%     plot(e1-e2)
%     
%     subplot(212)
%     plot(i1 + 2);
%     hold on
%     plot(i2 - 2)
%     plot(abs(i1)-abs(i2))
end
%%
ddir = 'D:\Sample Data\ensemble_ipos\segs_mod\';
nsegs = [1 5 10 25 50 100 250 500];
nsegs_e_mod = NaN(length(nsegs),1);
nsegs_i_mod = NaN(length(nsegs),1);
rep = 0;
for i = 1:length(nsegs)
    fn = sprintf('%sipos_nsegs%d_%d.mat', ddir, nsegs(i), rep);
    temp = load(fn);
    e1 = temp.ensem1;
    e2 = temp.ensem2;
    e = abs(e1-e2);
    i1 = nanmean(temp.ipos1, 1);
    i2 = nanmean(temp.ipos2, 1);
    ipos = abs(abs(i1) - abs(i2));
    s = temp.switch_binned;
    nsegs_e_mod(i) = nansum(e(s==1))/nansum(e);
    nsegs_i_mod(i) = nansum(ipos(s==1))/nansum(ipos);
    
%     figure(i); clf;
%     subplot(211)
%     plot(e1 + 1);
%     hold on
%     plot(e2 - 1)
%     plot(e)
%     
%     subplot(212)
%     plot(i1 + 2);
%     hold on
%     plot(i2 - 2)
%     plot(ipos)
end
%%
figure(101)
subplot(131); cla; hold on
plot(props, props_i_mod, 'ko-')
plot(props, props_e_mod, 'bo-')
axis([-10 110 -.1 1.1])
subplot(132); cla; hold on
plot(log2(nsegs), nsegs_i_mod, 'ko-')
plot(log2(nsegs), nsegs_e_mod, 'bo-')
axis([-1 log2(700) -.1 1.1])

