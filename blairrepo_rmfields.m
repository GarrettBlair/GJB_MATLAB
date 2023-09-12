badfields={'numgoodcells','infoLR','infoRL','isplaceLR_N','isplaceRL_N','isplace_N','isboth_N','isonlyLR_N','isonlyRL_N','isplaceRL_N','isrecurringplace','isplace'};    
cd 'D:\Sample Data\Blair et al 2023\tadblair\Blair_et_al_DATA\pretrn'
temp = dir('*_predata.mat');


for i = 1:length(temp)
    fname = temp(i).name;
    clearvars predata
    load(fname);
    if isfield(predata, badfields(1))
        disp(fname)
        predata = rmfield(predata,badfields);
        save(fname, 'predata');
    else
        disp('none')
    end
end


%%
clear
badfields={'C'};    
cd 'D:\Sample Data\Blair et al 2023\tadblair\Blair_et_al_DATA\sessiondata'
temp = dir('*_sess.mat');

for i = 1:length(temp)
    fname = temp(i).name;
    clearvars frame*
    load(fname);
    a = whos('frame*');
    ss = sprintf('x = %s;', a.name(:));
    eval(ss(:))
    if isfield(x, badfields(1))
        disp(fname)
        x = rmfield(x, badfields);
        ss = sprintf('%s = x;', a.name(:));
        eval(ss(:))
        save(fname, 'frame*');
    else
        disp('none')
    end
end
