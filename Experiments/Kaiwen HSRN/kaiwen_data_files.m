ddir = 'C:\Users\gjb326\Desktop\Kaiwen data\';
fs = {'2022_09_19_H16_37_23_TR14_@placecells.mat',...
    '2022_09_19_H17_07_47_DRK15_@placecells.mat',...
    '2022_09_21_H17_49_55_TR17_@placecells.mat',...
    '2022_09_21_H18_22_12_DRK18_@placecells.mat',...
    '2023_06_20_H13_11_53_TR15_@placecells_HPC_miniscope1.mat',...
    '2023_06_20_H13_33_34_CON16_@placecells_HPC_miniscope1.mat'}; % intentionally not light vs dark
outfs = {'light_1','dark_1','light_2','dark_2', 'light_3', 'dark_3'};

out_vars = {'arena_x', 'arena_y' 'filename' 'room_x' 'room_y' 'spks' 't'};
in_vars  = {'ms.arena.x', 'ms.arena.y' 'ms.filename' 'ms.room.x' 'ms.room.y' 'ms.spks' 'ms.timestamps'};

for ii = 1:length(fs)
    clearvars ms
    load(sprintf('%s%s', ddir, fs{ii}))
    raw = ms.neuron.YrA+ms.neuron.C;
    room_x = ms.room.x;
    room_y = ms.room.y;
    arena_x = ms.arena.x;
    arena_y = ms.arena.y;
    filename = ms.fileName;
    
    raw_traces = ms.neuron.C + ms.neuron.YrA;
    [dff_filt, ~] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, true);
    spks = normalize_rows(dff_filt);
    spks(isnan(spks)) = 0;
    t = ms.timestamps;
    save(sprintf('%s%s.mat', ddir, outfs{ii}), out_vars{:})
end