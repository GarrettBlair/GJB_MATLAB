dd = ["C:\Users\gjb326\Downloads\pkmz_hpcaccratsanalysis\datapoints_new_TF4_02.mat",...
"C:\Users\gjb326\Downloads\pkmz_hpcaccratsanalysis\datapoints_new_TF5_02.mat",...
"C:\Users\gjb326\Downloads\pkmz_hpcaccratsanalysis\datapoints_new_HPCACC24504_01.mat",...
"C:\Users\gjb326\Downloads\pkmz_hpcaccratsanalysis\datapoints_new_HPCACC24505_01.mat",...
"C:\Users\gjb326\Downloads\pkmz_hpcaccratsanalysis\datapoints_new_HPCACC24514_03.mat",...
"C:\Users\gjb326\Downloads\pkmz_hpcaccratsanalysis\datapoints_new_TF1_01.mat",...
"C:\Users\gjb326\Downloads\pkmz_hpcaccratsanalysis\datapoints_new_TF2_03.mat",...
"C:\Users\gjb326\Downloads\pkmz_hpcaccratsanalysis\datapoints_new_TF3_02.mat"];

regs = ["cx", "cc", "ori", "py", "srad", "lac", "ub", "lb", "thl"];
normreg = "cc";
% regs = ["cc", "thl"];
    figure(1); clf; hold on
for i = 1:length(dd)-2
    temp = load(dd(i));
    
    pkmz = NaN(length(regs),1);
    gfp  = NaN(length(regs),1);
    dapi = NaN(length(regs),1);
    eval(sprintf('norm_val = mean(temp.%s, 1);', normreg));
    for j = 1:length(regs)
        eval(sprintf('x = temp.%s;', regs(j)))
        pkmz(j) = mean(x(:,1)); % red channel pkmz
        gfp(j)  = mean(x(:,2)); % green channel gfp
        dapi(j) = mean(x(:,3)); % blue channel dapi
    end
    subplot(3,2,1); hold on
    plot(pkmz, 1:j, 'ro-')
    subplot(3,2,3); hold on
    plot(gfp, 1:j, 'go-')
    subplot(3,2,5); hold on
    plot(dapi, 1:j, 'bo-')
    
    subplot(3,2,2); hold on
    plot(pkmz/norm_val(1), 1:j, 'ro-')
    subplot(3,2,4); hold on
    plot(gfp/norm_val(2), 1:j, 'go-')
    subplot(3,2,6); hold on
    plot(dapi/norm_val(3), 1:j, 'bo-')
end
for i = length(dd)-2:length(dd)
    temp = load(dd(i));
    
    pkmz = NaN(length(regs),1);
    gfp  = NaN(length(regs),1);
    dapi = NaN(length(regs),1);
    eval(sprintf('norm_val = mean(temp.%s, 1);', normreg));
    for j = 1:length(regs)
        eval(sprintf('x = temp.%s;', regs(j)))
        pkmz(j) = mean(x(:,1)); % red channel pkmz
        gfp(j)  = mean(x(:,2)); % green channel gfp
        dapi(j) = mean(x(:,3)); % blue channel dapi
    end
    subplot(3,2,1); hold on
    plot(pkmz, 1:j, 'rx:')
    subplot(3,2,3); hold on
    plot(gfp, 1:j, 'gx:')
    subplot(3,2,5); hold on
    plot(dapi, 1:j, 'bx:')
    
    subplot(3,2,2); hold on
    plot(pkmz/norm_val(1), 1:j, 'rx:')
    subplot(3,2,4); hold on
    plot(gfp/norm_val(2), 1:j, 'gx:')
    subplot(3,2,6); hold on
    plot(dapi/norm_val(3), 1:j, 'bx:')
end
        subplot(3,2,1); hold on
set(gca, 'YTick', 1:length(regs), 'YTickLabel', regs(:), 'YDir', 'reverse')
        subplot(3,2,2); hold on
set(gca, 'YTick', 1:length(regs), 'YTickLabel', regs(:), 'YDir', 'reverse')
        subplot(3,2,3); hold on
set(gca, 'YTick', 1:length(regs), 'YTickLabel', regs(:), 'YDir', 'reverse')
        subplot(3,2,4); hold on
set(gca, 'YTick', 1:length(regs), 'YTickLabel', regs(:), 'YDir', 'reverse')
        subplot(3,2,5); hold on
set(gca, 'YTick', 1:length(regs), 'YTickLabel', regs(:), 'YDir', 'reverse')
        subplot(3,2,6); hold on
set(gca, 'YTick', 1:length(regs), 'YTickLabel', regs(:), 'YDir', 'reverse')

%%