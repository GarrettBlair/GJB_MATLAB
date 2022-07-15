function [kern] = gb_kernel(nsize, type)
% nsize - half width of kernel
switch type
    case 'ones'
        kern = ones(1, nsize*2 + 1);
    case 'step'
        kern = [zeros(1, nsize) ones(1, nsize+1)];
    case 'step1'
        a = linspace(1, 0, nsize+1);
        kern = [zeros(1, nsize) a];
    case 'step2'
        a = linspace(1, 0, nsize+1);
        w = (a.^2);
        kern = [zeros(1, nsize) w];
    case 'step3'
        a = linspace(1, 0, nsize+1);
        w = (a.^3);
        kern = [zeros(1, nsize) w];        
end
kern = kern./(sum(kern(:)));