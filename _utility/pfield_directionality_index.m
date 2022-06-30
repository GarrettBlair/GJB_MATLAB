function [D] = pfield_directionality_index(pfields1, pfields2)
%% Compute field directionality. High if a only has a field in one dir, low if field in two dir
% a la Ravassard, Mehta et al. 2013

D = NaN(size(pfields1,1),1);
for i = 1:size(pfields1, 1)
    p1 = pfields1(i,:); % field in one dir
    p2 = pfields2(i,:); % field in opp dir
    
    pdiff = abs(sum(p1-p2));
    psim  = sum(p1+p2);
    D(i)  = pdiff/psim;
%     cla
%     hold on
%     plot(p2)
%     plot(p1)
%     title(D(i))
%     drawnow
%     input('');
end
    