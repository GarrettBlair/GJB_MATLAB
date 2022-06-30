function [ xfix, yfix ] = smooth_jumps( xpos, ypos, kernel, plotting )
% WORK IN PROGRESS
% Takes in position vectors X and Y and smooths them with a kernal of size
% 2*K on either side of a large positional jump. Jumps are defined by a
% pixel difference greater than the mean plus 2.5 standard deviations.
%
% Outputs are vectors the same legnth as x and y inputs, with smoothed
% jumps.

num_stds = 2.5;
%% Check inputs
switch nargin == 4
    case 0
        fprintf(['ERROR : Number of input arguments must be 5:\n'...
            'Linear_maze_controller(animal_name, test_time, session_num, context, plotting)\n\n']);
        beep;
        return
    case 1
        if  ~isa(xpos, 'numeric')
            fprintf('ERROR : Improper X input (1st input; must be a number).\n')
            return
        elseif ~isa(ypos, 'numeric')
            fprintf('ERROR : Improper Y input (2st input; must be a number).\n')
            return
        elseif ~isa(kernel, 'numeric')
            fprintf('ERROR : Improper K (smoothing Kernel) (3rd input; must be numeric).\n')
            return
        elseif ~isa(plotting, 'logical')
            fprintf('ERROR : Improper ''plotting'' (4th input; must be logical).\n')
            return
        end
end

%%   Detailed explanation goes here
dist = sqrt(diff(xpos).^2 + diff(ypos).^2);
dist_threshold = mean(dist) + num_stds*std(dist);
jumps = dist >= dist_threshold;
jump_ind = find(jumps);
fprintf('Jumps detected = %d\n', sum(jumps))
xfix = xpos;
yfix = ypos;

for i = 1:2:length(jump_ind)-1
    try
    xfix(jump_ind(i)+1:jump_ind(i+1)) = mean([xpos(jump_ind(i)-kernel):xpos(jump_ind(i))...
        xpos(jump_ind(i+1)+1):xpos(jump_ind(i+1+kernel))]);
    yfix(jump_ind(i)+1:jump_ind(i+1)) = mean([yfix(jump_ind(i)-kernel):yfix(jump_ind(i))...
        yfix(jump_ind(i+1)+1):yfix(jump_ind(i+1+kernel))]);
    catch
        %
    end
end

if plotting == true
    xy = figure; clf
    hold on;
    plot(xpos, ypos, 'k.');
    plot(xfix, yfix, 'g.');
    hold off
    
    % figure(xy); hold on; plot(time(jump_ind), xfix(jump_ind), 'g.')
    %
    % figure(8); clf; plot(time, xpos, 'k.')
    % figure(8); hold on; plot(time(jump_ind), xpos(jump_ind), 'r.')
    % figure(8); hold on; plot(time(jump_ind), xfix(jump_ind), 'g.')
end
end

