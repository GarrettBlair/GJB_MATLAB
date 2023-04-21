% Set parameters (these could be arguments to a function)
rInner = 80;     % inner radius of the colour ring
rOuter = 200;    % outer radius of the colour ring
ncols = 12;      % number of colour segments
% Get polar coordinates of each point in the domain
[x, y] = meshgrid(-rOuter:rOuter);
[theta, rho] = cart2pol(x, y);
% Set up colour wheel in hsv space
hue = (theta + pi) / (2 * pi);     % hue into range (0, 1]
hue = ceil(hue * ncols) / ncols;   % quantise hue 
saturation = ones(size(hue));      % full saturation
brightness = double(rho >= rInner & rho <= rOuter);  % black outside ring
% Convert to rgb space for display
rgb = hsv2rgb(cat(3, hue, saturation, brightness));
% Check result
imshow(rgb);