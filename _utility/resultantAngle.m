function resultantAngle = resultantAngle(angles)
% resultanAngleFromDistribution calculates the resultant angle from a
% distribution of angles.
%
%   resultantAngle = resultantAngleFromDistribution(angles)
%
%   INPUT:
%       angles: A vector of angles in degrees or radians.
%
%   OUTPUT:
%       resultantAngle: The resultant angle in degrees.

    % Convert angles to radians
    anglesRad = deg2rad(angles);

    % Calculate the x and y components of each angle
    xComponents = cos(anglesRad);
    yComponents = sin(anglesRad);

    % Sum the x and y components
    sumX = sum(xComponents, 'omitnan');
    sumY = sum(yComponents, 'omitnan');

    % Calculate the resultant angle in radians
    resultantAngleRad = atan2(sumY, sumX);

    % Convert the resultant angle to degrees
    resultantAngle = rad2deg(resultantAngleRad);
end