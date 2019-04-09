%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function normalizes the points using a unique transformation. The
% normalized points have the center of mass at [0,0] and their mean
% distance from the center is sqrt(2).
% 
% Input:
%   points --> Points to be normalized
% 
% Output:
%   normPoints --> Normalized points
%            T --> Transformation matrix used to normalize the points
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [normPoints, T] = normalizePoints(points)

    % Compute centroid
    centroid = mean(points);
    
    % Translate points so that the centroid is at [0,0]
    translatedPoints = points - centroid;
    
    % Compute the scale to make mean distance from centroid equal to sqrt(2)
    meanDistance = mean(sqrt(sum(translatedPoints.^2,2)));
    % Protect against division by 0
    if meanDistance > 0
        scale = sqrt(2)/meanDistance;
    else
        scale = 1;
    end
    
    % Compute the matrix to scale
    T = [scale, 0, -scale*centroid(1); 0, scale, -scale*centroid(2); 0, 0, 1];
    
    % Compute the normalized homogeneous points
    normPoints = (T*points')';

end