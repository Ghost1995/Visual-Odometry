%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates absolute transformation matrix using MATLAB
% in-built functions.
% 
% Input:
%      oldPoints --> Points in frame 1
%      newPoints --> Corresponding points in frame 2
%   numIteration --> Maximum number of iterations for RANSAC
%              K --> Intrinsic matrix
%              H --> Absolute transformation matrix for frame i-1
% 
% Output:
%       H --> Absolute transformation matrix for frame i
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H,F] = getPose_MATLAB(oldPoints, newPoints, numIteration, K, H)

    % Estimate fundamental matrix
    [F,inliersIndex] = estimateFundamentalMatrix(oldPoints,newPoints,'Method','RANSAC','NumTrials',numIteration);
    
    % Get the relative rotation and translation
    [relativeRot,relativeTransl] = relativeCameraPose(F,K,oldPoints(inliersIndex),newPoints(inliersIndex));
    
    % Compute absolute rotation and translation
    H = H*[relativeRot' relativeTransl'; 0 0 0 1];

end