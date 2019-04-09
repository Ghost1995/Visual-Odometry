%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates absolute transformation matrix using user-defined
% functions.
% 
% Input:
%      oldPoints --> Matching points in the frame 1
%      newPoints --> Matching points in the frame 2
%        imgSize --> Size of each frame
%   numIteration --> Number of iterations for RANSAC
%              K --> Lower triangular intrinsic matrix
%              H --> Absolute transformation matrix for frame i-1
% 
% Output:
%       H --> Absolute transformation matrix for frame i
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = getPose_User(oldPoints, newPoints, imgSize, numIteration, K, H)

    % Estimate fundamental matrix
    [F,inliersIndex] = estFundamentalMatrix(oldPoints,newPoints,imgSize,numIteration);
    
    % Get the relative rotation and translation
    [relativeRot,relativeTransl] = relativeCamPose(F,K',oldPoints(inliersIndex,:),newPoints(inliersIndex,:));
    
    % Compute absolute rotation and translation
    H = H*[relativeRot relativeTransl; 0 0 0 1];

end