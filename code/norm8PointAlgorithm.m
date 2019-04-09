%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates the fundamental matrix using Normalized 8-Point
% algorithm.
% 
% Input:
%   oldPoints --> Points in frame 1
%   newPoints --> Corresponding points in frame 2
% 
% Output:
%         F --> Fundamental matrix
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = norm8PointAlgorithm(oldPoints, newPoints)

    % Normalize points
    [oldPoints,oldT] = normalizePoints(oldPoints);
    [newPoints,newT] = normalizePoints(newPoints);

    % Estimate the fundamental matrix for normalized coordinates
    A = [newPoints(:,1)*ones(1,3), newPoints(:,2)*ones(1,3), ones(length(oldPoints),3)];
    A(:,1:3:7) = A(:,1:3:7).*oldPoints(:,1);
    A(:,2:3:8) = A(:,2:3:8).*oldPoints(:,2);
    [~,~,V] = svd(A,0);
    F = reshape(V(:,end),[3 3])';
    
    % Enforce rank 2 constraint on the fundamental matrix
    [U,S,V] = svd(F);
    S = diag([(S(1,1)+S(2,2))/2 (S(1,1)+S(2,2))/2 0]);
    F = U*S*V';

    % De-normalize the fundamental matrix
    F = newT'*F*oldT;
    
    % Normalize the fundamental matrix
    F = F/norm(F);
    if F(end) < 0
        F = -F;
    end

end