%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function gives the relative rotation and translation from the
% fundamental matrix computed in the last step
% 
% Input:
%           F --> Fundamental matrix
%           K --> Intrinsic matrix
%   oldPoints --> Points in frame 1
%   newPoints --> Corresponding points in frame 2
% 
% Output:
%      Rotation --> Relative rotation
%   Translation --> Relative translation
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rotation, Translation] = relativeCamPose(F, K, oldPoints, newPoints)

    % Compute essential matrix
    E = K'*F*K;
    
    % Enforce rank 2 constraint
    [U,S,V] = svd(E);
    S = diag([(S(1,1)+S(2,2))/2 (S(1,1)+S(2,2))/2 0]);
    E = U*S*V';

    [U,~,V] = svd(E);
    % Get all possible rotations
    W = [0 -1 0; 1 0 0; 0 0 1];
    R = cat(3,U*W'*V',U*W*V');
    % Force rotations to have determinant of +1
    for i = 1:2
        if det(R(:,:,i))<0
            R(:,:,i) = -R(:,:,i);
        end
    end
    % Get all possible translations
    T = repmat([U(:,3), -U(:,3)],[1,2]);
    T(:,T(3,:)>0) = [];
    
    % Get calibrated coordinates for old frame
    oldPoints = K\[oldPoints'; ones(1,length(oldPoints))];
    oldPoints = (oldPoints./oldPoints(3,:))';

    % Get calibrated coordinates for new frame
    newPoints = K\[newPoints'; ones(1,length(newPoints))];
    newPoints = (newPoints./newPoints(3,:))';

    % Get the correct rotation and translation
    for i=1:length(oldPoints)
        index = false(1,2);
        for j = 1:2
            % Compute camera matrix
            oldP = [eye(3) zeros(3,1)];
            newP = inv([R(:,:,j) T(:,j); 0 0 0 1]);

            % Compute point in 3-D
            A = [oldPoints(i,1)*oldP(3,:) - oldP(1,:); oldPoints(i,2)*oldP(3,:) - oldP(2,:);...
                 newPoints(i,1)*newP(3,:) - newP(1,:); newPoints(i,2)*newP(3,:) - newP(2,:)];
            [~,~,V] = svd(A,0);
            X = V(:,end);
            newX = newP\X;
            if X(4)~=0
                X = X(1:3)/X(4);
            end
            if newX(4)~=0
                newX = newX(1:3)/newX(4);
            end
            % Check if the point is in front of both cameras
            if (X(3)>0)&&(newX(3)>0)
                index(j) = true;
            end
        end
        if sum(index)==1
            newP = inv([R(:,:,index) T(:,index); 0 0 0 1]);
            Rotation = newP(1:3,1:3);
            Translation = newP(1:3,4);
            return;
        end
    end
    Rotation = eye(3);
    Translation = zeros(3,1);

end