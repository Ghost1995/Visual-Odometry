%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates the fundamental matrix using Normalized 8-Point
% algorithm where RANSAC is used based on Zhang's Robust Estimation Method.
% 
% Input:
%      oldPoints --> Points in frame 1
%      newPoints --> Corresponding points in frame 2
%        imgSize --> Size of each frame
%   numIteration --> Number of iterations for RANSAC
% 
% Output:
%         F --> Fundamental matrix
%   inliers --> Inliers used to generate the fundamental matrix
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, inliers] = estFundamentalMatrix(oldPoints, newPoints, imgSize, numIteration)

    % Convert points to homogeneous coordinates
    oldPoints = [oldPoints, ones(length(oldPoints),1)];
    newPoints = [newPoints, ones(length(newPoints),1)];

    % Divide image into a 8x8 grid
    imgSize(1) = 0.8*imgSize(1); % To reduce outliers on hood of car
    res = imgSize/8;
    segment = cell(8,8);
    segmentLength = zeros(8);
    for i = 1:8
        for j = 1:8
            segment{i,j} = find((oldPoints(:,1)>((j-1)*res(2)))&(oldPoints(:,1)<=(j*res(2)))&...
                                   (oldPoints(:,2)>((i-1)*res(1)))&(oldPoints(:,2)<=(i*res(1))));
            segmentLength(i,j) = numel(segment{i,j});
        end
    end
    filledSegments = find(segmentLength~=0);
    
    % Get which index of which grid to be used
    count = 0;
    gridComb = zeros(numIteration,8);
    indexComb = zeros(numIteration,8);
    while count < numIteration
        count = count+1;
        gridComb(count,:) = filledSegments(randperm(length(filledSegments),8));
        for index = 1:8
            i = rem(gridComb(count,index),8);
            if i == 0
                i = 8;
            end
            j = ceil(gridComb(count,index)/8);
            indexComb(count,index) = randi(segmentLength(i,j));
        end
        if count~=1
            sameGrid = find(all(gridComb(1:count-1,:)==gridComb(count,:),2));
            if ~isempty(sameGrid)
                if any(all(indexComb(sameGrid,:)==indexComb(count,:),2))
                    count = count - 1;
                end
            end
        end
    end
    i = rem(gridComb,8);
    i(i==0) = 8;
    j = ceil(gridComb/8);
    
    % Run RANSAC
    points1 = zeros(8,3);
    points2 = zeros(8,3);
    numInliers = 0;
    inliers = [];
    for count = 1:numIteration
        % Get the points to be used
        for index = 1:8
            points1(index,:) = oldPoints(segment{i(count,index),j(count,index)}(indexComb(count,index)),:);
            points2(index,:) = newPoints(segment{i(count,index),j(count,index)}(indexComb(count,index)),:);
        end
        
        % Estimate the fundamental matrix
        F = norm8PointAlgorithm(points1, points2);
        
        % Calculate error
        epipole1 = F*oldPoints';
        epipole2 = F'*newPoints';
        error = sum(newPoints*F.*oldPoints,2).^2./(sum(epipole1(1:2,:).^2)'+sum(epipole2(1:2,:).^2)');
        currInliers = error <= 0.01;
        if numInliers < sum(currInliers)
            numInliers = sum(currInliers);
            inliers = currInliers;
        end
    end
    
    % Estimate the fundamental matrix using inliers
    F = norm8PointAlgorithm(oldPoints(inliers,:), newPoints(inliers,:));

end