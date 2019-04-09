%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates the solution for Project 2
% 
% Output:
%   Translation --> Consists of the locations of the camera center from 
%                   both the user-defined and the MATLAB code
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify folders of interest
imgFolder = '..\input\stereo\centre\';
modelFolder = '..\input\model\';

% Get image file names
imgFiles = dir([imgFolder '*.png']);

% Get intrinsic parameters
[fx, fy, cx, cy, ~, LUT] = ReadCameraModel(imgFolder,modelFolder);
K = cameraParameters('IntrinsicMatrix',[fx 0 0; 0 fy 0; cx cy 1]);

% Read the first frame
for i = 1:length(imgFiles)
    Iold = histeq(rgb2gray(demosaic(imread([imgFolder imgFiles(i).name]),'gbrg')));
    if sum(sum(Iold > 200))<0.25*numel(Iold)
        firstIndex = i;
        break;
    end
end

% Read the first frame
Iold_color = UndistortImage(demosaic(imread([imgFolder imgFiles(firstIndex).name]),'gbrg'),LUT);
Iold = histeq(rgb2gray(Iold_color));
points1 = detectSURFFeatures(Iold,'ROI',[1 1 size(Iold,2) 0.8*size(Iold,1)]);
[features1,points1] = extractFeatures(Iold,points1);

% Initiate the figure
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1) % Image frame
imshow(Iold_color);
hold on
plot(points1.Location(:,1),points1.Location(:,2),'gx');
hold off
title('Camera Frames with Detected Features')

subplot(2,1,2) % 2-D plot from vehicle view
plot(0,0,'b');
hold on
plot(0,0,'r');
hold off
xlabel('Motion in x-direction')
ylabel('Motion in z-direction')
title('2-D Motion of the Camera')
legend('Result from User-defined Functions','Result from MATLAB in-built Functions')
legend('boxoff')

% Create a video
vidObj = VideoWriter('../output/trajectory_plot.avi');
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

% Initiate transformation matrices
H = repmat(eye(4),[1,1,length(imgFiles)]);
trueH = H; % for comparison
C = zeros(3,length(imgFiles));
trueC = C; % for comparison
% Define maximum number of iterations for RANSAC
numIteration = 500;

% Compute rotation and translation
for i = firstIndex+1:length(imgFiles)
    % Read the new frame
    Inew_color = UndistortImage(demosaic(imread([imgFolder imgFiles(i).name]),'gbrg'),LUT);
    Inew = histeq(rgb2gray(Inew_color));
    points2 = detectSURFFeatures(Inew,'ROI',[1 1 size(Inew,2) 0.8*size(Inew,1)]);
    [features2,points2] = extractFeatures(Inew,points2);
    
    % Match features
    indexPairs = matchFeatures(features1,features2);
    oldPoints = points1(indexPairs(:,1));
    newPoints = points2(indexPairs(:,2));
    
    % Compute the camera center pose
    H(:,:,i) = getPose_User(oldPoints.Location,newPoints.Location,size(Iold),numIteration,K.IntrinsicMatrix,H(:,:,i-1));
    C(:,i) = H(1:3,4,i);
    trueH(:,:,i) = getPose_MATLAB(oldPoints,newPoints,numIteration,K,trueH(:,:,i-1)); % for comparison
    trueC(:,i) = trueH(1:3,4,i);
    
    % Update the figure
    subplot(2,1,1) % Image frame
    imshow(Inew_color);
    hold on
    plot(newPoints.Location(:,1),newPoints.Location(:,2),'gx');
    hold off
    title('Camera Frames with Detected Features')

    subplot(2,1,2) % 2-D plot from vehicle view
    hold on
    plot(C(1,i-1:i),C(3,i-1:i),'b');
    plot(trueC(1,i-1:i),trueC(3,i-1:i),'r');
    hold off
    pbaspect([2 1 1])
    daspect([1 1 1])
    legend('Result from User-defined Functions','Result from MATLAB in-built Functions')
    legend('boxoff')
    
    % Add updated figure to the video
    writeVideo(vidObj, getframe(gcf));
    
    % Update the reference frame
    features1 = features2;
    points1 = points2;
end
close(vidObj);

% Save the camera centers
T = C(:,firstIndex:end);
T(:,:,2) = trueC(:,firstIndex:end);
save('../output/cameraCenters.mat','T')

% Show the final trajectory generated
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1) % 2-D plot from vehicle view for user-defined functions
plot(T(1,:,1),T(3,:,1),'b')
xlabel('Motion in x-direction')
ylabel('Motion in z-direction')
title('Motion of the Camera Generated using User-defined Functions')


subplot(1,2,2) % 2-D plot from vehicle view for MATLAB in-built functions
plot(T(1,:,2),T(3,:,2),'r')
xlabel('Motion in x-direction')
ylabel('Motion in z-direction')
title('Motion of the Camera Generated using MATLAB in-built Functions')

% Save the motion generated
saveas(gcf,'../output/Motion_of_the_Camera.jpg')