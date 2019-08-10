% ENPM 673 Project 2 -  Visual Odometry
% Author : Pranav Inani
clc;
clear all;
% Extract and intrinsic Camera Parametes
[fx, fy, cx, cy, G_camera_image, LUT] = ReadCameraModel('..\input\stereo\centre','../input/model');
K = [fx 0 cx;0 fy cy; 0 0 1];
R = eye(3);
t = [0;0;0];
cordinates = [0,0,0]';
% Set paths to read images sequentially
image_dir = '..\input\stereo\centre';
filePattern = sprintf('%s/*.png', image_dir);
baseFileNames = dir(filePattern);
file_names = {baseFileNames.name};
numberOfImageFiles = length(baseFileNames);
Img_1_name = fullfile(image_dir,cell2mat(file_names(10)));
% Read 1st image and perform pre-processing
img1 = imread(Img_1_name);
rgb_img1 = demosaic(img1,'gbrg');
undist_img1 = UndistortImage(rgb_img1, LUT);
gray_img1 = rgb2gray(undist_img1);
% Extrct SURF features for image 1
points1 = detectSURFFeatures(gray_img1);
[features1, valid_points1] = extractFeatures(gray_img1, points1,'Upright', true);
for p = 11:(numberOfImageFiles)
    % Read 2nd image and perform pre-processing
    Img_2_name = fullfile(image_dir,cell2mat(file_names(p)));
    img2 = imread(Img_2_name);
    rgb_img2 = demosaic(img2,'gbrg');
    undist_img2 = UndistortImage(rgb_img2, LUT);
    gray_img2 = rgb2gray(undist_img2);
    % Extrct SURF features for image 2
    points2 = detectSURFFeatures(gray_img2);
    [features2, valid_points2] = extractFeatures(gray_img2, points2,'Upright', true);
    % Match Features and Extract points
    indexPairs = matchFeatures(features1,features2, 'Unique', true, 'MaxRatio', 0.3);
    matchedPoints1 = valid_points1(indexPairs(:,1),:);
    matchedPoints2 = valid_points2(indexPairs(:,2),:);
    % Randomly choose 8 point correspondences 
    k = randperm(length(indexPairs));
    sPts1 = matchedPoints1(k(1:8),:).Location;
    sPts2 = matchedPoints2(k(1:8),:).Location;
    
    % Estimate Fundamental and Essential Matrix
    E = EstimateEssentialMatrix(sPts1,sPts2,K);
    % Extract Rotation and Translation 
    [Rnew,tnew] = EstimateRotTrans(sPts1,sPts2,E,K);
    
    % Update rotation and translation
    t = t + R*tnew;
    R = R*Rnew;
    % Store and plot cordinates
    cordinates = [cordinates t];
    plot(cordinates(1,:),cordinates(3,:)), drawnow

    % Make img2 -> img1
    features1 = features2;
    valid_points1 = valid_points2;
end