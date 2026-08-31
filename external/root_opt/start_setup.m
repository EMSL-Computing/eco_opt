clear all;
close all;
clc;
%%
addpath(genpath(pwd));

% List of zip files
zipFiles = {'casadi-windows-matlabR2016a-v3.5.5.zip', 'yop-master.zip'};

for i = 1:length(zipFiles)
    zipFile = zipFiles{i};
    
    % Remove .zip extension to get folder name
    [~, folderName, ~] = fileparts(zipFile);
    
    % Destination folder = current directory + folderName
    destFolder = fullfile(pwd, folderName);
    
    % Create destination folder if it doesn't exist
    if ~exist(destFolder, 'dir')
        mkdir(destFolder);
    end
    
    % Unzip into the destination folder
    unzip(zipFile, destFolder);
    
    % Add unzipped folder (and its subfolders) to path
    addpath(genpath(destFolder));
end


%% run test
test_yop(true)
%% run script to generate figure used in the manuscript
Figure2_and_3
Figure_4_soilCN_soil_N
Figure_5_soil_plant
Figure_6_soil_plant
%%


for i = 1:length(zipFiles)
    zipFile = zipFiles{i};
    % Remove .zip extension to get folder name
    [~, folderName, ~] = fileparts(zipFile);
    
    % Destination folder = current directory + folderName
    destFolder = fullfile(pwd, folderName);

    % Add unzipped folder (and its subfolders) to path
    rmpath(genpath(destFolder));
end
