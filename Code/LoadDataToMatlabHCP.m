clc
clear all
close all

% Read data from text files and save in Matlab formats

site = 'HCP';
%pipeline = 'ICAFIX';
pipeline = 'STANDARD';

numberOfSubjects = 100;
numberOfTimepoints = 300;
numberOfParcels = 200;

allData = zeros(numberOfSubjects,numberOfParcels,numberOfTimepoints);

currentDirectory = ['rois_cc200/HCP/' pipeline '/'];

% Read data
files = dir([currentDirectory '*.1D']);
delimiterIn = '\t';
headerlinesIn = 1;

disp('Reading data')
for subject = 1:numberOfSubjects
    subject
    filename = files(subject).name
    data = importdata([currentDirectory filename],delimiterIn,headerlinesIn);
    temp = (data.data)';
    
    allData(subject,:,:) = temp(:,1:300);
end

save([site '_' pipeline '_cc200.mat'],'allData')

