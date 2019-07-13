clc
clear all
close all

% Read data from text files and save in Matlab formats

for pipeline_ = 1:1
    
    if pipeline_ == 1
        pipeline = 'ccs';
    elseif pipeline_ == 2
        pipeline = 'cpac';
    elseif pipeline_ == 3
        pipeline = 'dparsf';
    elseif pipeline_ == 4
        pipeline = 'niak';
    end
    
    for preprocessing_ = 1:4
        
        if preprocessing_ == 1
            preprocessing = 'filt_global';
        elseif preprocessing_ == 2
            preprocessing = 'filt_noglobal';
        elseif preprocessing_ == 3
            preprocessing = 'nofilt_global';
        elseif preprocessing_ == 4
            preprocessing = 'nofilt_noglobal';
        end
        
        for site_ = 2:2
            
            if site_ == 1
                site = 'NYU';
                numberOfSubjects = 171;
                if pipeline_ == 1
                    numberOfTimepoints = 175;
                elseif pipeline_ == 2
                    numberOfTimepoints = 176;
                elseif pipeline_ == 3
                    numberOfTimepoints = 175;
                elseif pipeline_ == 4
                    numberOfTimepoints = 176;
                end
            elseif site_ == 2
                site = 'UM_1';
                numberOfSubjects = 82;
                if pipeline_ == 1
                    numberOfTimepoints = 295;
                elseif pipeline_ == 3
                    numberOfTimepoints = 295;
                else
                    numberOfTimepoints = 296;
                end
            end
            
            numberOfParcels = 200;
            
            allData = zeros(numberOfSubjects,numberOfParcels,numberOfTimepoints);
            
            currentDirectory = ['rois_cc200/Outputs/' pipeline '/' preprocessing '/rois_cc200/'];
            
            % Read data
            files = dir([currentDirectory site '*']);
            delimiterIn = '\t';
            headerlinesIn = 1;
            
            disp('Reading data')
            for subject = 1:numberOfSubjects
                subject
                filename = files(subject).name;
                data = importdata([currentDirectory filename],delimiterIn,headerlinesIn);
                temp = (data.data)';  
                if pipeline_ == 4
                    [np,nt] = size(temp); % For niak
                    allData(subject,:,1:nt) = temp; % 1:nt for niak, time points vary by subject
                else
                    allData(subject,:,:) = temp;
                end
            end
            
            save([site '_' pipeline '_' preprocessing '_cc200.mat'],'allData')
            
        end
    end
end









