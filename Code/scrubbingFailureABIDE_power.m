clear all
close all
clc

addpath('Packages/FSLNets')
addpath(sprintf('%s/etc/matlab',getenv('FSLDIR')))

correlationIncrease = 0.09;

% Number of random group analyses per parameter combination
simulations = 1000;

% Assume equal variance or not for parametric statistics
unequalVariance = 0;

% Number of parcels to analyze, do not use all parcels
% to reduce processing time and 
% to make it easier to do partial correlation 
% (due to low number of timepoints)
numberOfParcels = 50;

% Use permutation test or not
permutation = 0;
nPerms = 1000;

% Regularization for ridge p (Tikhonov)
partialCorrelationRegularization = 0.1;

% Number of subjects in each random group
numberOfControls = 20;
numberOfDiseased = 20;

meanScrubbingControls = 10; % Percent

stdScrubbingControls = 5; % Percent
stdScrubbingDiseased = 15; % Percent

numberOfTests = 0;
for p1 = 1:numberOfParcels
    for p2 = 1:numberOfParcels
        if p2 < p1
            numberOfTests = numberOfTests + 1;
        end
    end
end

myPool = parpool(25);

% pipeline, preprocessing, site, convert to z,
% correlation type, mean scrubbing diseased
FWEs1 = zeros(4,4,2,3,2,10);
FWEs2 = zeros(4,4,2,3,2,10);

errors1_total = zeros(4,4,2,3,2,10);
errors2_total = zeros(4,4,2,3,2,10);

numberOfNansControls = zeros(4,4,2,3,2,10);
numberOfNansDiseased = zeros(4,4,2,3,2,10);

%-------

for pipeline_ = 1:3
    
    if pipeline_ == 1
        pipeline = 'ccs';
    elseif pipeline_ == 2
        pipeline = 'cpac';
    elseif pipeline_ == 3
        pipeline = 'dparsf';
    elseif pipeline_ == 4
        pipeline = 'niak';
    end
    
    for preprocessing_ = [1 2 3 4]
        
        if preprocessing_ == 1
            preprocessing = 'filt_global';
        elseif preprocessing_ == 2
            preprocessing = 'filt_noglobal';
        elseif preprocessing_ == 3
            preprocessing = 'nofilt_global';
        elseif preprocessing_ == 4
            preprocessing = 'nofilt_noglobal';
        end
        
        for site_ = 1:2
            
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
            
            % Load current data
            load([site '_' pipeline '_' preprocessing '_cc200.mat'])
            
            % Find subjects with data for all parcels of interest
            goodSubjects = 0;
            for subject = 1:numberOfSubjects
                tempData = allData(subject,1:numberOfParcels,:);
                if (sum(tempData(:) == 0) == 0)
                    goodSubjects = goodSubjects + 1;
                end
            end
            
            tempAllData = zeros(goodSubjects,200,numberOfTimepoints);
            
            goodSubjects = 0;
            for subject = 1:numberOfSubjects
                tempData = allData(subject,1:numberOfParcels,:);
                if (sum(tempData(:) == 0) == 0)
                    goodSubjects = goodSubjects + 1;
                    tempAllData(goodSubjects,:,:) = allData(subject,:,:);
                end
            end
            
            allData = tempAllData;
            clear tempAllData;
            numberOfSubjects = goodSubjects;
            
            for fisher = 2:3
                
                if fisher == 1
                    convert_to_z = 0; % No conversion to z
                elseif fisher == 2
                    convert_to_z = 1; % Conversion with autocorrelation correction
                elseif fisher == 3 
                    convert_to_z = -1; % Conversion without autocorrelation correction
                end
                
                for correlationType = 1:2
                    
                    if correlationType == 1
                        full_correlation = 1;
                    else
                        full_correlation = 0;
                    end
                    
                    % Loop over different means for scrubbing of diseased
                    iteration = 1;
                    
                    for meanScrubbingDiseased = [5 10 15 20 25 30 35 40 45 50]
                        
                        FWE1 = zeros(simulations,1);
                        FWE2 = zeros(simulations,1);
                        
                        errors1 = zeros(simulations,1);
                        errors2 = zeros(simulations,1);
                        
                        nans1 = zeros(simulations,1);
                        nans2 = zeros(simulations,1);
                        
                        parfor simulation = 1:simulations
                            
                            simulation
                            
                            % Randomly shuffle subjects
                            randomizedSubjects = randperm(numberOfSubjects);
                            
                            % Reset correlation matrices
                            correlationMatrixControls = zeros(numberOfControls,numberOfParcels,numberOfParcels);
                            correlationMatrixDiseased = zeros(numberOfDiseased,numberOfParcels,numberOfParcels);
                            
                            % Create a control group and a diseased group from the
                            % randomly shuffled data, pick the numberOfParcels
                            % first parcels from the data
                            controlData = squeeze(allData(randomizedSubjects(1:numberOfControls),1:numberOfParcels,:));
                            diseasedData = squeeze(allData(randomizedSubjects((numberOfControls+1):(numberOfControls+numberOfDiseased)),1:numberOfParcels,:));
                            
                            % Randomly apply scrubbing to each control subject, and
                            % calculate correlation matrix
                            for subject = 1:numberOfControls
                            
								% Randomize how many time points to save
                                scrubbingPortion = abs(round((meanScrubbingControls + stdScrubbingControls * randn)));
                                savedTimepoints = round(numberOfTimepoints * (100 - scrubbingPortion)/100);
                                savedTimepoints = min(savedTimepoints,numberOfTimepoints); % We cannot save more time points than we start with
                                savedTimepoints = max(savedTimepoints,50); %Always save at least 50 time points

                                % Do random scrubbing
                                %keep = [ones(savedTimepoints,1); zeros(numberOfTimepoints - savedTimepoints,1)];
                                %keep = keep(randperm(numberOfTimepoints));
                                keep = generateRandomScrubbing(numberOfTimepoints,savedTimepoints);
                                scrubbedData = squeeze(controlData(subject,:,keep == 1));
                                
                                if full_correlation == 1
                                    allCorrelations = nets_netmats(scrubbedData',convert_to_z,'corr');
                                else
                                    allCorrelations = nets_netmats(scrubbedData',convert_to_z,'ridgep',partialCorrelationRegularization);
                                end
                                
                                nans1(simulation) = nans1(simulation) + sum(isnan(allCorrelations(:)));
                                
                                correlationMatrixControls(subject,:,:) = allCorrelations;
                            end
                            
                            
                            % Randomize where to increase correlation
                            firstROI = 1;
                            
                            % Randomly apply scrubbing to each diseased subject,
                            % and calculate correlation matrix
                            for subject = 1:numberOfDiseased
                                
                                % Make correlations stronger
                                diseasedData(subject,firstROI:(firstROI+9),:) = diseasedData(subject,firstROI:(firstROI+9),:) + correlationIncrease * diseasedData(subject,(firstROI+10):(firstROI+10+9),:);
                                                              
								% Randomize how many time points to save
                                scrubbingPortion = abs(round((meanScrubbingDiseased + stdScrubbingDiseased * randn)));
                                savedTimepoints = round(numberOfTimepoints * (100 - scrubbingPortion)/100);
                                savedTimepoints = min(savedTimepoints,numberOfTimepoints); % We cannot save more time points than we start with
                                savedTimepoints = max(savedTimepoints,50); %Always save at least 50 time points

                                % Do random scrubbing
                                %keep = [ones(savedTimepoints,1); zeros(numberOfTimepoints - savedTimepoints,1)];
                                %keep = keep(randperm(numberOfTimepoints));
                                keep = generateRandomScrubbing(numberOfTimepoints,savedTimepoints);
                                scrubbedData = squeeze(diseasedData(subject,:,keep == 1));
                                
                                if full_correlation == 1
                                    allCorrelations = nets_netmats(scrubbedData',convert_to_z,'corr');
                                else
                                    allCorrelations = nets_netmats(scrubbedData',convert_to_z,'ridgep',partialCorrelationRegularization);
                                end

                                nans2(simulation) = nans2(simulation) + sum(isnan(allCorrelations(:)));
                                
                                correlationMatrixDiseased(subject,:,:) = allCorrelations;
                            end
                            
                            %% Perform two sample t-test for each parcel
                            
                            if permutation == 1
                                
                                netmats = zeros(numberOfControls + numberOfDiseased, numberOfParcels * numberOfParcels);
                                
                                for subject = 1:numberOfControls
                                    temp = correlationMatrixControls(subject,:,:);
                                    netmats(subject,:) = temp(:);
                                end
                                
                                for subject = 1:numberOfDiseased
                                    temp = correlationMatrixDiseased(subject,:,:);
                                    netmats(subject + numberOfControls,:) = temp(:);
                                end
                                
                                try
                                    [p_uncorrected,p_corrected1] = nets_glm(netmats,'design.mat','design1.con',0,nPerms);
                                    
                                    if max(p_corrected1(:)) > 0.95
                                        FWE1(simulation) = 1;
                                    end
                                    
                                    size(p_corrected1)
                                    
                                catch
                                    errors1(simulation) = 1;
                                end
                                
                                %try
                                %    [p_uncorrected,p_corrected2] = nets_glm(netmats,'design.mat','design2.con',0,nPerms);
                                    
                                %    if max(p_corrected2(:)) > 0.95
                                %        FWE2(simulation) = 1;
                                %    end
                                    
                                %catch
                                %    errors2(simulation) = 1;
                                %end
                                
                            else
                                
                                t_scores = zeros(numberOfParcels);
                                
                                for p1 = 1:numberOfParcels
                                    for p2 = 1:numberOfParcels
                                        if p2 < p1
                                            
                                            controls = correlationMatrixControls(:,p1,p2);
                                            diseased = correlationMatrixDiseased(:,p1,p2);
                                            
                                            try
                                                if unequalVariance == 0
                                                    [H,P,CI,STATS] = ttest2(controls,diseased);
                                                elseif unequalVariance == 1
                                                    [H,P,CI,STATS] = ttest2(controls,diseased,'vartype','unequal');
                                                end
                                                t_scores(p1,p2) = STATS.tstat;
                                            catch
                                                errors1(simulation) = 1;
                                            end
                                            
                                        end
                                    end
                                end
                                
                                threshold = icdf('t',1-0.05,numberOfControls + numberOfDiseased - 2);
                                FWEthreshold = icdf('t',1-0.05/numberOfTests,numberOfControls + numberOfDiseased - 2);
                                
                                if max(t_scores(:)) > FWEthreshold
                                    FWE1(simulation) = 1;
                                end
                                
                                if min(t_scores(:)) < -1*FWEthreshold
                                    FWE2(simulation) = 1;
                                end
                                
                            end
                            
                        end
                        
                        FWEs1(pipeline_,preprocessing_,site_,fisher,correlationType,iteration) = sum(FWE1)/simulations * 100;
                        
                        FWEs2(pipeline_,preprocessing_,site_,fisher,correlationType,iteration) = sum(FWE2)/simulations * 100;
                        
                        errors1_total(pipeline_,preprocessing_,site_,fisher,correlationType,iteration) = sum(errors1);
                        
                        errors2_total(pipeline_,preprocessing_,site_,fisher,correlationType,iteration) = sum(errors2);
                        
                        numberOfNansControls(pipeline_,preprocessing_,site_,fisher,correlationType,iteration) = sum(nans1);
                        
                        numberOfNansDiseased(pipeline_,preprocessing_,site_,fisher,correlationType,iteration) = sum(nans2);
                        
                        
                        % Clean temp directory
                        if permutation == 1
                            system('rm /tmp/*.nii.gz')
                        end
                        
                        iteration = iteration + 1;
                        
                    end
                    
                end
                
            end
        end
    end
end

FWEs1

FWEs2

save(['power_correlationIncrease_' num2str(correlationIncrease)  '_newscrubbingfunction_10modifiedparcels_constant.mat'],'FWEs1','FWEs2','numberOfNansControls','numberOfNansDiseased');


%save(['FWEs_' num2str(numberOfParcels) 'parcels.mat'],'FWEs1','FWEs2')


delete(myPool)




