clear all
close all
clc

addpath('Packages/FSLNets')
addpath(sprintf('%s/etc/matlab',getenv('FSLDIR')))

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
partialCorrelationRegularization = 1;

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
FWEs1 = zeros(2,1,1,3,2,10);
FWEs2 = zeros(2,1,1,3,2,10);
FWEs_twosided = zeros(2,1,1,3,2,10);

errors1_total = zeros(1,1,1,3,2,10);
errors2_total = zeros(1,1,1,3,2,10);

numberOfNansControls = zeros(1,1,1,3,2,10,numberOfControls,simulations);
numberOfNansDiseased = zeros(1,1,1,3,2,10,numberOfDiseased,simulations);

%-------
site = 'HCP';
numberOfSubjects = 100;
numberOfTimepoints = 300;

site_ = 1;
preprocessing_ = 1;

for pipeline_ = 1:2
    
    if pipeline_ == 1
        pipeline = 'ICAFIX';
    elseif pipeline_ == 2
        pipeline = 'STANDARD';
    end
    
    % Load current data
    load([site '_' pipeline '_cc200.mat'])
    
    % Change number of time points
    allData = allData(:,:,1:numberOfTimepoints);
    
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
                FWE_twosided = zeros(simulations,1);
                
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
                        % More realistic scrubbing
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
                    
                    % Randomly apply scrubbing to each diseased subject,
                    % and calculate correlation matrix
                    for subject = 1:numberOfDiseased
                        
                        % Randomize how many time points to save
                        scrubbingPortion = abs(round((meanScrubbingDiseased + stdScrubbingDiseased * randn)));
                        savedTimepoints = round(numberOfTimepoints * (100 - scrubbingPortion)/100);
                        savedTimepoints = min(savedTimepoints,numberOfTimepoints); % We cannot save more time points than we start with
                        savedTimepoints = max(savedTimepoints,50); %Always save at least 50 time points
                        
                        % Do random scrubbing
                        %keep = [ones(savedTimepoints,1); zeros(numberOfTimepoints - savedTimepoints,1)];
                        %keep = keep(randperm(numberOfTimepoints));
                        % More realistic scrubbing
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
                        
                        FWEthreshold_twosided = icdf('t',1-0.025/numberOfTests,numberOfControls + numberOfDiseased - 2);
                        
                        if max(t_scores(:)) > FWEthreshold
                            FWE1(simulation) = 1;
                        end
                        
                        if min(t_scores(:)) < -1*FWEthreshold
                            FWE2(simulation) = 1;
                        end
                        
                        if max(t_scores(:)) > FWEthreshold_twosided
                            FWE_twosided(simulation) = 1;
                        elseif min(t_scores(:)) < -1*FWEthreshold_twosided
                            FWE_twosided(simulation) = 1;
                        end
                        
                    end
                    
                end
                
                FWEs1(pipeline_,preprocessing_,site_,fisher,correlationType,iteration) = sum(FWE1)/simulations * 100;
                
                FWEs2(pipeline_,preprocessing_,site_,fisher,correlationType,iteration) = sum(FWE2)/simulations * 100;
                
                FWEs_twosided(pipeline_,preprocessing_,site_,fisher,correlationType,iteration) = sum(FWE_twosided)/simulations * 100;
                
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

FWEs1

FWEs2

save(['FWEs_HCP_' num2str(numberOfParcels) 'parcels_partialregularisation_' num2str(partialCorrelationRegularization) '_newscrubbingfunction.mat'],'FWEs1','FWEs2','FWEs_twosided')

delete(myPool)


