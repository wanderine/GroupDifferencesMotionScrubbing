close all
clear all
clc

addition = '_realisticscrubbingfunction_twosided';

load(['FWEs_50parcels_partialregularisation_0.1' addition '.mat']);

meanScrubbingDiseased = [5 10 15 20 25 30 35 40 45 50]; % Percent

twosided = 0;

% pipeline, preprocessing, site, convert to z,
% correlation type, mean scrubbing diseased
%FWEs1 = zeros(4,4,2,2,2,10);

% if preprocessing_ == 1
%     preprocessing = 'filt_global';
% elseif preprocessing_ == 2
%     preprocessing = 'filt_noglobal';
% elseif preprocessing_ == 3
%     preprocessing = 'nofilt_global';
% elseif preprocessing_ == 4
%     preprocessing = 'nofilt_noglobal';
% end

for site_ = 1:2
    
    if site_ == 1
        site = 'New York';
        siteFilename = 'NewYork';
    elseif site_ == 2
        site = 'Michigan';
        siteFilename = 'Michigan';
    end
    
    for pipeline_ = 1:3
        
        if pipeline_ == 1
            pipeline = 'CCS';
        elseif pipeline_ == 2
            pipeline = 'CPAC';
        elseif pipeline_ == 3
            pipeline = 'DPARSF';
        elseif pipeline_ == 4
            pipeline = 'NIAK';
        end
        
        for fisher = 2:3
            
            if fisher == 1
                convert_to_z = 0;
                fisherText = 'no Fisher r-to-z';
                fisherFilename = 'noFisher';
            elseif fisher == 2
                convert_to_z = 1;
                fisherText = 'corrected Fisher r-to-z';
                fisherFilename = 'Fisher';
            elseif fisher == 3
                convert_to_z = -1;
                fisherText = 'raw Fisher r-to-z';
                fisherFilename = 'rawFisher';
            end
            
            for correlationType = 1:2
                
                if correlationType == 1
                    full_correlation = 1;
                    correlationTypeText = 'full correlation';
                    correlationTypeFilename = 'fullcorrelation';
                elseif correlationType == 2
                    full_correlation = 0;
                    correlationTypeText = 'partial correlation';
                    correlationTypeFilename = 'partialcorrelation';
                end
                
                figure
                
                if twosided == 0
                    plot(meanScrubbingDiseased,squeeze(FWEs1(pipeline_,1,site_,fisher,correlationType,:)),'r','LineWidth',3)
                    hold on
                    plot(meanScrubbingDiseased,squeeze(FWEs1(pipeline_,2,site_,fisher,correlationType,:)),'g','LineWidth',3)
                    hold on
                    plot(meanScrubbingDiseased,squeeze(FWEs1(pipeline_,3,site_,fisher,correlationType,:)),'b','LineWidth',3)
                    hold on
                    plot(meanScrubbingDiseased,squeeze(FWEs1(pipeline_,4,site_,fisher,correlationType,:)),'y','LineWidth',3)
                    hold on
                    plot(meanScrubbingDiseased,5*ones(size(meanScrubbingDiseased)),'k','LineWidth',3)
                    hold off
                elseif twosided == 1
                    plot(meanScrubbingDiseased,squeeze(FWEs_twosided(pipeline_,1,site_,fisher,correlationType,:)),'r','LineWidth',3)
                    hold on
                    plot(meanScrubbingDiseased,squeeze(FWEs_twosided(pipeline_,2,site_,fisher,correlationType,:)),'g','LineWidth',3)
                    hold on
                    plot(meanScrubbingDiseased,squeeze(FWEs_twosided(pipeline_,3,site_,fisher,correlationType,:)),'b','LineWidth',3)
                    hold on
                    plot(meanScrubbingDiseased,squeeze(FWEs_twosided(pipeline_,4,site_,fisher,correlationType,:)),'y','LineWidth',3)
                    hold on
                    plot(meanScrubbingDiseased,5*ones(size(meanScrubbingDiseased)),'k','LineWidth',3)
                    hold off
                end
                
                title(sprintf('%s, %s, 20 subjects per group, \n%s, %s, 50 parcels\nMean scrubbing for first group is 10 percent',site,pipeline,correlationTypeText,fisherText))
                
                xlabel('Mean scrubbing for second group (%)')
                ylabel('Familywise error rate (%)')
                if correlationType == 1 && fisher == 2
                    legend('Bandpass filtering, global signal regression','Bandpass filtering, no global signal regression', 'No bandpass filtering, global signal regression', 'No bandpass filtering, no global signal regression','Expected','Location','NorthWest')
                end
                set(gca,'FontSize',15)
                axis([5 50 0 100])
                
                print([siteFilename '_' pipeline '_' correlationTypeFilename '_' fisherFilename  '_FWEs.png'],'-dpng')
                
            end
        end
    end
end




% pipeline, preprocessing, site, convert to z,
% correlation type, mean scrubbing diseased
%FWEs1 = zeros(4,4,2,2,2,10);

% Site

temp = FWEs1(1:3,:,1,2:3,:,:);
meanNewYork = mean(temp(:))

temp = FWEs1(1:3,:,2,2:3,:,:);
meanMichigan = mean(temp(:))

% Bandpass filtering or not

temp = FWEs1(1:3,1:2,:,2:3,:,:);
meanFiltering = mean(temp(:))

temp = FWEs1(1:3,3:4,:,2:3,:,:);
meanNoFiltering = mean(temp(:))

% GSR or not

temp = FWEs1(1:3,[1 3],:,2:3,:,:);
meanGSR = mean(temp(:))

temp = FWEs1(1:3,[2 4],:,2:3,:,:);
meanNoGSR = mean(temp(:))

% Full vs partial correlation

temp = FWEs1(1:3,:,:,2:3,1,:);
meanFull = mean(temp(:))

temp = FWEs1(1:3,:,:,2:3,2,:);
meanPartial = mean(temp(:))


% Pipelines

temp = FWEs1(1,:,:,2:3,:,:);
meanCCS = mean(temp(:))

temp = FWEs1(2,:,:,2:3,:,:);
meanCPAC = mean(temp(:))

temp = FWEs1(3,:,:,2:3,:,:);
meanDPARSF = mean(temp(:))






