%% statistical significance in mvpa
% non-parametric technique by combining permutations and bootstrapping

% step1:
% For each subject, the labels of the different conditions (eg. motion_vertical and motion_horizontal) were permuted,
% and the same decoding analysis was performed.
% The previous step was repeated 100 times for each subject.

% DONE in our decoding scripts

% step2:
% A bootstrap procedure was applied in order to obtain a group-level null distribution
% that was representative of the whole group.
% From each subject’s null distribution, one value was randomly chosen (with replacement)
% and averaged across all participants.
% This step was repeated 100,000 times resulting in a group-level null distribution of 100,000 values.

% step3:
% The statistical significance of the MVPA results was estimated by comparing the observed result
% to the group-level null distribution. This was done by calculating the proportion of observations
% in the null distribution that had a classification accuracy higher than the one obtained in the real test.

% step4:
% To account for the multiple comparisons, all p values were corrected using false discovery rate (FDR) correction

clc;
clear;

%% STEP 1: SET THINGS UP

% load the .mat file with decoding accuracies
tic

% Add FDR function from MATLAB FORUM, no mafdr, toolbox licence missing :(
addpath('/Volumes/Slim_Reaper/Projects/Language_MVPA/code/lib/fdr_bh');

decodTitle = 'ReadSpeech_JuBrain';


% This has to be separated from the modality in unimodal! 
%decodingCondition = 'WordPseudoword'; %'WordPseudoword', 'WordControl', 'PseudowordControl';
%decodingModality = 'reading'; % 'reading' 'speech'

decodingCondition = {'WordPseudoword', 'WordControl', 'PseudowordControl'};
decodingModality = {'reading', 'speech'};

roiList = {'Broca', 'FG2', 'F4', 'MTG', 'V1'};
%roiList = {'V1'};

%List the sub numbers
subNum = {'01', '02','03','04','05',...
    '06','07','08','09','10','11','12',...
    '13','14','15','16','17','18','19','20'}; 

%List the groups
subGroup = {'blind', 'sighted'};


%paste numbers with group to create full subjects list
subList = {};
for i = 1:length(subGroup)
    subList = [subList, strcat(subGroup{i}, subNum)];
end

nbRowsInAccu = 30; % HARDCODED: NB ROI * NB MODALITIES * CONDITIONS

im = 'beta'; %'tmap', 'beta'

smooth = 2;

%featureRatio = 200;


% number of iterations for group level null distribution
nbIter = 1000;

% Setup the structure for p values! 
% assign subObsPval to each specific field
for g = 1:numel(subGroup)
    for m = 1:numel(decodingModality)
        for r = 1:numel(roiList)
            for c = 1:numel(decodingCondition)
                
                Pvalues.(subGroup{g}).(decodingModality{m}).(roiList{r}).(decodingCondition{c}) = [];
                
            end
        end
    end
end



% GROUP LOOP SHOULD START HERE? 

% Load acuu from one group only!

accuFile = dir(['/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/cosmo-mvpa/task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas/JuBrain/unimodal/permutations/','sub-',subGroup{1},'*.mat']);

count = 1;

    
accuGroup = struct( ...
        'sub', [], ...
        'roiArea', [], ...
        'roiDimension', [], ...
        'roiNbVoxels', [], ...
        'mvpaFeatures',[], ...
        'ffxResults', [], ...
        'conditions', [], ...
        'modality', [], ...
        'accuracy', [], ...
        'permutation', [], ...
        'pValue', []);
 
%Go through all the loaded file and join their content to a group structure  
% Expected size = 
    
    for iFile = 1:length(accuFile)
        
        load(fullfile(accuFile(iFile).folder, accuFile(iFile).name));
        
        for iRow = 1:nbRowsInAccu
            % store results
            accuGroup(count).sub = accu(iRow).sub;
            accuGroup(count).roiArea = accu(iRow).roiArea;
            accuGroup(count).roiDimension = accu(iRow).roiDimension;
            accuGroup(count).roiNbVoxels = accu(iRow).roiNbVoxels;
            accuGroup(count).mvpaFeatures = accu(iRow).mvpaFeatures;
            accuGroup(count).ffxResults = accu(iRow).ffxResults;
            accuGroup(count).conditions = accu(iRow).conditions;
            accuGroup(count).modality = accu(iRow).modality;
            accuGroup(count).accuracy = accu(iRow).accuracy;
            accuGroup(count).permutation = accu(iRow).permutation;
            accuGroup(count).pValue = accu(iRow).pValue;

                count = count + 1;
                
        end
    end

accu = accuGroup;





for iRoi = 1:length(roiList)
    
    roiLabel = roiList(iRoi);
    
    fprintf('roi: %s \n\n', roiLabel{1})
    
    
    
    %% STEP 2: create group null distribution
    timeStart = datestr(now,'HH:MM');
    
    subSamp = zeros(length(subList), nbIter);
    
% Now go into the group accu and randomly draw nbIter values from
% permutations for each subject in that given modality/decoding condition
    
    for iIter = 1:nbIter
        
        for iAccu = 1:length(accu)
            
            for iSub = 1:length(subList)
                
                subID = subList(iSub);
                                
                if strcmp(char({accu(iAccu).sub}.'), char(subID)) == 1
                    
                    if strcmp(char({accu(iAccu).ffxResults}.'), im) == 1 
                    
                    %check if all the parameters and conditions match
                        
                        if strcmp(string({accu(iAccu).modality}.'), decodingModality)==1 %CHANGED THIS TO MODALITY
                            
                            if strcmp(string({accu(iAccu).conditions}.'), decodingCondition)==1
                            
                                if strcmp(string({accu(iAccu).roiArea}.'),roiLabel) == 1
                                
                                %read the subject level permutations = accu(iAccu).permutation;
                                %pick one decoding accuracy randomly with replacement
                                
                                 if iIter == 1
                                    fprintf('sapmling from: %s sub-%s %s %s \n\n', ...
                                        accu(iAccu).roiArea, ...
                                        accu(iAccu).sub, ...
                                        accu(iAccu).ffxResults, ...
                                        accu(iAccu).modality)
                                 end
                                
                                
                                 subSamp(iSub, iIter) = datasample(accu(iAccu).permutation, 1);
                                 
                                end
                                
                             end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
    %subSamp should now have a size of N subjects x N Iterations 
    %(e.g. 20 x 10000)
    
    timeEnd = datestr(now,'HH:MM');
    
    % This averages all subjects PER iteration, resulting in a 1 x NbIter
    % distribution
    groupNullDistr = mean(subSamp);
    
    disp('step 2 done')
    
    %% STEP 3: check where does the avg accu of the group falls in the group level null ditribution
    % calculate the proportion of values in the group null ditribution which are above the actual decoding
    % accuracy for a one-tailed test. accordingly change for two-tailed test.
    % p = sum(accuracy < acc0) / nbIter; %from Ceren
    % pValue = (sum(Pooled_Perm>accuracy)+1)/(1+NrPermutations); % from Mohamed
    
    subAccu=zeros(length(subList),1);
    
    %subObsPVal=zeros(length(subList),1);
    
    for iAccu=1:length(accu)
        
        for iSub=1:length(subList)
            
            subID=subList(iSub);
            
            if strcmp(char({accu(iAccu).sub}.'),char(subID)) == 1
                
                %check if all the parameters and conditions match
                
                if strcmp(string({accu(iAccu).modality}.'),decodingModality)==1
                    
                    if strcmp(string({accu(iAccu).conditions}.'),decodingCondition)==1
                        
                        if strcmp(string({accu(iAccu).roiArea}.'),roiLabel)==1
                            
                            %read the actual decoding accuracy
                            subAccu(iSub)=[accu(iAccu).accuracy].';
                            
                        end
                    end
                end
            end
        end
    end
    
    % Check how many times group average Accuracy (OBSERVED) is lower than
    % the random one? 
    
    subObsPVal(iRoi) = sum(mean(subAccu)<groupNullDistr)/nbIter;
    
    disp('step 3 done')

    
end

% So P values for a single condition, group and modality make sense! 
% Now what you need to do is to figure out how to make it work for: 
% - multiple groups
% - multiple decoding conditions
% - multiple sensory modalities

% in a way, this is a 4 dimension structure (roi x condition x modality x
% group) 

% So you can either include SAVING in a loop to create separate stats for
% group & modality with 2D p values (ROI x Condition) 
% OR
% Figure out how to pack them into struct arrays, e.g. 
%Pvalues.group.modality.roi.condition

%% STEP 4: correct the obtained p-value

% function mafdr([vector of pvalues], BHFDR, 'true') % from Stefania
%fdrCorrPVal = mafdr(subObsPVal, 'BHFDR', 'true');

% try using fdr_bh function from matlab file exchange 
% https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh
% added to code/lib! 
% Usage:
%  >> [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);


%% save the outout
% save the
% decoding type
% decoding condition
% roi List
% group Null distribution
% subObsPVal
% FDR corrected values

% set output folder/name
pathOutput='/Volumes/MICK/analisys_high-re_multiSensSpatFreq/outputs/derivatives/cosmo-mvpa';

savefileMat = fullfile(pathOutput, ...
    ['stats', '_',  decodingCondition , '_featureRatio-', num2str(featureRatio), '_smooth-', num2str(smooth), '_ffx-', im, '_', datestr(now, 'yyyymmddHHMM'), '.mat']);

% set structure array for keeping the results
% mvpaStats = struct( ...
%             'decodTitle', [], ...
%             'decodCondition', [], ...
%             'roiList', [], ...
%             'groupNullDis', [], ...
%             'obsPVal', [], ...
%             'fdrCorPVal', []);

%store output
mvpaStats.decodTitle = decodTitle;
mvpaStats.decodCondition = decodingCondition;
mvpaStats.roiList = roiList; % this tells the order of corresponding p-values
mvpaStats.groupNullDistr = groupNullDistr;
mvpaStats.obsPVal = subObsPVal;
mvpaStats.fdrCorPVal = fdrCorrPVal;

% mat file
save(savefileMat, 'mvpaStats');

toc