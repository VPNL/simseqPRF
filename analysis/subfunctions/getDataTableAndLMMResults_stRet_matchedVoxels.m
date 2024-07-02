function [allDS, allLMMResults, allLMMResults_model] = getDataTableAndLMMResults_stRet_matchedVoxels(projectDir,...
    subjnrs, sModels, tModels, stRetParamsFlag, saveFolderNameDS)

%% Define model params
if ~exist('sModels','var') || isempty(sModels)
    sModels = {'oneGaussianFit','oneGaussianFit','oneGaussianFit'};
end
if ~exist('tModels','var') || isempty(tModels)
    tModels = {'3ch-stLN','3ch-stLN','1ch-dcts'};
end
if ~exist('stRetParamsFlag','var') || isempty(stRetParamsFlag)
    stRetParamsFlag = true;
end

cmapROIs = getROISummaryColors(0);

saveFigs = true;

% Create data table with matched voxels between CST_fix, CST_opt & DN-ST ROIs
allDS = createDataTableMatchedSTRet(projectDir);

% Reorder ROI names
if ~exist('roisToPlot','var')
    allRoisToPlot = unique(allDS{1}.ds.ROI,'stable');
    roiNames = ["V1","V2","V3","hV4","VO1/VO2","V3AB","IPS0/IPS1","LO1/LO2","TO1/TO2"];
    for ii = 1:length(roiNames)
        newROIOrder(ii) = find(ismember(string(allRoisToPlot),roiNames{ii}));
    end
    roisToPlot = allRoisToPlot(newROIOrder);
end

for ii = 1:length(roisToPlot)
    for mm = 1:3
        tmp = allDS{mm}.ds(allDS{mm}.ds.ROI==roisToPlot(ii),:);
        subjMatch{mm,ii} = unique(tmp.Subject);
    end
end
%% Fit LMM to data
allLMMResults = cell(1,3);

for mm = 1:3
    % Random intercepts and slopes per subject and condition
    fLMM = cell(1,length(roisToPlot));
    LMM_params.LMMlabel = 'feIntcptSlopeAmplCond_reIntcptSlopeConditionInteraction_GroupSubj';
    LMM_params.fitStr_Data = 'MeanSimAmp ~ MeanSeqAmp * Condition + (1 + MeanSeqAmp * Condition | Subject)';
    
    LMM_params.nrConditions    = length(unique(allDS{mm}.ds.Condition));
    LMM_params.nrSubjects      = length(unique(allDS{mm}.ds.Subject));
    LMM_params.nrFixedEffectLevels = [LMM_params.nrConditions,LMM_params.nrConditions];
    LMM_params.nrRandomEffectLevels = [LMM_params.nrConditions,LMM_params.nrConditions];
    
    fixedIntercepts   = fLMM;
    fixedIntercepts_CI = fLMM;
    fixedIntercepts_SE = fLMM;
    fixedSlopes       = fLMM;
    fixedSlopes_CI    = fLMM;
    fixedSlopes_SE    = fLMM;
    subjIntercepts    = fLMM;
    subjIntercepts_CI = fLMM;
    subjIntercepts_SE = fLMM;
    subjSlopes        = fLMM;
    subjSlopes_CI     = fLMM;
    subjSlopes_SE     = fLMM;
    
    for ii = 1:length(roisToPlot)
        [fLMM{ii},fixedIntercepts{ii},fixedIntercepts_CI{ii}, fixedIntercepts_SE{ii},...
            fixedSlopes{ii}, fixedSlopes_CI{ii}, fixedSlopes_SE{ii}, ...
            subjIntercepts{ii}, subjIntercepts_CI{ii}, subjIntercepts_SE{ii},...
            subjSlopes{ii}, subjSlopes_CI{ii}, subjSlopes_SE{ii}] = ...
            getLMMcoefficients(allDS{mm}.ds(allDS{mm}.ds.ROI==roisToPlot(ii),:),LMM_params.fitStr_Data, ...
            LMM_params.nrFixedEffectLevels,LMM_params.nrRandomEffectLevels);
    end
    
    lmmResults = struct();
    lmmResults.fixedIntercepts   = fixedIntercepts;
    lmmResults.fixedIntercepts_CI = fixedIntercepts_CI;
    lmmResults.fixedIntercepts_SE = fixedIntercepts_SE;
    lmmResults.fixedSlopes       = fixedSlopes;
    lmmResults.fixedSlopes_CI    = fixedSlopes_CI;
    lmmResults.fixedSlopes_SE    = fixedSlopes_SE;
    lmmResults.subjIntercepts    = subjIntercepts;
    lmmResults.subjIntercepts_CI = subjIntercepts_CI;
    lmmResults.subjIntercepts_SE = subjIntercepts_SE;
    lmmResults.subjSlopes        = subjSlopes;
    lmmResults.subjSlopes_CI     = subjSlopes_CI;
    lmmResults.subjSlopes_SE     = subjSlopes_SE;
    
    allLMMResults{mm} = lmmResults;
end


%% Fit LMM to Model Amplitudes
modelNames = {'CST','CST_ST', 'DN_ST'};
LMMlabel_Model{1} = 'fe_IntcptSlopeAmplCST_Cond_reIntcptSlopeConditionInteraction_GroupSubj';
LMMlabel_Model{2} = 'fe_IntcptSlopeAmplCST_STCond_reIntcptSlopeConditionInteraction_GroupSubj';
LMMlabel_Model{3} = 'fe_IntcptSlopeAmplDN_STCond_reIntcptSlopeConditionInteraction_GroupSubj';
fitStr_Model{1} = 'MeanSimAmpModelCST ~ MeanSeqAmpModelCST * Condition + (1 + MeanSeqAmpModelCST * Condition | Subject)';
fitStr_Model{2} = 'MeanSimAmpModelCST_ST ~ MeanSeqAmpModelCST_ST * Condition + (1 + MeanSeqAmpModelCST_ST * Condition | Subject)';
fitStr_Model{3} = 'MeanSimAmpModelDN_ST ~ MeanSeqAmpModelDN_ST * Condition + (1 + MeanSeqAmpModelDN_ST * Condition | Subject)';

allLMMResults_model = {1,3};

for mm = 1:3
    lmmResults_Model = struct();
    fLMM_Model = cell(1,length(roisToPlot));
    
    LMM_params.modelNames = modelNames(mm);
    LMM_params.LMMlabel_Model = LMMlabel_Model{mm};
    LMM_params.fitStr_Model = fitStr_Model{mm};
    
    for ii = 1:length(roisToPlot)
        
        [fLMM_Model{ii},...
            fixedIntercepts_Model{ii},fixedIntercepts_CI_Model{ii}, fixedIntercepts_SE_Model{ii},...
            fixedSlopes_Model{ii}, fixedSlopes_CI_Model{ii}, fixedSlopes_SE_Model{ii}, ...
            subjIntercepts_Model{ii}, subjIntercepts_CI_Model{ii}, subjIntercepts_SE_Model{ii},...
            subjSlopes_Model{ii}, subjSlopes_CI_Model{ii}, subjSlopes_SE_Model{ii}] = ...
            getLMMcoefficients(allDS{mm}.ds(allDS{mm}.ds.ROI==roisToPlot(ii),:),LMM_params.fitStr_Model, ...
            LMM_params.nrFixedEffectLevels,LMM_params.nrRandomEffectLevels);
        
        lmmResults_Model.fixedIntercepts   = fixedIntercepts_Model;
        lmmResults_Model.fixedIntercepts_CI = fixedIntercepts_CI_Model;
        lmmResults_Model.fixedIntercepts_SE = fixedIntercepts_SE_Model;
        lmmResults_Model.fixedSlopes       = fixedSlopes_Model;
        lmmResults_Model.fixedSlopes_CI    = fixedSlopes_CI_Model;
        lmmResults_Model.fixedSlopes_SE    = fixedSlopes_SE_Model;
        lmmResults_Model.subjIntercepts    = subjIntercepts_Model;
        lmmResults_Model.subjIntercepts_CI = subjIntercepts_CI_Model;
        lmmResults_Model.subjIntercepts_SE = subjIntercepts_SE_Model;
        lmmResults_Model.subjSlopes        = subjSlopes_Model;
        lmmResults_Model.subjSlopes_CI     = subjSlopes_CI_Model;
        lmmResults_Model.subjSlopes_SE     = subjSlopes_SE_Model;
        
        
        
    end
    
    allLMMResults_model{mm} = lmmResults_Model;
end

% Save data table and lmm results if foldername is defined
fNameDS = fullfile(saveFolderNameDS,...
    sprintf('SIMSEQ_dataTable_stRetParams_matchVoxels%s.mat',datestr(now,'yyyymmdd')));

save(fNameDS,'allDS',...
    'allLMMResults','allLMMResults_model','LMM_params')

end
