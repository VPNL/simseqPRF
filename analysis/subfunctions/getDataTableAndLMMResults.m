function [ds, fLMM, fLMM_Model, lmmResults, lmmResults_Model] = ...
    getDataTableAndLMMResults(projectDir, lmmModelFlag, saveFolderNameDS)
% Function to create and run main linear mixed model (LMM) on data and
% model mean stimulus block amplitudes. For comparing LMMs to more
% conservative models, see s_compareLMM.m. Note: Running model LMMs take
% ~20 min, so ideally should be run once and then stored
% 
% INPUTS:
% projectDir       : (string) base of project
% lmmModelFlag     : (boolean) fit LMM to predicted BOLD amplitudes or not,
%                       for each pRF model (LSS, CSS, CST)
% saveFolderNameDS : (string) foldername where to store data table and LMM
%                       results. If undefined, no data will be saved.
%
% OUTPUTS:
% ds               : (dataset) Table with mean amplitudes, pRF params and
%                       more for each subject's voxel
% fLMM             : (cell) LMM MATLAB object, fitted to Data (1 x 9 ROIs)
% fLMM_Model       : (cell) LMM MATLAB object, fitted to Model predictions (3 pRF models x 9 ROIs)
% lmmResults       : (struct) LMM Observed slopes and intercepts (1 x 9 ROIs), 
%                             each slope/intercept field has 4 conditions
% lmmResults_Model : (struct) LMM Predicted slopes and intercepts (3 models x 9 ROIs), 
%                             each slope/intercept field has 4 conditions
% LMM_params       : (struct) Structure with fitted model string, nr of
%                             fixed and random variables and model names

%% Check inputs
if ~exist('lmmModelFlag','var') || isempty(lmmModelFlag)
    lmmModelFlag = false;
end

if ~exist('saveFolderNameDS','var') || isempty(saveFolderNameDS)
    saveLMMResults = false;
end

% Create data table
postFix = '_variableBlockOnset_gridFit2';
ds = createDataTable(projectDir,postFix);

% Reorder ROIs
roisToPlot = unique(ds.ROI,'stable');
newROIOrder = [1,2,3,4,5,8,9,6,7];
roisToPlot = roisToPlot(newROIOrder);

%% Fit LMM to data
% Random intercepts and slopes per subject and condition
fLMM = cell(1,length(roisToPlot));
LMM_params.LMMlabel = 'feIntcptSlopeAmplCond_reIntcptSlopeConditionInteraction_GroupSubj';
LMM_params.fitStr_Data = 'MeanSimAmp ~ MeanSeqAmp * Condition + (1 + MeanSeqAmp * Condition | Subject)';

LMM_params.nrConditions    = length(unique(ds.Condition));
LMM_params.nrSubjects      = length(unique(ds.Subject));
LMM_params.nrFixedEffectLevels = [nrConditions,nrConditions];
LMM_params.nrRandomEffectLevels = [nrConditions,nrConditions];

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
        getLMMcoefficients(ds(ds.ROI==roisToPlot(ii),:),LMM_params.fitStr_Data, ...
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

if lmmModelFlag % Fit the same model to predicted mean Sim and Seq amplitudes
    
    %% Fit LMM to Model Amplitudes
    fLMM_Model = cell(3,length(roisToPlot));
    
    LMM_params.LMMlabel_Model{1} = 'fe_IntcptSlopeAmplLSSCond_reIntcptSlopeConditionInteraction_GroupSubj';
    LMM_params.fitStr_Model{1} = 'MeanSimAmpModelLSS ~ MeanSeqAmpModelLSS * Condition + (1 + MeanSeqAmpModelLSS * Condition | Subject)';
    
    LMM_params.LMMlabel_Model{2} = 'fe_IntcptSlopeAmplCSSCond_reIntcptSlopeConditionInteraction_GroupSubj';
    LMM_params.fitStr_Model{2} = 'MeanSimAmpModelCSS ~ MeanSeqAmpModelCSS * Condition + (1 + MeanSeqAmpModelCSS * Condition | Subject)';
    
    LMM_params.LMMlabel_Model{3} = 'fe_IntcptSlopeAmplCSTCond_reIntcptSlopeConditionInteraction_GroupSubj';
    LMM_params.fitStr_Model{3} = 'MeanSimAmpModelCST ~ MeanSeqAmpModelCST * Condition + (1 + MeanSeqAmpModelCST * Condition | Subject)';
    
    LMM_params.modelNames = {'LSS','CSS','CST'};
    for ii = 1:length(roisToPlot)
        for mm = 1:length(LMM_params.modelNames)
            [fLMM_Model{mm,ii},...
                fixedIntercepts_Model{mm,ii},fixedIntercepts_CI_Model{mm,ii}, fixedIntercepts_SE_Model{mm,ii},...
                fixedSlopes_Model{mm,ii}, fixedSlopes_CI_Model{mm,ii}, fixedSlopes_SE_Model{mm,ii}, ...
                subjIntercepts_Model{mm,ii}, subjIntercepts_CI_Model{mm,ii}, subjIntercepts_SE_Model{mm,ii},...
                subjSlopes_Model{mm,ii}, subjSlopes_CI_Model{mm,ii}, subjSlopes_SE_Model{mm,ii}] = ...
                getLMMcoefficients(ds(ds.ROI==roisToPlot(ii),:),LMM_params.fitStr_Model{mm}, ...
                LMM_params.nrFixedEffectLevels,LMM_params.nrRandomEffectLevels);
        end
    end
    
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
else
    fLMM_Model = {};
    lmmResults_Model = struct();
end

% Save data table and lmm results if foldername is defined
if saveLMMResults
    if lmmModelFlag
        fNameDS = fullfile(sprintf(saveFolderNameDS,...
            sprintf('SIMSEQ_dataTable_%s.mat',datestr(now,'yyyymmdd'))));
        save(fNameDS,'ds','fLMM','fLMM_Model',...
            'lmmResults','lmmResults_Model','LMM_params')
    else
        fNameDS = fullfile(saveFolderNameDS,...
            sprintf('SIMSEQ_dataOnlyTable_%s.mat',datestr(now,'yyyymmdd')));
        save(fNameDS,'ds','fLMM','LMMlabel','lmmResults',...
            'LMM_params')
    end
end