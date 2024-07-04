function fH = makeFigure7_PredictedSuppressionLevels(ds,lmmResults, lmmResults_Model,roisToPlot,cmapROIs,saveFigs,saveFigDir)
% Function to reproduce main manuscript figure 7: 
% Suppression level (LMM  slopes) vs model-based (predicted) suppression
% level by each pRF model (LSS,CSS,CST)
%
% From the paper:
% Title:   Rethinking simultaneous suppression in visual cortex via 
%          compressive spatiotemporal population receptive fields.
% Authors: Kupers, Kim, Grill-Spector (2024).
% Journal: Nature Communications
% DOI:     XXX
%
% Requires getting MRI data from OSF (see downloadDataTableFromOSF.m)
%
% Code written by E.R. Kupers (2024) Stanford University
% 
% INPUTS (required):
% - lmmResults       : DATA cell (1x number of ROIs), containing a struct with fields:
%                       fixedIntercepts, fixedSlopes, 
%                       fixedIntercepts_CI, fixedSlopes_CI
% - lmmResults_Model : MODELS cell (3x number of ROIs), containing a struct with fields:
%                       fixedIntercepts, fixedSlopes, 
%                       fixedIntercepts_CI, fixedSlopes_CI
% - roisToPlot       : cell with ROI names
% - cmapROIs         : color map for ROIs
% - saveFigs         : save figures or not?
% - saveFigDir       : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle
%
%% Plot predicted regression slopes
all_lmmResults = cat(1,lmmResults_Model,lmmResults);

LMMOrder = {'LSS','CSS','CST','Data'};
useSTRetParams = false;
fH = plotLMMfittedRegressionSlopes(ds,all_lmmResults,LMMOrder,roisToPlot,cmapROIs, useSTRetParams, saveFigs, saveFigDir);

return

%% To plot line fits to single voxel predicted BOLD amplitudes (alike Figure 4A)
% modelNames = = {'LSS','CSS','CST'};
% for  mm = 1:3
%     fH = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM_Model(mm,:), ...
%         'fixedIntercepts', lmmResults_Model.fixedIntercepts(mm,:), ...
%         'fixedSlopes', lmmResults_Model.fixedSlopes(mm,:), ...
%         'fixedIntercepts_CI', lmmResults_Model.fixedIntercepts_CI(mm,:), ...
%         'fixedSlopes_CI', lmmResults_Model.fixedSlopes_CI(mm,:), ...
%         'plotAllSubjectsTogether', true, ...
%         'plotModelAmpl', true, ...
%         'whichModelToPlot',modelNames{mm});
% end
