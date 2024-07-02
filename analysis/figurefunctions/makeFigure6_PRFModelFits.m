function fH = makeFigure6_PRFModelFits(ds,roisToPlot,spatialModel,temporalModel,saveFigs, saveFigDir)
% Function to reproduce main manuscript figure 6: 
% panel a/b: V1 & VO1 time series + LSS, CSS, CST model prediction for each stimulus condition 
% panel c: violin plots showing cross-validated R^2 for LSS, CSS, CST models
% panel d: bar graph showing difference in cross-validated R^2 for LSS, CSS, CST models
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
% - ds              : dataset
% - roisToPlot      : cell with ROI names
% - spatialModel    : spatial components of pRF models
% - temporalModel   : temporal components of pRF models
% - saveFigs        : save figures or not?
% - saveFigDir      : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle

%% Panel A: Plot V1 voxel time series + predicted pRF response per condition
projectDir  = fullfile(simseqRootPath);
subjnr      = 3;
fH(1)       = plotDataModelFit_singleVoxel_StimulusConditionBlocks(projectDir,subjnr,...
                'roisToPlot',{'V1'},'selectedDataVoxels',941, 'plotModelFlag',true);

%% Panel B: Plot VO1/2 voxel time series + predicted pRF response per condition
fH(2)       = plotDataModelFit_singleVoxel_StimulusConditionBlocks(projectDir,subjnr,...
    'roisToPlot',{'VO1'},'selectedDataVoxels',220, 'plotModelFlag',true);

%% Get resampled data
subjnrs = [1:3,7:13];
useSTRetParams = false;
output = resamplePRFParams_wReplacement(ds, roisToPlot, useSTRetParams, temporalModel, spatialModel, subjnrs);
        
%% Panel C: Plot distribution pRF model cross-validated R2 (violin graph)   
fH(3) = plotPRFModelCVRSQ_violinPlot(ds,output.resampledCVR2,output.maxNC,roisToPlot,saveFigs, saveFigDir);

%% Panel D: Plot pRF model difference cross-validated R2 (bar graph)
fH(4) = plotDiffPRFModelCVR2_bargraph(output.mean_diffCVR2,roisToPlot,saveFigs, saveFigDir);

%% Do some statistics (ANOVA)

useAdjustedR2 = false; % Do we want to adjust for the nr of regressors in pRF models 
                       % (CST having 1 more than LSS and CSS).
nrTimePoints = 648;

% Preallocate arrays
allROINames = {};
allModelNames = [];
allSubjectsR2 = [];

for ii = 1:length(roisToPlot)
    groupLSS = ds.R2LSS(ds.ROI==roisToPlot(ii) & ds.Condition==nominal(1));
    groupCSS = ds.R2CSS(ds.ROI==roisToPlot(ii) & ds.Condition==nominal(1));
    groupCST = ds.R2CST(ds.ROI==roisToPlot(ii) & ds.Condition==nominal(1));

    if useAdjustedR2
        groupLSS = 100.*adjustedR2forNumParams(groupLSS./100, nrTimePoints,1);
        groupCSS = 100.*adjustedR2forNumParams(groupCSS./100, nrTimePoints,1);
        groupCST = 100.*adjustedR2forNumParams(groupCST./100, nrTimePoints,2);
    end
    
    allSubjectsR2 = cat(1,allSubjectsR2,[groupLSS;groupCSS;groupCST]);
    modelNames1   = cat(1,repmat('LSS',size(groupLSS)),repmat('CSS',size(groupCSS)),repmat('CST',size(groupCST)));
    roiNames1     = repmat({char(roisToPlot(ii))},[size(modelNames1,1),1]);
    allROINames   = cat(1,allROINames,roiNames1);
    allModelNames = cat(1,allModelNames,modelNames1);
end

Tanova_full = table(allSubjectsR2,allROINames,allModelNames, 'VariableNames',{'R2','ROI','pRFModel'});
WithinDesignT = table([1 2 3],'VariableNames',{'Models'});
rm = fitrm(Tanova_full,'R2~ROI*pRFModel','WithinDesign',WithinDesignT);
twoWayAnovaResults = rm.anova;
ttestResults = rm.multcompare('pRFModel','By','ROI','ComparisonType','Bonferroni');

grpstat_ROI = grpstats(rm,'ROI','pRFModel',{'mean','std','meanci','var'});

