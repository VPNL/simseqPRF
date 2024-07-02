function fH = makeFigure2_exampleV1Data(ds, fLMM, lmmResults, projectDir,subjnr, saveFigs, saveFigDir) 
% Function to reproduce main manuscript figure 2: 
% panel a: V1 time series for each stimulus condition 
% panel b: V1 SEQ vs SIM amplitude voxel scatter plot for each stimulus condition
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
% - ds         : dataset
% - lmmResults : cell (1x number of ROIs), containing a struct with fields:
%                       fixedIntercepts, fixedSlopes, 
%                       fixedIntercepts_CI, fixedSlopes_CI
% - projectDir  : base folder of project
% - subjnr      : Subject nr to plot
% - saveFigs    : save figures or not?
% - saveFigDir  : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle

%% Panel A: Plot V1 voxel time series per condition
fH(1) = plotDataModelFit_singleVoxel_StimulusConditionBlocks(projectDir,subjnr,...
    'roisToPlot',{'V1'},'selectedDataVoxels',941, 'plotModelFlag',false, ...
    'saveFigDir',saveFigDir, 'saveFigs', saveFigs);

%% Panel B: Scatter plot V1 all voxel data, subject S3
fH(2) = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, lmmResults, ...
    'subjnrs',subjnr,'roisToPlot',1,'saveFigDir',saveFigDir,'saveFigs',false);

return