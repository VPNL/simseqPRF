function fH = makeFigure3_exampleVO12Data(ds, fLMM, lmmResults, projectDir,subjnr, saveFigs, saveFigDir)
% Function to reproduce main manuscript figure 3: 
% panel a: VO1 time series for each stimulus condition 
% panel b: VO1 SEQ vs SIM amplitude voxel scatter plot for each stimulus condition  
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
% - fLMM       : (struct) fitted Linear Mixed Model results for each ROI
% - lmmResults : cell (1x number of ROIs), containing a struct with fields:
%                       fixedIntercepts, fixedSlopes, 
%                       fixedIntercepts_CI, fixedSlopes_CI
% - projectDir  : base folder of project
% - subjnr      : Subject nr to plot
% - saveFigs    : save figures or not?
% - saveFigDir  : folder to save figures

%% Panel A: Plot VO1/2 voxel time series per condition
fH(1) = plotDataModelFit_singleVoxel_StimulusConditionBlocks(projectDir,subjnr,...
    'roisToPlot',{'VO1'},'selectedDataVoxels',220, 'plotModelFlag',false, ...
    'saveFigDir',saveFigDir, 'saveFigs', saveFigs);

%% Panel B: Scatter plot VO1/2 all voxel data, subject S3
fH(2) = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, lmmResults, ...
    'subjnrs',subjnr,'roisToPlot',5, ...
    'saveFigDir',saveFigDir, 'saveFigs', saveFigs);

return