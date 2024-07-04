function fH = makeSupplFigure8_STRetParams_TimeSeriesPredictionsRsquared(projectDir,...
    allDS,roisToPlot,temporalModel, spatialModel, subjnrs, saveFigs, saveFigDir)
% Function to reproduce supplementary manuscript figure 8: 
% panel a/b: V1 time series + CSTfix, CSTopt, DN-ST model prediction for each stimulus condition 
% panel b: violin plots showing cross-validated R^2 for CSTfix, CSTopt, DN-ST pRF models
% panel c: Difference in cv-R^2 between CSTfix, CSTopt, DN-ST models, 
% for all ROIs except IPS0/1 (due to insufficient nr of subjects with data in this ROI)
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
% INPUTS (required):
% - projectDir
% - allDS           : cell (3x) of datasets for each pRF model
% - roisToPlot      : cell with ROI names
% - spatialModel    : spatial components of pRF models
% - temporalModel   : temporal components of pRF models
% - subjnrs         : Subjects to plot
% - saveFigs        : save figures or not?
% - saveFigDir      : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle  
%
% Code written by E.R. Kupers (2024) Stanford University

%% Set params
cmapModels = getColormapPRFModels(3);

% Get resampled data
useSTRetParams = true;

for jj = 1:length(temporalModel)
    output{jj} = resamplePRFParams_wReplacement(allDS{jj}.ds, roisToPlot, useSTRetParams, ...
                    temporalModel{jj},spatialModel{jj}, subjnrs);
end

%% Suppl Fig 8A: Plot single pRF time series
subjnr     = 3;
exampleVox = 217;
fH(1) = plotDataModelFitSTRet_singleVoxel_StimulusConditionBlocks(projectDir,subjnr,...
    'roisToPlot',{'VO1'},'selectedDataVoxels',exampleVox, 'plotModelFlag',true,'saveFigs', true, 'saveFigDir',saveFigDir);

%%  Suppl Fig 8B: Plot violin plots
fH(2) = plotDNST_CSTopt_CViolinPlots(output,roisToPlot,...
                            cmapModels, saveFigs, saveFigDir);
                        
%%  Suppl Fig 8C: Diff cv-R^2 bar graphs
fH(3) = plotDNST_CSTopt_DiffCV_BarPlot(output, roisToPlot, cmapModels, saveFigs, saveFigDir);
end