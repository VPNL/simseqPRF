function fH = makeFigure8_SuppressionLevels_v_pRFParams(ds,lmmResults,roisToPlot,...
                        cmapROIs,temporalModel, spatialModel,saveFigs, saveFigDir)
% Function to reproduce main manuscript figure 8: 
% panel a: Suppression level (LMM  slopes) vs effective pRF size for each ROI
% panel b: Suppression level (LMM  slopes) vs pRF CST exponent for each ROI
% panel c: Sustained & combined transient CST channel beta weights for each ROI
% panel d: CST vs CSS pRF exponent for each ROI
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
% - roisToPlot       : cell with ROI names
% - cmapROIs         : color map for ROIs
% - spatialModel     : spatial components of pRF models
% - temporalModel    : temporal components of pRF models
% - saveFigs         : save figures or not?
% - saveFigDir       : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle
%
%% Get resampled data
subjnrs = [1:3,7:13];
useSTRetParams = false;
output = resamplePRFParams_wReplacement(ds, roisToPlot, useSTRetParams, temporalModel, spatialModel, subjnrs);

%% Panel A: Plot Suppression slopes vs CST pRF size
fH(1) = plotLMMregressionSlopes_v_pRFSize(lmmResults, ...
        output.median_resampledPRFSz, roisToPlot, cmapROIs, saveFigs,saveFigDir);

%% Panel B: Plot Suppression slopes vs CST pRF exponent
fH(2) = plotLMMRegressionSlopes_v_pRFCSTExp(lmmResults,...
        output.median_resampledCSTExp,roisToPlot, cmapROIs, saveFigs,saveFigDir);

%% Panel C:Plot CST pRF exponent vs CSS pRF exponent
fH(3) = plotMeanCSTExponent_v_MeanCSSExponent(output.median_resampledCSTExp,...
        output.median_resampledCSSExp,roisToPlot,cmapROIs,saveFigs,saveFigDir);

%% Panel D:Plot CST pRF model sustained and transient beta weights 
fH(4) = plotMeanCSTSustainedTransientBetaWeights(output.mean_resampledBetavalSust, ...
        output.mean_resampledBetavalTrans,roisToPlot,saveFigs,saveFigDir);


