function fH = makeFigure8_SuppressionLevels_v_pRFParams(ds,lmmResults,roisToPlot,...
                        cmapROIs,temporalModel, spatialModel,saveFigs, saveFigDir)

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


