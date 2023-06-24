function fH = makeFigure8_SuppressionLevels_v_pRFParams(ds,lmmResults,roisToPlot,cmapROIs,saveFigs)

% Get resampled data
[median_resampledPRFSz,median_resampledCSSExp,median_resampledCSTExp, ...
          mean_resampledBetavalSust, mean_resampledBetavalTrans] = ...
            resamplePRFParams_wReplacement(ds);

%% Panel A: Plot Suppression slopes vs CST pRF size
fH(1) = plotLMMregressionSlopes_v_pRFSize(lmmResults, ...
    median_resampledPRFSz, roisToPlot, cmapROIs, saveFigs);

%% Panel B: Plot Suppression slopes vs CST pRF exp
fH(2) = plotLMMRegressionSlopes_v_pRFCSTExp(lmmResults,...
    median_resampledCSTExp,roisToPlot, cmapROIs, saveFigs);

%% Panel C:Plot CST pRF exp vs CSS pRF exp
fH(3) = plotMeanCSTExponent_v_MeanCSSExponent(median_resampledCSTExp,...
    median_resampledCSSExp,roisToPlot,cmapROIs,saveFigs);

%% Panel D:Plot CST pRF model sustained and transient beta weights 
fH(4) = plotMeanCSTSustainedTransientBetaWeights(mean_resampledBetavalSust, ...
    mean_resampledBetavalTrans,roisToPlot,saveFigs);


