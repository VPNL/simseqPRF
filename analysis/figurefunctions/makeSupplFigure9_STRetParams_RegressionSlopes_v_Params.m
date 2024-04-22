function fH = makeSupplFigure9_STRetParams_RegressionSlopes_v_Params( ...
    allDS,allLMMResults,allLMMResults_model,cmapModels,roisToPlot,cmapROIs,...
    temporalModel, spatialModel, subjnrs, saveFigs, saveFigDir)
%% Function to make supplementary figure 9 of the paper:
% 
%
%
% Note that Pearson's correlation values may be slightly vary from values
% published (±0.01), due to random sampling of parameter values. 

%% Plot predicted regression slopes (alike Figure 7)
LMMOrder       = {'CST_{fix}','CST_{opt}','DN_ST'};
fH(1) = plotLMMfittedRegressionSlopes_stRetParams(allDS,allLMMResults,allLMMResults_model,LMMOrder,...
                roisToPlot,cmapModels, saveFigs, saveFigDir);
            
%% Plot suppression slopes against average pRF size / exponent / time window (alike Figure 8)

% Get resampled data
useSTRetParams = true;
clear output
for jj = 1:length(temporalModel)
    output{jj} = resamplePRFParams_wReplacement(allDS{jj}.ds, roisToPlot, useSTRetParams, ...
                    temporalModel{jj},spatialModel{jj}, subjnrs);
end

%% Panel B: Plot Suppression slopes vs CST pRF size
fH(2) = plotLMMregressionSlopes_v_pRFSize_stRetParams(allLMMResults, ...
    output, 'CST_ST', roisToPlot, cmapROIs, saveFigs, saveFigDir);

%% Panel C: Plot Suppression slopes vs CST pRF exp
fH(3) = plotLMMregressionSlopes_v_pRFExp_stRetParams(allLMMResults,...
    output,'CST_ST', roisToPlot, cmapROIs, saveFigs,saveFigDir);

%% Panel D: Plot Suppression slopes vs CST pRF tau
fH(4) = plotLMMregressionSlopes_v_pRFtau_stRetParams(allLMMResults,...
    output,'CST_ST', roisToPlot, cmapROIs, saveFigs,saveFigDir);

%% Panel E:Plot CST pRF model sustained and transient beta weights
fH(5) = plotMeanCSTSustainedTransientBetaWeights_stRetParams(allDS, ...
    roisToPlot, saveFigs, saveFigDir);

%% Panel F: Plot Suppression slopes vs DNST pRF size
fH(6) = plotLMMregressionSlopes_v_pRFSize_stRetParams(allLMMResults, ...
    output, 'DN_ST', roisToPlot, cmapROIs, saveFigs, saveFigDir);
   
%% Panel G: Plot Suppression slopes vs DNST pRF exponent
fH(7) = plotLMMregressionSlopes_v_pRFExp_stRetParams(allLMMResults, ...
    output, 'DN_ST', roisToPlot, cmapROIs, saveFigs, saveFigDir);
   
%% Panel H: Plot Suppression slopes vs DNST pRF tau
fH(8) = plotLMMregressionSlopes_v_pRFtau_stRetParams(allLMMResults, ...
    output, 'DN_ST', roisToPlot, cmapROIs, saveFigs, saveFigDir);
   
%% Panel I: Plot Suppression slopes vs DNST pRF semisaturation constant
fH(9) = plotLMMregressionSlopes_v_pRFsemisat_stRetParams(allLMMResults, ...
    output, roisToPlot, cmapROIs, saveFigs, saveFigDir);

end

