function fH = makeSupplFigure9_STRetParams_RegressionSlopes_v_Params( ...
    allDS,allLMMResults,allLMMResults_model,cmapModels,roisToPlot,cmapROIs,...
    temporalModel, spatialModel, subjnrs, saveFigs, saveFigDir)
%% Function to reproduce supplementary figure 9:
% Note that Pearson's correlation values may be slightly vary from values
% published (±0.01), due to random sampling of parameter values. 
%
% panel a: Suppression level (LMM  slopes) vs model-based (predicted) suppression
%           level by each pRF model (CSTfix, CSTopt, DN-ST)
% panel b: CSTopt & CST fix Suppression level (LMM  slopes) vs effective pRF size for each ROI
% panel c: CSTopt & CST fix Suppression level (LMM  slopes) vs pRF CST exponent for each ROI
% panel d: CSTopt & CST fix Suppression level (LMM  slopes) vs pRF tau for each ROI
% panel e: CSTopt Sustained & combined transient CST channel beta weights for each ROI
% panel f: DN-ST Suppression level (LMM  slopes) vs effective pRF size for each ROI
% panel g: DN-ST Suppression level (LMM  slopes) vs pRF CST exponent (n) for each ROI
% panel h: DN-ST Suppression level (LMM  slopes) vs pRF tau1,tau2 for each ROI
% panel i: DN-ST Suppression level (LMM  slopes) vs pRF semi-saturation constant for each ROI
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
% - projectDir
% - allDS               : cell (3x) of datasets for each pRF model
% - allLMMResults       : DATA cell (3x number of ROIs), containing a struct with fields:
%                           fixedIntercepts, fixedSlopes, 
%                           fixedIntercepts_CI, fixedSlopes_CI
% - allLMMResults_model : MODELS cell (3x number of ROIs), containing a struct with fields:
%                           fixedIntercepts, fixedSlopes, 
%                           fixedIntercepts_CI, fixedSlopes_CI
% - cmapROIs            : color map for ROIs 
% - roisToPlot          : cell with ROI names
% - cmapModels          : color map for pRF models
% - spatialModel        : spatial components of pRF models
% - temporalModel       : temporal components of pRF models
% - subjnrs             : Subjects to plot
% - saveFigs            : save figures or not?
% - saveFigDir          : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle  
%


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

