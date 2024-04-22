function fH = makeSupplFigure8_STRetParams_TimeSeriesPredictionsRsquared(projectDir,...
    allDS,modelName,roisToPlot,temporalModel, spatialModel, subjnrs, saveFigs, saveFigDir)

%%% Set params
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