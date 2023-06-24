function fH = makeFigure2_exampleV1Data(ds, fLMM, lmmResults, projectDir,subjnr) 

%% Panel A: Plot V1 voxel time series per condition
fH(1) = plotDataModelFit_singleVoxel_StimulusConditionBlocks(projectDir,subjnr,...
    'roisToPlot',{'V1'},'selectedDataVoxels',941, 'plotModelFlag',false);

%% Panel B: Scatter plot V1 all voxel data, subject S3
fH(2) = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, lmmResults, 'subjnrs',subjnr,'saveFigs',false, 'roisToPlot',1);

return