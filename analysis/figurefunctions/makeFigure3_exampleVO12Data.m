function fH = makeFigure3_exampleVO12Data(ds, fLMM, lmmResults, projectDir,subjnr) 

%% Panel A: Plot VO1/2 voxel time series per condition
fH(1) = plotDataModelFit_singleVoxel_StimulusConditionBlocks(projectDir,subjnr,...
    'roisToPlot',{'VO1'},'selectedDataVoxels',220, 'plotModelFlag',false);

%% Panel B: Scatter plot VO1/2 all voxel data, subject S3
fH(2) = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, lmmResults, 'subjnrs',subjnr,'saveFigs',false, 'roisToPlot',5);

return