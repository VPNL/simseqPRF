%% s_makeAllManuscriptFigures.m

%% Define params
subjnrs        = [1,2,3,7,8,9,10,11,12,13];
projectDir     = fullfile(simseqRootPath);
spatialModels  = {'onegaussianFit','cssFit','onegaussianFit'};
temporalModels = {'1ch-glm','1ch-glm','3ch-stLN'};
modelNames     = {'LSS','CSS','CST'};

pths           = getSubjectPaths(projectDir,1);
saveFigDir     = fullfile(simseqRootPath,'figures');
saveFigs       = true; % Save figures or not
loadLMMResults = true; % Load presaved results or recompute from scratch

% Define colors
cmapROIs             = getROISummaryColors(0);

% Set defaults for prettier figure plotting
makeprettyfigures;

%% Load or run creation of data table & LMM on the fly
if loadLMMResults
    d = dir(fullfile(pths.dataDirSimSeq,'group',sprintf('SIMSEQ_dataTable_*.mat')));
    load(fullfile(d(end).folder,d(end).name))
    
     % Reorder ROIs
    roisToPlot = unique(ds.ROI,'stable');
    newROIOrder = [1,2,3,4,5,8,9,6,7];
    roisToPlot = roisToPlot(newROIOrder);
    cmapROIs = cmapROIs(newROIOrder,:);
else    
    saveFolderNameDS = fullfile(projectDir,'simseq/results/group');
    lmmModelFlag = true;
    [ds, fLMM, fLMM_Model] = getDataTableAndLMMResults(projectDir, lmmModelFlag, saveFolderNameDS)
end


%% FIGURE 2:
subjnr = 3;
fH2 = makeFigure2_exampleV1Data(ds, fLMM, lmmResults, projectDir,subjnr); 

%% FIGURE 3:
subjnr = 3;
fH3 = makeFigure3_exampleVO12Data(ds, fLMM, lmmResults, projectDir,subjnr);

%% FIGURE 4:
fH4 = makeFigure4_SuppressionLevelsAcrossVisualAreas(ds,fLMM,lmmResults,roisToPlot,cmapROIs);

%% FIGURE 6 & Suppl Table 3:
fH6 = makeFigure6_PRFModelFits(ds,roisToPlot);

%% FIGURE 7:
fH7 = makeFigure7_PredictedSuppressionLevels(lmmResults, lmmResults_Model,roisToPlot,cmapROIs);

%% FIGURE 8:
fH8 = makeFigure8_SuppressionLevels_v_pRFParams(ds,lmmResults,roisToPlot,cmapROIs,saveFigs);

%% SUPPLEMENTARY FIGURE 1: Sim vs Seq scatter plots: All visual areas & all conditions
fHS1 = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, lmmResults, ...
    'plotAllSubjectsTogether', true);

%% SUPPLEMENTARY FIGURE 2: Simulated model predictions (and Fig 5: spatial/temporal filters)
fHS2 = makeSupplFigure2_SimulateModelPredictions(saveFigs);

%% SUPPLEMENTARY FIGURE 3: CST exponent distributions
fHS3 = makeSupplFigure3_DistributionCSTExponents(ds,roisToPlot,cmapROIs,saveFigs);

%% SUPPLEMENTARY FIGURE 4: Eye tracking
fHS4 = makeSupplFigure4_eyefixation(projectDir,[],saveFigs);

%% SUPPLEMENTARY FIGURE 5: LMM comparisons
fHS5 = makeSupplFigure5_compareLMMs(projectDir, saveFigs);

%% Extra figures
% TO PLOT INDIVIDUAL SUBJECTS SIM vs SEQ VOXEL AMPLITUDES (SCATTERPLOTS)
fHX1 = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, lmmResults, ...
        'plotAllSubjectsTogether', false);
    
fHX2 = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM,lmmResults, ...
        'plotAllSubjectsTogether', true,...
        'plotModelAmplSubj', true);





