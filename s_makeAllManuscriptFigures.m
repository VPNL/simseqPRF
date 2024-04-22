%% s_makeAllManuscriptFigures.m

%% Define params
subjnrs          = [1,2,3,7,8,9,10,11,12,13];
projectDir       = fullfile(simseqRootPath);
modelName        = 'suppl_stRet'; % Choose from: 'main' (i.e., 'CST','CSS','LSS'),'suppl_DoG','suppl_stRet';
saveFigs         = true; % Save figures or not
loadLMMResults   = true; % Load presaved results or recompute from scratch
pths             = getSubjectPaths(projectDir,1);
saveFigDir       = fullfile(simseqRootPath,'figures');
saveFolderNameDS = fullfile(projectDir,'simseq/results/group');

% Define colors, labels, plot order
cmapROIs             = getROISummaryColors(0);
conditionNamesSimSeq = {'Small-0.2s','Small-1s','Large-0.2s','Large-1s'};
nrConditions         = length(conditionNamesSimSeq);
conditionOrderSimSeq = [1:4];
conditionNamesSimSeq = conditionNamesSimSeq(conditionOrderSimSeq);
labelsComb           = {'V1','V2','V3','hV4','VO1/2','LO1/2','TO1/2','V3AB','IPS0/1'};
roiOrder             = [1:6,8:11,7,12,13]; %V1-V02, LO1/2 + TO1/2, V3AB, IPS0/1
combROIs             = [5,8,10,12]; %VO1/2 (5,6), LO1/2(8,9), TO1/2 (10,11), IPS0/1 (12,13)
nrows                = 4; % number of rows in figures
ncols                = 9; % number of columns in figures
szMarker             = 5; % size of markers in figures

% Define noise and variance explained thresholds:
nc_thresh            = 0.1; % threshold in percentage noise ceiling from simseq exp
ve_thresh            = 0.2; % threshold in percentage of pRF var expl. from toonotopy

% Set defaults for prettier figure plotting
makeprettyfigures;


switch modelName

    case {'main'}
        if strcmp(modelName,'CST')
            spatialModel   = 'onegaussianFit';
            temporalModel  = '3ch-stLN';
        elseif strcmp(modelName,'CSS')
            spatialModel   = 'cssFit';
            temporalModel  = '1ch-glm';
        elseif strcmp(modelName,'LSS')
            spatialModel   = 'onegaussianFit';
            temporalModel  = '1ch-glm';
        elseif strcmp(modelName,'main')
            modelName = {'CST','CSS','LSS'};
            spatialModel   = {'onegaussianFit','cssFit','onegaussianFit'};
            temporalModel  = {'3ch-stLN','1ch-glm','1ch-glm'};
        end
        useSTRetParams     = false;
        cmapModels          = getColormapPRFModels(0);
    
    case 'suppl_DoG'
        spatialModel       = 'differenceOfGaussiansFit';
        temporalModel      = '1ch-glm';
        useSTRetParams     = false;
        cmapModels          = getColormapPRFModels(2);
    
    case {'suppl_stRet'}
        subjnrs            = [1,2,3,7,8,9,10];
        spatialModel       = {'oneGaussianFit','oneGaussianFit','oneGaussianFit'};
        temporalModel      = {'3ch-stLN','3ch-stLN','1ch-dcts'};
        useSTRetParams     = true;
        cmapModels          = getColormapPRFModels(3);
end

                
%% Load or run creation of data table & LMM on the fly
if loadLMMResults
     if strcmp(modelName,'suppl_stRet')
         d = dir(fullfile(pths.dataDirSimSeq,'group',sprintf('SIMSEQ_dataTable_stRetParams_matchVoxels_*.mat')));
     elseif strcmp(modelName,'main')
         d = dir(fullfile(pths.dataDirSimSeq,'group',sprintf('SIMSEQ_dataTable_20230426.mat')));  
     elseif strcmp(modelName,'suppl_DoG')
         d = dir(fullfile(pths.dataDirSimSeq,'group',sprintf('SIMSEQ_dataTable_differenceOfGaussiansFit_20231128.mat')));
     end
    load(fullfile(d(end).folder,d(end).name))
else    
    if useSTRetParams
        [allDS, allLMMResults, allLMMResults_model] = ...
            getDataTableAndLMMResults_stRet_matchedVoxels(projectDir,spatialModel, temporalModel, useSTRetParams, saveFolderNameDS, postFix);
    else
        lmmModelFlag = true;
        [ds, fLMM, fLMM_Model, lmmResults, lmmResults_Model] = ...
            getDataTableAndLMMResults(projectDir, subjnrs, saveFolderNameDS, postFix, lmmModelFlag, useSTRetParams, spatialModel, temporalModel);
    end
end

% Reorder ROIs
if useSTRetParams
    allRoisToPlot = unique(allDS{1}.ds.ROI,'stable');
else
    allRoisToPlot = unique(ds.ROI,'stable');
end
if strcmp(spatialModel,'differenceOfGaussiansFit')
    roisToPlot = allRoisToPlot;
else
    roiNames = ["V1","V2","V3","hV4","VO1/VO2","V3AB","IPS0/IPS1","LO1/LO2","TO1/TO2"];
    for ii = 1:length(roiNames)
        newROIOrder(ii) = find(ismember(string(allRoisToPlot),roiNames{ii}));
    end
    
    roisToPlot = allRoisToPlot(newROIOrder);
end

if strcmp(modelName,'main')
    %% FIGURE 2:
    subjnr = 3;
    fH2 = makeFigure2_exampleV1Data(ds, fLMM, lmmResults, projectDir,subjnr, saveFigs, saveFigDir);

    %% FIGURE 3:
    subjnr = 3;
    fH3 = makeFigure3_exampleVO12Data(ds, fLMM, lmmResults, projectDir,subjnr, saveFigs, saveFigDir);

    %% FIGURE 4:
    fH4 = makeFigure4_SuppressionLevelsAcrossVisualAreas(ds,lmmResults,roisToPlot,cmapROIs, ...
                                                            useSTRetParams, saveFigs, saveFigDir);

    %% FIGURE 6 & Suppl Table 3:
    fH6 = makeFigure6_PRFModelFits(ds,roisToPlot,spatialModel,temporalModel,saveFigs, saveFigDir);

    %% FIGURE 7:
    fH7 = makeFigure7_PredictedSuppressionLevels(lmmResults, lmmResults_Model,roisToPlot,cmapROIs,saveFigs,saveFigDir);

    %% FIGURE 8:
    fH8 = makeFigure8_SuppressionLevels_v_pRFParams(ds,lmmResults,roisToPlot,...
                        cmapROIs,spatialModel,temporalModel,saveFigs,saveFigDir);

    %% SUPPLEMENTARY FIGURE 1: Sim vs Seq scatter plots: All visual areas & all conditions
    fHS1 = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, lmmResults, ...
        'plotAllSubjectsTogether', true);

    %% SUPPLEMENTARY FIGURE 2: Delta suppresion levels across the visual hierarchy
    fHS2 = makeSupplFigure2_SuppressionDeltas(ds,lmmResults,roisToPlot, ...
                    cmapROIs,subjnrs,useSTRetParams,temporalModel, spatialModel, ...
                    saveFigs,saveFigDir);

    %% SUPPLEMENTARY FIGURE 3: Simulated model predictions (and Fig 5: spatial/temporal filters)
    fHS3 = makeSupplFigure2_SimulateModelPredictions(saveFigs);

    %% SUPPLEMENTARY FIGURE 4: CST exponent distributions
    fHS4 = makeSupplFigure3_DistributionCSTExponents(ds,roisToPlot,cmapROIs,saveFigs);

    %% SUPPLEMENTARY FIGURE 5: Eye tracking
    fHS5 = makeSupplFigure4_eyefixation(projectDir,[],saveFigs);

    %% SUPPLEMENTARY FIGURE 6: LMM comparisons
    fHS6 = makeSupplFigure5_compareLMMs(projectDir, saveFigs);
    
elseif strcmp(modelName,'suppl_DoG')
    %% SUPPLEMENTARY FIGURE 7: Difference of Gaussian model
    fHS7 = makeSupplFigure7_DoGPredictions(projectDir,ds,lmmResults, ...
        lmmResults_Model, modelName, roisToPlot,cmapModels, useSTRetParams, ...
        temporalModel, spatialModel, subjnrs, saveFigs, saveFigDir);
    
elseif strcmp(modelName,'suppl_stRet')
    %% SUPPLEMENTARY FIGURE 8: Spatiotemporal pRF models (time series and cv-R^2)
    fHS8 = makeSupplFigure8_STRetParams_TimeSeriesPredictionsRsquared(projectDir, ...
    allDS,modelName,roisToPlot,temporalModel, spatialModel, subjnrs, saveFigs, saveFigDir);

    %% SUPPLEMENTARY FIGURE 9: Spatiotemporal pRF models (predicted suppression and parameters)
    fHS9 = makeSupplFigure9_STRetParams_RegressionSlopes_v_Params( ...
    allDS, allLMMResults,allLMMResults_model,cmapModels,roisToPlot,cmapROIs,...
    temporalModel, spatialModel, subjnrs, saveFigs, saveFigDir);

end


return
% 
% %% Extra figures
% % TO PLOT INDIVIDUAL SUBJECTS SIM vs SEQ VOXEL AMPLITUDES (SCATTERPLOTS)
% fHX1 = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, lmmResults, ...
%         'plotAllSubjectsTogether', false);
%     
% fHX2 = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM,lmmResults, ...
%         'plotAllSubjectsTogether', true,...
%         'plotModelAmplSubj', true);





