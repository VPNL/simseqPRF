function fH = makeFigure6_PRFModelFits(ds,roisToPlot)

saveFigs = false;

%% Panel A: Plot V1 voxel time series + predicted pRF response per condition
projectDir  = fullfile(simseqRootPath);
subjnr      = 3;
fH(1)       = plotDataModelFit_singleVoxel_StimulusConditionBlocks(projectDir,subjnr,...
                'roisToPlot',{'V1'},'selectedDataVoxels',941, 'plotModelFlag',true);

%% Panel B: Plot VO1/2 voxel time series + predicted pRF response per condition
fH(2)       = plotDataModelFit_singleVoxel_StimulusConditionBlocks(projectDir,subjnr,...
    'roisToPlot',{'VO1'},'selectedDataVoxels',220, 'plotModelFlag',true);

%% Get resampled data
[~,~,~,~,~,resampledCVR2, ~, mean_diffCVR2, maxNC] = ...
            resamplePRFParams_wReplacement(ds);
        
%% Panel C: Plot distribution pRF model cross-validated R2 (violin graph)   
fH(3) = plotPRFModelCVRSQ_violinPlot(ds,resampledCVR2,maxNC,roisToPlot,saveFigs);

%% Panel D: Plot pRF model difference cross-validated R2 (bar graph)
fH(4) = plotDiffPRFModelCVR2_bargraph(mean_diffCVR2,roisToPlot,saveFigs);


%% Do some statistics (ANOVA)

% Do we want to adjust for the nr of regressors in pRF models (CST having 1
% more than LSS and CSS).
useAdjustedR2 = false;
nrTimePoints = 648;

% Preallocate arrays
allROINames = {};
allModelNames = [];
allSubjectsR2 = [];

for ii = 1:length(roisToPlot)
    groupLSS = ds.R2LSS(ds.ROI==roisToPlot(ii) & ds.Condition==nominal(1));
    groupCSS = ds.R2CSS(ds.ROI==roisToPlot(ii) & ds.Condition==nominal(1));
    groupCST = ds.R2CST(ds.ROI==roisToPlot(ii) & ds.Condition==nominal(1));

    if useAdjustedR2
        groupLSS = 100.*adjustedR2forNumParams(groupLSS./100, nrTimePoints,1);
        groupCSS = 100.*adjustedR2forNumParams(groupCSS./100, nrTimePoints,1);
        groupCST = 100.*adjustedR2forNumParams(groupCST./100, nrTimePoints,2);
    end
    
    allSubjectsR2 = cat(1,allSubjectsR2,[groupLSS;groupCSS;groupCST]);
    modelNames1   = cat(1,repmat('LSS',size(groupLSS)),repmat('CSS',size(groupCSS)),repmat('CST',size(groupCST)));
    roiNames1     = repmat({char(roisToPlot(ii))},[size(modelNames1,1),1]);
    allROINames   = cat(1,allROINames,roiNames1);
    allModelNames = cat(1,allModelNames,modelNames1);
end

Tanova_full = table(allSubjectsR2,allROINames,allModelNames, 'VariableNames',{'R2','ROI','pRFModel'});
WithinDesignT = table([1 2 3],'VariableNames',{'Models'});
rm_an2 = fitrm(Tanova_full,'R2~ROI*pRFModel','WithinDesign',WithinDesignT);
rm_an2.anova
rm_an2.multcompare('pRFModel','By','ROI','ComparisonType','Bonferroni')

