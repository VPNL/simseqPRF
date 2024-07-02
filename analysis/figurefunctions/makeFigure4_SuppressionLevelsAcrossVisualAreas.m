function fH = makeFigure4_SuppressionLevelsAcrossVisualAreas(ds,lmmResults, ...
            roisToPlot, cmapROIs, useSTRetParams, saveFigs, saveFigDir)
% Function to reproduce main manuscript figure 4: 
% panel a: All 9 ROIs SEQ vs SIM amplitude voxel scatter plot for small & short stimulus condition
% panel b: Suppression level (LMM  slopes) vs ROI vs stimulus condition
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
% - ds              : dataset
% - lmmResults      : cell (1x number of ROIs), containing a struct with fields:
%                       fixedIntercepts, fixedSlopes, 
%                       fixedIntercepts_CI, fixedSlopes_CI
% - roisToPlot      : cell with ROI names
% - cmapROIs        : color map for ROIs
% - useSTRetParams  : (boolean) are we using supplementary spatiotemporal
%                               retinotopy data or not?
% - saveFigs        : save figures or not?
% - saveFigDir      : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle

%% Check inputs
if isempty(cmapROIs) || ~exist('cmapROIs','var')
    cmapROIs = getROISummaryColors(0);
end

if isempty(useSTRetParams) || ~exist('useSTRetParams','var')
    useSTRetParams = false;
end

%% Panel A: PLOT ALL SUBJECTS DATA AND FIT IN ONE FIGURE
fH(1) = plotMeanSeqVsSimAmplitude_voxel(ds, lmmResults,...
    'plotAllSubjectsTogether', true, ...
    'saveFigs', saveFigs, ...
    'saveFigDir', saveFigDir);

%% Plot LMM fitted regression slopes
LMMOrder = {'Data'};
fH(2) = plotLMMfittedRegressionSlopes(ds,lmmResults,LMMOrder, ...
        roisToPlot, cmapROIs, useSTRetParams, saveFigs,saveFigDir);

%% Do some stats (ANOVA and post-hoc multiple comparisons)

nrSubjects   = length(unique(ds.Subject));
nrROIs       = length(roisToPlot);
conditionNamesSimSeq = {'Small & Short (0.2s)','Small & long (1s)', ...
                        'Big & Short (0.2s)','Big & Long (1s)'};
nrConditions = length(conditionNamesSimSeq);

% Accumulate all slopes in one array
allROISlopes = NaN(nrSubjects,nrConditions,nrROIs); 
for ii = 1:nrROIs
    subjIdx = unique(ds.Subject(ds.ROI==roisToPlot(ii) & ds.Condition==nominal(1)));
    allROISlopes(subjIdx,:,ii) = lmmResults.subjSlopes{ii}; 
end

subjectT   = repmat([1:nrSubjects],[1,nrConditions*nrROIs]);
conditionT = repmat(conditionNamesSimSeq,[nrSubjects,nrROIs,1]);
roiT       = repmat(string(roisToPlot)',[nrSubjects*nrConditions,1]);

% Create table for repeated measures ANOVA (subjects x condition x ROI) 
T_ANOVA_slopeXroi = table(allROISlopes(:),conditionT(:),roiT(:),subjectT(:),...
                                       'VariableNames',{'Slope','Condition','ROI','Subject'});
WithinSubjectsT = table([1:10],'VariableNames',{'Subjects'});
rmANOVA_slopeXroi_within = fitrm(T_ANOVA_slopeXroi,'Slope~Condition*ROI',...
                'WithinDesign',WithinSubjectsT,'WithinModel','separatemeans');
rm_anova  = rmANOVA_slopeXroi_within.anova
rm_multcp = rmANOVA_slopeXroi_within.multcompare('Condition','By','ROI','ComparisonType','Bonferroni')

% grpstat_ROI = grpstats(rmANOVA_slopeXroi_within,'ROI',{'mean','std','meanci','var'});

