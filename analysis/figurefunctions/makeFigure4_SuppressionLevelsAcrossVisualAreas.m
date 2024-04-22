function fH = makeFigure4_SuppressionLevelsAcrossVisualAreas(ds,lmmResults, ...
            roisToPlot, cmapROIs, useSTRetParams, saveFigs, saveFigDir)

% Check inputs
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
rmANOVA_slopeXroi_within.anova
rmANOVA_slopeXroi_within.multcompare('Condition','By','ROI','ComparisonType','Bonferroni')
