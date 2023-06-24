function fH = makeFigure4_SuppressionLevelsAcrossVisualAreas(ds,fLMM,lmmResults,roisToPlot,cmapROIs)


%% Panel A: PLOT ALL SUBJECTS DATA AND FIT IN ONE FIGURE
fH(1) = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, lmmResults,...
    'fixedIntercepts', lmmResults.fixedIntercepts, ...
    'fixedSlopes', lmmResults.fixedSlopes, ...
    'fixedIntercepts_CI', lmmResults.fixedIntercepts_CI, ...
    'fixedSlopes_CI', lmmResults.fixedSlopes_CI, ...
    'plotAllSubjectsTogether', true);

%% Plot LMM fitted regression slopes
LMMOrder = {'Data'};
saveFigs = false;
fH(2) = plotLMMfittedRegressionSlopes(lmmResults,LMMOrder,roisToPlot,cmapROIs, saveFigs);

%% Do some stats (ANOVA and post-hoc multiple comparisons)

nrSubjects   = 10;
nrConditions = 4;
nrROIs       = length(roisToPlot);
conditionNamesSimSeq = {'Small & short (0.2s)','Small & long (1s)', ...
                        'Large & Short (0.2s)','Large & Long (1s)'};
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
