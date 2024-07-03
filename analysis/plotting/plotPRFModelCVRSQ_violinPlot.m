function fH = plotPRFModelCVRSQ_violinPlot(ds,resampledCVR2,maxNC,roisToPlot,saveFigs, saveFigDir)
% Function to generate Violin plots with cross-validated variance explained
% (R^2) for each pRF model,

% Set colors
cmapModels = getColormapPRFModels(0);
modelNames = {'LSS','CSS','CST'};

fH = figure; clf; set(gcf,'Position',[500 1 965 950],'color','w');
plIdx = 1;

% Get max noise ceiling across subjects
mdGroupNC = max(maxNC, [], 'omitnan').*100;

modelOrder = [3,2,1]; % CST, CSS, LSS
cmapModels = cmapModels(modelOrder,:);
modelNames = modelNames(modelOrder);
for idx = 1:length(roisToPlot)
    clear dataToPlot
  
    for mm = 1:length(modelOrder)
        dataToPlot{mm} = squeeze(resampledCVR2(:,idx,modelOrder(mm),:)); 
        median_dataToPlot{mm} = median(dataToPlot{mm}(:),'omitnan');
    end
    
    if ~isempty(dataToPlot)
        ax = subplot(3,3,plIdx); hold all;
        violinPlot(dataToPlot,'color',cmapModels(1,:),'xNames',modelNames, 'showMM',3)
        title(sprintf('%s md: %2.0f %2.0f %2.0f',string(roisToPlot(idx)),...
            median_dataToPlot{1},median_dataToPlot{2},median_dataToPlot{3}))
        plIdx = plIdx+1;
        ax.Children(3).FaceColor = cmapModels(2,:);
        ax.Children(2).FaceColor = cmapModels(3,:);
        plot([0.5 3.5],[mdGroupNC(idx) mdGroupNC(idx)],'k:');
        ylim([-10 100]); set(gca,'YTick',[0:20:100])
        
    end
    if plIdx ==5
        ylabel('Cross-validated R^2 (%)')
    end
end

if saveFigs
    if ~exist('saveFigDir','var')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir, 'fig6');
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    fName = 'R2Group_Gridfit2_ViolinPlotAllVoxels';
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2')
end

%% STATS
%%
useAdjustedR2 = false;

% Preallocate space
allROINames    = {};
allGroupR2     = [];
allGroupR2_col = [];
allROIIndx_col = [];
allModelNames  = [];
allROIIndx     = [];

% Loop over ROIs,
for ii = 1:length(roisToPlot)
    groupCVR2LSS = ds.R2LSS(ds.ROI==roisToPlot(ii) & ds.Condition==nominal(1));
    groupCVR2CSS = ds.R2CSS(ds.ROI==roisToPlot(ii) & ds.Condition==nominal(1));
    groupCVR2CST = ds.R2CST(ds.ROI==roisToPlot(ii) & ds.Condition==nominal(1));

    if useAdjustedR2 % Inputs: R2, numDataPoints, numExplanatoryParams
        groupCVR2LSS = 100.*adjustedR2forNumParams(groupCVR2LSS./100, 648,1);
        groupCVR2CSS = 100.*adjustedR2forNumParams(groupCVR2CSS./100, 648,1);
        groupCVR2CST = 100.*adjustedR2forNumParams(groupCVR2CST./100, 648,2);
    end
    
    allGroupR2     = cat(1,allGroupR2,[groupCVR2LSS;groupCVR2CSS;groupCVR2CST]);
    allGroupR2_col = cat(1,allGroupR2_col,[groupCVR2LSS,groupCVR2CSS,groupCVR2CST]);
    allROIIndx_col = cat(1,allROIIndx_col,repmat(ii-1,[size(groupCVR2LSS,1),1]));

    modelNames1   = cat(1,repmat('LSS',size(groupCVR2LSS)),...
                          repmat('CSS',size(groupCVR2CSS)),...
                          repmat('CST',size(groupCVR2CST)));
    roiZeroIdx    = repmat(ii-1,[size(modelNames1,1),1]);
    roiNames1     = repmat({char(roisToPlot(ii))},[size(modelNames1,1),1]);
    allROINames   = cat(1,allROINames,roiNames1);
    allModelNames = cat(1,allModelNames,modelNames1);
    allROIIndx    = cat(1,allROIIndx,roiZeroIdx);
    
end
 
Tanova_full = table(allGroupR2,allROINames,allModelNames, 'VariableNames',{'R2','ROI','pRFModel'});
t2 = table([1 2 3]','VariableNames',{'Models'});
rm_an2 = fitrm(Tanova_full,'R2~ROI*pRFModel');
rm_an2.anova
rm_an2.multcompare('pRFModel','By','ROI','ComparisonType','Bonferroni')

%% Other statistical tests for comparison:
% 2-WAY repeated measures ANOVAN
% [p1,t1,stats1] = anovan(allGroupR2,{allModelNames,allROIIndx},'varnames',{'pRF_model','ROI idx'}, ...
%     'model','full');
% [results1,means1] = multcompare(stats1,'CType','bonferroni');

% % LMM Anova
% T = dataset2table(ds);
% Tanova = T(T.Condition=="1",:);
% Tanova = Tanova(:,[1,2,10,11,12]);
% t = table([1 2 3]','VariableNames',{'Models'});
% rm_an = fitrm(Tanova,'R2LSS,R2CSS,R2CST~ROI','WithinDesign',t,'WithinModel','separatemeans');

