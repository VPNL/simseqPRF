function fH = makeSupplFigure7_DoGPredictions(projectDir,ds,lmmResults, ...
                lmmResults_Model,modelName, roisToPlot,cmapModels, useSTRetParams, ...
                temporalModel, spatialModel, subjnrs, saveFigs, saveFigDir)
% Check inputs
if isempty(cmapModels) || ~exist('cmapModels','var')
    cmapModels = getColormapPRFModels(2);
end

if isempty(useSTRetParams) || ~exist('useSTRetParams','var')
    useSTRetParams = false;
end

%% Panel A: Predicted and observed time series (alike Figure 6A)
subjnr = 3;
fH(1) = plotDataModelFitDoG_singleVoxel_StimulusConditionBlocks(...
    projectDir,subjnr,...
    'roisToPlot',{'V1'}, ...
    'selectedDataVoxels',[313],'plotModelFlag',true, ...
    'spatialModel',{'onegaussianFit','cssFit','onegaussianFit','differenceOfGaussiansFit'},...
    'temporalModel',{'1ch-glm','1ch-glm','3ch-stLN','1ch-glm'}, ...
    'mdllbl',{'LSS','CSS','CST','DoG'}, ...
    'saveFigs',saveFigs, ...
    'saveFigDir',fullfile(saveFigDir,sprintf('subj%02d',subjnr)));


%% Panel B: Predicted and observed regression slopes (alike Figure 7)
all_lmmResults = cat(1,lmmResults_Model,lmmResults);
LMMOrder       = {modelName,'Data'};
fH(2)          = plotLMMfittedRegressionSlopes(ds, all_lmmResults,LMMOrder, ...
                    roisToPlot,cmapModels,useSTRetParams, saveFigs, saveFigDir);

%% Panel C: Model cv-R^2 (alike Figure 6C)

% Get resampled data
out = resamplePRFParams_wReplacement(ds, roisToPlot, useSTRetParams, temporalModel, spatialModel, subjnrs);

fH(3) = figure; clf; set(gcf,'Position',[121 631 965 273],'color','w');

% % Get max noise ceiling across subjects
% mdGroupNC = max(out.NC, [], 'omitnan');

for idx = 1:4
    clear dataToPlot
  
    dataToPlot{1} = squeeze(out.resampledCVR2(:,idx,:));
    median_dataToPlot = median(dataToPlot{1}(:),'omitnan');
    
    if ~isempty(dataToPlot)
        ax = subplot(1,4,idx); hold all;
        violinPlot(dataToPlot,'color',cmapModels,'xNames',modelName, 'showMM',3)
        title(sprintf('%s md: %2.0f',string(roisToPlot(idx)),median_dataToPlot))
        ylim([-10 100]); set(gca,'YTick',[0:20:100])
        
    end
    if idx ==1
        ylabel('Cross-validated R^2 (%)')
    end
end

if saveFigs
    fName = 'R2Group_DoG_ViolinPlotAllVoxels';
    subDir = 'SupplFig2';
    if ~exist('saveFigDir','dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir,subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose')
end


end