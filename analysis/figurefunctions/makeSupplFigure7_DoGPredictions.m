function fH = makeSupplFigure7_DoGPredictions(projectDir,ds,lmmResults, ...
                lmmResults_Model,modelName, roisToPlot,cmapModels, useSTRetParams, ...
                temporalModel, spatialModel, subjnrs, saveFigs, saveFigDir)
% Function to reproduce supplementary manuscript figure 7: 
% panel a/b: V1 time series + LSS, CSS, CST, DOG model prediction for each stimulus condition 
% panel b: violin plots showing cross-validated R^2 for LSS, CSS, CST models
% panel c: Data and modelbased (predicted) suppression slopes for V1-hV4
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
% - projectDir
% - ds              : dataset
% - lmmResults      : cell (1x number of ROIs), containing a struct with fields:
%                       fixedIntercepts, fixedSlopes, 
%                       fixedIntercepts_CI, fixedSlopes_CI
% - lmmResults_Model : MODELS cell (3x number of ROIs), containing a struct with fields:
%                       fixedIntercepts, fixedSlopes, 
%                       fixedIntercepts_CI, fixedSlopes_CI
% - modelName       : Abbreviated names of pRF models
% - roisToPlot      : cell with ROI names
% - cmapModels      : color map for pRF models
% - useSTRetParams  : (boolean) are we using supplementary spatiotemporal
%                               retinotopy data or not?
% - spatialModel    : spatial components of pRF models
% - temporalModel   : temporal components of pRF models
% - subjnrs         : Subjects to plot
% - saveFigs        : save figures or not?
% - saveFigDir      : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle            
            
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
LMMOrder       = {modelName{1},'Data'};
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
    fName = 'SupplFig7_R2Group_DoG_ViolinPlotAllVoxels';
    subDir = 'SupplFig7';
    if ~exist('saveFigDir','dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir,subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%      print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose')
end


end