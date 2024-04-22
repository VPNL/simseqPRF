function fH = plotDNST_CSTopt_CViolinPlots(params,roisToPlot,...
    cmapModels, saveFigs, saveFigDir)

%% CV-R^2 violin Plots
modelName      = {'CST_{fix}','CST_{opt}','DNST'};
cmapModelsCell = {cmapModels(1,:);cmapModels(2,:);cmapModels(3,:)};
mnR2           = NaN(9,3,7);

fH = figure; set(gcf,'Position',[500 1 965 950],'color','w');
plIdx = 1; 

ipsIdx = cellfind(regexp(string(roisToPlot),'IPS0/IPS1','once'));
roisNoIPS = setdiff([1:length(roisToPlot)],ipsIdx);

for idx = roisNoIPS
    clear dataToPlot
    
    dataToPlot{1} = squeeze(params{1}.resampledCVR2(:,idx,:));
    dataToPlot{2} = squeeze(params{2}.resampledCVR2(:,idx,:));
    dataToPlot{3} = squeeze(params{3}.resampledCVR2(:,idx,:));
    
    median_dataToPlot{1} = median(dataToPlot{1}(:),'omitnan');
    median_dataToPlot{2} = median(dataToPlot{2}(:),'omitnan');
    median_dataToPlot{3} = median(dataToPlot{3}(:),'omitnan');
    
    for subj = 1:7
        mnR2(idx,1,subj) = mean(params{1}.R2CST{subj,idx},'omitnan');
        mnR2(idx,2,subj) = mean(params{2}.R2CST{subj,idx},'omitnan');
        mnR2(idx,3,subj) = mean(params{3}.R2DNST{subj,idx},'omitnan');
    end
    
    % Get max noise ceiling across subjects
    tmp = squeeze(params{1}.maxNC(:,idx));
    mdGroupNC = max(tmp(:),[],'omitnan');
    
    ax = subplot(3,3,plIdx); hold all;
    violinPlot(dataToPlot,'color',cmapModelsCell,'xNames',modelName, 'showMM',3)
    title(sprintf('%s md: %2.0f %2.0f %2.0f',string(roisToPlot(idx)),...
        median_dataToPlot{1},median_dataToPlot{2},median_dataToPlot{3}))
    plIdx = plIdx+1;
    plot([0.5 1.5],[mdGroupNC mdGroupNC],'k:');
    plot([1.5 2.5],[mdGroupNC mdGroupNC],'k:');
    plot([2.5 3.4],[mdGroupNC mdGroupNC],'k:');
    ylim([-10 100]); set(gca,'YTick',[0:20:100])
    
    if plIdx ==5
        ylabel('Cross-validated R^2 (%)')
    end
end


if saveFigs
    fName = 'SupplFig8_R2_CSTfix_CSTopt_DNST_stRetParams_ViolinPlots';
    subDir = 'SupplFig8';
    if ~exist('saveFigDir','dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end

%% Do some stats
% Do we want to adjust for the nr of regressors in pRF models (CST having 1
% more than LSS and CSS).
useAdjustedR2 = false;
nrTimePoints = 648;

% Preallocate arrays
allROINames = {};
allModelNames = [];
allvoxelsR2 = [];
allSubjectsR2 = [];
allSubjectNames = [];

roisNoIPS = [1:6,8,9];
for ii = roisNoIPS
    groupCST = cell2mat(params{1}.R2CST(:,ii));
    groupCST_ST = cell2mat(params{2}.R2CST(:,ii));
    groupDNST = cell2mat(params{3}.R2DNST(:,ii));
    
    if useAdjustedR2
        groupCST = 100.*adjustedR2forNumParams(groupCST./100, nrTimePoints,1);
        groupCST_ST = 100.*adjustedR2forNumParams(groupCST_ST./100, nrTimePoints,1);
        groupDNST = 100.*adjustedR2forNumParams(groupDNST./100, nrTimePoints,2);
    end
    
    allvoxelsR2 = cat(1,allvoxelsR2,[groupCST;groupCST_ST;groupDNST]);
    modelNames1   = cat(1,repmat('CSTfix',size(groupCST)),repmat('CSTopt',size(groupCST_ST)),repmat('DSTopt',size(groupDNST)));
    roiNames1     = repmat({char(roisToPlot(ii))},[size(modelNames1,1),1]);
    allROINames   = cat(1,allROINames,roiNames1);
    allModelNames = cat(1,allModelNames,modelNames1);
    allSubjectNames = cat(1,allSubjectNames,repmat([1:7]',[3,1]));
end

% 2-way ANOVA (ROI x Model) across voxels
Tanova = table(allvoxelsR2,allROINames,allModelNames, 'VariableNames',{'R2','ROI','pRFModel'});
withinDesignT = table([1 2 3],'VariableNames',{'Models'});
rm_vox = fitrm(Tanova,'R2~ROI*pRFModel','WithinDesign',withinDesignT);
anovatbl_vox = rm_vox.anova;
rm_vox_ttest = rm_vox.multcompare('pRFModel','By','ROI','ComparisonType','Bonferroni');

