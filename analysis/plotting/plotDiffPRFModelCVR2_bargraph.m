function fH = plotDiffPRFModelCVR2_bargraph(mnDiffR2Subj,roisToPlot,saveFigs, saveFigDir)

% Set up figure params
cmapModels = getColormapPRFModels(0); % Model colors
barWidth   = 0.75;
xposBar    = [1:length(roisToPlot)];
xposError  = [-0.215,0,0.215]; % X position for errorbars
modelDiffNames = {'CSS-LSS','CST-LSS','CST-CSS'};

% Change order of models, start with CSS first
newOrder = [2,3,1]; 
mnDiffR2SubjNewOrder   = mnDiffR2Subj(:,:,newOrder);
modelDiffNamesNewOrder = modelDiffNames(newOrder);
cmapModelsNewOrder     = cmapModels(newOrder,:);

% Summarize group data (mean and SEM)
for idx = 1:length(roisToPlot)
    for mm = 1:3
        nrSubjWithData        = sum(~isnan(mnDiffR2SubjNewOrder(:,idx,mm)));
        mnDiffR2Group(idx,mm) = mean(mnDiffR2SubjNewOrder(:,idx,mm),1,'omitnan');
        seDiffR2Group(idx,mm) = std(mnDiffR2SubjNewOrder(:,idx,mm),[],1,'omitnan')./sqrt(nrSubjWithData);
    end
end

% Start plotting
fH = figure; clf; set(gcf,'Position',[633    84   558   657],'color','w'); hold all;
b = bar(xposBar,mnDiffR2Group, barWidth); hold on;
b(1).FaceColor = cmapModels(3,:);
b(2).FaceColor = cmapModels(3,:);
b(3).FaceColor = cmapModels(2,:);

b(1).EdgeColor = cmapModels(1,:);
b(2).EdgeColor = cmapModels(2,:);
b(3).EdgeColor = cmapModels(1,:);
b(1).LineWidth = 5;
b(2).LineWidth = 5;
b(3).LineWidth = 5;

% Plot errorbars
errorbar2(xposBar'+xposError,mnDiffR2Group, seDiffR2Group, 1,'k'); hold on;

% Add individual subjects
for idx = 1:length(roisToPlot)
    for mm = 1:3
        currIndivSubj = mnDiffR2SubjNewOrder(:,idx,mm);
        scatter(idx+xposError(mm)+(0.02*randn(1,size(mnDiffR2SubjNewOrder,1)))',currIndivSubj,30,cmapModelsNewOrder(mm,:), 'filled', ...
            'MarkerEdgeAlpha',0.7, 'MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0 0 0]); hold on;
    end
end

% Add axis labels, ticks, legend, and title
xlim([0 10]); ylim([-7 20]); box off
set(gca,'XTick',[1:length(roisToPlot)], 'XTickLabel',string(roisToPlot),'XTickLabelRotation',45)
l = gca;
legend(l.Children([length(l.Children):-1:1]),modelDiffNamesNewOrder, 'FontSize',9, 'Location','NorthEast'); 
legend boxoff
sgtitle('Group Average cv-R^2 mean +/- SEM per visual area')
ylabel('Cross-validated variance explained (%)')

% Save if requested
if saveFigs
    if ~exist('saveFigDir','var')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir, 'fig6');
    fName = sprintf('SummaryGroupR2_gridFit2_diff_CST');   
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose');
end

