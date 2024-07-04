function fH = plotMeanCSTExponent_v_MeanCSSExponent(...
    resampled_medianCSTExp,resampled_medianCSSExp,roisToPlot,cmapROIs,saveFigs,saveFigDir)

fH = figure; clf; set(gcf,'Position',[889,294,963,604],'color','w');

for idx = 1:length(roisToPlot)

    groupMeanOfMedianCSTExp = mean(resampled_medianCSTExp(:,idx),1,'omitnan'); 
    groupSEMOfMedianCSTExp  = std(resampled_medianCSTExp(:,idx),[],1,'omitnan')./sqrt(sum(~isnan(resampled_medianCSTExp(:,idx))));
    groupMeanOfMedianCSSExp = mean(resampled_medianCSSExp(:,idx),1,'omitnan'); 
    groupSEMOfMedianCSSExp  = std(resampled_medianCSSExp(:,idx),[],1,'omitnan')./sqrt(sum(~isnan(resampled_medianCSSExp(:,idx))));
  
    plot([0, 1], [0, 1], 'k:');   hold all;
    scatter(groupMeanOfMedianCSTExp,groupMeanOfMedianCSSExp,100,cmapROIs(idx,:),'filled','MarkerEdgeColor',[0 0 0]);
    plot([groupMeanOfMedianCSTExp,groupMeanOfMedianCSTExp],groupMeanOfMedianCSSExp+[-groupSEMOfMedianCSSExp, groupSEMOfMedianCSSExp],'color',cmapROIs(idx,:));
    plot(groupMeanOfMedianCSTExp+[-groupSEMOfMedianCSTExp,groupSEMOfMedianCSTExp],[groupMeanOfMedianCSSExp,groupMeanOfMedianCSSExp],'color',cmapROIs(idx,:));
    
    box off; axis square
    ylabel('median CSS pRF exponent'); xlabel('median CST pRF exponent')
    xlim([0, 1]); ylim([0 1]);
end
l = gca;
legend(l.Children([length(l.Children)-1:-4:1]),string(roisToPlot), 'FontSize',9, 'Location','None', 'Position',[0.8 0.5 0.1 0.10],'FontSize',10); legend boxoff
sgtitle(sprintf('Group average (+/- SEM) (N=%d)', size(resampled_medianCSTExp,1)));

if saveFigs
    if ~exist(saveFigDir,'dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    subDir     = 'fig8';
    fName = sprintf('SummaryGroupAverageSEM_MedianCSSvsCSTExpPRF');
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end