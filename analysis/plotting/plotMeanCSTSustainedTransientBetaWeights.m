function fH = plotMeanCSTSustainedTransientBetaWeights(mean_resampledBetavalSust, ...
    mean_resampledBetavalTrans,roisToPlot,saveFigs)

% Set colors
roiColors = getROISummaryColors(0);
newROIOrder = [1,2,3,4,5,8,9,6,7];
roiColors = roiColors(newROIOrder,:);
white = [1 1 1];

% Compute grand mean + SEM across subjects
for rn = 1:length(roisToPlot)
    grandMn_S(rn) = squeeze(mean(mean_resampledBetavalSust(:,rn),1,'omitnan'));
    grandSE_S(rn) = squeeze(std(mean_resampledBetavalSust(:,rn),[],1,'omitnan'))./sqrt(sum(~isnan(mean_resampledBetavalSust(:,rn))));
    grandMn_T(rn) = squeeze(mean(mean_resampledBetavalTrans(:,rn),1,'omitnan'));
    grandSE_T(rn) = squeeze(std(mean_resampledBetavalTrans(:,rn),[],1,'omitnan'))./sqrt(sum(~isnan(mean_resampledBetavalTrans(:,rn))));
end

fH = figure; set(gcf,'Position',[441   195   643   602]);
b = bar([1:length(roisToPlot)],[grandMn_S;grandMn_T],'FaceColor','flat'); hold all;
b(1).CData = roiColors;
b(2).CData = white;
eb1 = errorbar([1:length(roisToPlot)]-0.15,grandMn_S, grandSE_S, grandSE_S,'k','linewidth',0.5);
eb1.LineStyle = 'none';
eb1.CapSize = 0;
eb2 = errorbar([1:length(roisToPlot)]+0.15,grandMn_T, grandSE_T, grandSE_T,'k','linewidth',0.5);
eb2.LineStyle = 'none';
eb2.CapSize = 0;

title('Beta values for CST Sustained & Transient channels');
set(gca,'XTickLabel',string(roisToPlot)', 'XTickLabelRotation',45)
set(gca,'FontSize',12, 'TickDir','out')
box off;
ylabel('Beta weight (a.u.)')
ylim([0 1.8])
legend('Sustained','Transient'); legend box off

if saveFigs
    saveFigDir = fullfile(simseqRootPath,'results','group');
    subDir     = 'fig8';
    fName = sprintf('SummaryGroupAverageSEM_BetaCSTSustainedTransient_resampled');
    thisSaveFigDir = fullfile(saveFigDir,subDir);
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end