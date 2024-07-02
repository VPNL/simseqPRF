function fH = plotLMMRegressionSlopes_v_pRFCSTExp(lmmResults,resampled_medianCSTExp,roisToPlot, cmapROIs, saveFigs,saveFigDir)

fH = figure; set(gcf,'Position',[889,294,963,604]);  hold all;

for idx = 1:length(roisToPlot)
    % Compute group average and SEM - CST exp
    currGroupMeanOfMedianExp = mean(resampled_medianCSTExp(:,idx),1,'omitnan');
    currGroupSEMOfMedianExp = std(resampled_medianCSTExp(:,idx),[],1,'omitnan')./sqrt(sum(~isnan(resampled_medianCSTExp(:,idx))));
    % Compute group average and SEM - LMM Slopes
    grandMeanSlope = mean(lmmResults.fixedSlopes{idx}); % same as mean(mean(subjSlopes{idx},2)));
    grandSEMSlope  = std(mean(lmmResults.subjSlopes{idx},2,'omitnan'))./sqrt(length(lmmResults.subjSlopes{idx}(:,1)));
       
    scatter(currGroupMeanOfMedianExp,grandMeanSlope,100,cmapROIs(idx,:),'o', 'filled', 'MarkerEdgeColor',[0 0 0]);
    plot(currGroupMeanOfMedianExp+[-currGroupSEMOfMedianExp,currGroupSEMOfMedianExp],[grandMeanSlope,grandMeanSlope],'k-'); %'color',cmapROIs(idx,:));
    plot([currGroupMeanOfMedianExp,currGroupMeanOfMedianExp],grandMeanSlope+[-grandSEMSlope,grandSEMSlope],'k-');%'color',cmapROIs(idx,:));
end
xlim([0 1]); ylim([0 1]); box off; axis square;
set(gca,'XTick',[0.1:0.1:1],'YTick',[0.1:0.1:1])
xlabel('Median pRF CST exponent (a.u.)');
ylabel('Suppression level (fitted regression slope)')


%% Compute Pearson's correlation (rho)
mnCondSlopes = NaN(length(roisToPlot),size(lmmResults.subjSlopes{1},1));
mnCSTExp     = mnCondSlopes;
for idx = 1:length(roisToPlot)
    mnCondSlopes(idx,~isnan(resampled_medianCSTExp(:,idx))) = mean(lmmResults.subjSlopes{idx},2);
    mnCSTExp(idx,:) = resampled_medianCSTExp(:,idx);
end
[rho, rhoPval, rhoLower, rhoUpper] = corrcoef(mnCondSlopes,mnCSTExp, 'Rows','complete');
fprintf('Pearson''s correlation group average LMM slopes and pRF CST exponent:\n')
fprintf('rho: %1.2f\tp-val: %1.4f\t95%% CI: [%1.2f, %1.2f]\n',...
    rho(1,2), rhoPval(1,2),rhoLower(1,2),rhoUpper(1,2))
title({'Group Average (+/- SEM): median prf CST exponent vs suppression slope',...
        sprintf('rho: %1.2f\tp-val: %1.4f\t95%% CI: [%1.2f, %1.2f]\n',...
        rho(1,2), rhoPval(1,2),rhoLower(1,2),rhoUpper(1,2))},'FontSize',13);

if saveFigs
    saveFigDir = fullfile(simseqRootPath,'results','group');
    subDir     = 'fig8';
    fName = sprintf('SummaryGroupAverageSEM_LMMfitSlopeVsPRFExp');
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end
