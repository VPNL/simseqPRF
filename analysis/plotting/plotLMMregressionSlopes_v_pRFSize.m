function fH = plotLMMregressionSlopes_v_pRFSize(lmmResults, resampled_medianPRFSz, roisToPlot, cmapROIs,saveFigs,saveFigDir)

fH = figure; set(gcf,'Position',[889,294,963,604]); hold all;

for idx = 1:length(roisToPlot)
    currMnSz     = mean(resampled_medianPRFSz(:,idx),1,'omitnan');
    currSEMSz    = std(resampled_medianPRFSz(:,idx),[],1,'omitnan')./sqrt(sum(~isnan(resampled_medianPRFSz(:,idx))));
    grandMean = mean(lmmResults.fixedSlopes{idx}); % same as mean(mean(subjSlopes{idx},2)));
    grandSEM  = std(mean(lmmResults.subjSlopes{idx},2,'omitnan'))./sqrt(length(mean(lmmResults.subjSlopes{idx},2,'omitnan')));
    scatter(currMnSz,grandMean,100,cmapROIs(idx,:),'o','filled','MarkerEdgeColor',[0 0 0]); hold on;
    plot(currMnSz+[-currSEMSz,currSEMSz],[grandMean,grandMean],'-','color',cmapROIs(idx,:));
    plot([currMnSz,currMnSz],grandMean+[-grandSEM,grandSEM],'-','color',cmapROIs(idx,:));
end
xlim([0.2 10]); ylim([0 1]); box off; axis square
xlabel('Median pRF size (deg)');
ylabel('Suppression (fitted regression slope)')

%% Compute Pearson's correlation (rho) 
% between group average condition slopes and pRF size
mnCondSlopes = NaN(9,10);
mnSz = mnCondSlopes;
for idx = 1:length(roisToPlot)
    mnCondSlopes(idx,~isnan(resampled_medianPRFSz(:,idx))) = mean(lmmResults.subjSlopes{idx},2);
    mnSz(idx,:) = resampled_medianPRFSz(:,idx);
end


%% Print stats
[rho, rhoPval, rhoLower, rhoUpper] = corrcoef(mnCondSlopes,mnSz, 'Rows','complete');
fprintf('Pearson''s correlation group average LMM slopes and pRF size:\n')
fprintf('rho: %1.3f\tp-val: %1.4f\t95%% CI: [%1.2f, %1.2f]\n',...
    rho(1,2), rhoPval(1,2),rhoLower(1,2),rhoUpper(1,2))
    
title({'Group Average (+/- SEM): median pRF size vs suppression slope',...
    sprintf('CST model: Pearson''s r=%1.3f CI=[%1.2f,%1.2f] (p=%1.5f)',...
    rho(1,2),rhoPval(1,2),rhoUpper(1,2),rhoPval(1,2))},...
    'FontSize',13);

if saveFigs
    saveFigDir = fullfile(simseqRootPath,'results','group');
    subDir     = 'fig8';
    fName = sprintf('SummaryGroup_LMMfitSlopeVsMedianPRFSize');
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
    print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end



