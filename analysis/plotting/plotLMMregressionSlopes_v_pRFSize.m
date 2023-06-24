function fH = plotLMMregressionSlopes_v_pRFSize(lmmResults, resampled_medianPRFSz, roisToPlot, cmapROIs,saveFigs)

conditionNamesSimSeq = {'Small & short (0.2s)','Small & long (1s)',...
                        'Big & Short (0.2s)','Big & Long (1s)'};

fH = figure; set(gcf,'Position',[889,294,963,604]); hold all;

for idx = 1:length(roisToPlot)
     
    currMnSz     = mean(resampled_medianPRFSz(:,idx),1,'omitnan');
    currSEMSz    = std(resampled_medianPRFSz(:,idx),[],1,'omitnan')./sqrt(sum(~isnan(resampled_medianPRFSz(:,idx))));
    grandMean = mean(lmmResults.fixedSlopes{idx}); % same as mean(mean(subjSlopes{idx},2)));
    grandSEM  = std(mean(lmmResults.subjSlopes{idx},2,'omitnan'))./sqrt(length(mean(lmmResults.subjSlopes{idx},2,'omitnan')));
    scatter(currMnSz,grandMean,100,cmapROIs(idx,:),'o','filled','MarkerEdgeColor',[0 0 0]); hold on;
    plot(currMnSz+[-currSEMSz,currSEMSz],[grandMean,grandMean],'k-');
    plot([currMnSz,currMnSz],grandMean+[-grandSEM,grandSEM],'k-');
end
sgtitle('Group Average +/- SEM: prf size vs suppression slope','FontSize',13)
xlim([0.2 10]); ylim([0 1]); box off; axis square
xlabel('Median pRF size (deg)');
ylabel('Suppression (fitted regression slope)')

if saveFigs
    saveFigDir = fullfile(simseqRootPath,'results','group');
    subDir     = 'fig8';
    fName = sprintf('SummaryGroup_LMMfitSlopeVsMedianPRFSize');
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
    print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end

%% Compute Pearson's correlation (rho) 
% between group average condition slopes and pRF size
mnCondSlopes = NaN(9,10);
mnSz = mnCondSlopes;
for idx = 1:9
    mnCondSlopes(idx,~isnan(resampled_medianPRFSz(:,idx))) = mean(lmmResults.subjSlopes{idx},2);
    mnSz(idx,:) = resampled_medianPRFSz(:,idx);
end
[rho, rhoPval, rhoLower, rhoUpper] = corrcoef(mnCondSlopes,mnSz, 'Rows','complete');
fprintf('Pearson''s correlation group average LMM slopes and pRF size:\n')
fprintf('rho: %1.2f\tp-val: %1.4f\t95%% CI: [%1.2f, %1.2f]\n',...
    rho(1,2), rhoPval(1,2),rhoLower(1,2),rhoUpper(1,2))
