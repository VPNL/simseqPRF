function fH = plotDeltaLMMregressionSlopes_v_pRFSize(deltasGroup,deltasSubj, ...
                        resampled_medianPRFSz, deltaIdx, cmapROIs,saveFigs,saveFigDir)

fH = figure; set(gcf,'Position',[1 1 800 941]); hold all;

for idx = 1:size(deltaIdx,1)
    prfSizeDiff = resampled_medianPRFSz(:,deltaIdx(idx,1))-resampled_medianPRFSz(:,deltaIdx(idx,2));
    currMnSz     = mean(prfSizeDiff,1,'omitnan');
    currSEMSz    = std(prfSizeDiff,[],1,'omitnan')./sqrt(sum(~isnan(prfSizeDiff)));
    grandMean    = mean(deltasGroup(:,idx),1); % same as mean(mean(deltasSubj{idx},2)));
    grandSEM     = std(mean(deltasSubj{idx},2,'omitnan'),[],'omitnan')./sqrt(length(mean(deltasSubj{idx},2,'omitnan')));
    scatter(currMnSz,grandMean,100,cmapROIs(idx+1,:),'o','filled','MarkerEdgeColor',[0 0 0]); hold on;
    plot(currMnSz+[-currSEMSz,currSEMSz],[grandMean,grandMean],'color',cmapROIs(idx+1,:));
    plot([currMnSz,currMnSz],grandMean+[-grandSEM,grandSEM],'color',cmapROIs(idx+1,:));
end

title('Group Average +/- SEM: prf size vs Delta suppression slope','FontSize',13)
xlim([0 5]); ylim([-0.25 0.05]); box off; axis square
plot([0 10],[0 0],'k:','linewidth',2)
xlabel('Delta Median pRF size (deg)');
ylabel('Delta suppression (difference in fitted regression slope)')

if saveFigs
    fName = sprintf('SupplFig2b_DeltaLMMfitSlopeVsDeltaMedianPRFSize');
    subDir = 'SupplFig2';
    if ~exist('saveFigDir','dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir,subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose')
end