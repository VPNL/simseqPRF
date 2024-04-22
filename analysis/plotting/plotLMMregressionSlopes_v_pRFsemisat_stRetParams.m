function fH = plotLMMregressionSlopes_v_pRFsemisat_stRetParams(allLMMResults, ...
    params, roisToPlot, cmapROIs, saveFigs, saveFigDir)

ipsIdx    = cellfind(regexp(string(roisToPlot),'IPS0/IPS1','once'));
roisNoIPS = setdiff([1:length(roisToPlot)],ipsIdx);

fH = figure; set(gcf,'Position',[544   512   687   465]); hold all;

for idx = roisNoIPS
    % Compute group average and SEM - DNST semisat
    currGroupMeanOfMedianDNsemisat = mean(params{3}.median_resampledDNST_semisat(:,idx),1,'omitnan');
    currGroupSEMOfMedianDNsemisat = std(params{3}.median_resampledDNST_semisat(:,idx),[],1,'omitnan')./...
                                    sqrt(sum(~isnan(params{3}.median_resampledDNST_semisat(:,idx))));
    
    % Compute group average and SEM - LMM Slopes
    grandMeanSlope = mean(mean(allLMMResults{3}.subjSlopes{idx},2,'omitnan'),'omitnan'); 
    grandSEMSlope  = std(mean(allLMMResults{3}.subjSlopes{idx},2,'omitnan'))./sqrt(size(allLMMResults{3}.subjSlopes{idx},1));
    
    scatter(currGroupMeanOfMedianDNsemisat,grandMeanSlope,100,cmapROIs(idx,:),'o', 'filled', 'MarkerEdgeColor',cmapROIs(idx,:));
    plot(currGroupMeanOfMedianDNsemisat+[-currGroupSEMOfMedianDNsemisat,currGroupSEMOfMedianDNsemisat],[grandMeanSlope,grandMeanSlope],'color',cmapROIs(idx,:));
    plot([currGroupMeanOfMedianDNsemisat,currGroupMeanOfMedianDNsemisat],grandMeanSlope+[-grandSEMSlope,grandSEMSlope],'color',cmapROIs(idx,:));
end

% title('Mean across conditions','FontSize',13)
xlim([0 0.1]); ylim([0 1]); box off; axis square;
set(gca,'XTick',[0:0.01:1],'YTick',[0.1:0.1:1])
xlabel('Median pRF DN-ST semisaturation (a.u.)');
ylabel('Suppression level (fitted regression slope)')
l = gca;
legend(l.Children([(length(l.Children)-1):-3:1]),string(roisToPlot), ...
    'FontSize',9, 'Location','SouthWest'); legend boxoff

% Compute Pearson's correlation (rho)
mnCondSlopes = NaN(length(roisToPlot),size(allLMMResults{3}.subjSlopes{1},1));
mnDNSTsigma     = mnCondSlopes;
for idx = roisNoIPS
    mnCondSlopes(idx,~isnan(params{3}.median_resampledDNST_semisat(:,idx))) = mean(allLMMResults{3}.subjSlopes{idx},2,'omitnan');
    mnDNSTsigma(idx,:) = params{3}.median_resampledDNST_semisat(:,idx);
end
[rho, rhoPval, rhoLower, rhoUpper] = corrcoef(mnCondSlopes,mnDNSTsigma, 'Rows','complete');
fprintf('Pearson''s correlation group average LMM slopes and pRF DN-ST semisat:\n')
fprintf('rho: %1.2f\tp-val: %1.4f\t95%% CI: [%1.2f, %1.2f]\n',...
    rho(1,2), rhoPval(1,2),rhoLower(1,2),rhoUpper(1,2))

title({'Group Average (+/- SEM): median prf DN-ST semisaturation vs suppression slope',...
    sprintf('rho=%1.2f CI=[%1.2f,%1.2f] (p=%1.3f)',rho(1,2),rhoLower(1,2),rhoUpper(1,2),rhoPval(1,2))},...
    'FontSize',13);

if saveFigs
    fName = 'FigS9I_GroupAverageSlopes_v_pRFsemisatconstant_DNST';

    subDir = 'SupplFig9';
    if ~exist('saveFigDir','dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
    print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end