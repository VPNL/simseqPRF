function fH = plotLMMregressionSlopes_v_pRFSize_stRetParams(allLMMResults, ...
    params, modelToPlot, roisToPlot, cmapROIs, saveFigs, saveFigDir)

ipsIdx    = cellfind(regexp(string(roisToPlot),'IPS0/IPS1','once'));
roisNoIPS = setdiff([1:length(roisToPlot)],ipsIdx);

fH = figure; set(gcf,'Position',[544   512   687   465]); hold all;
for idx = roisNoIPS
    
    groupMnSzFix     = mean(params{1}.median_resampledPRFSz(:,idx),1,'omitnan');
    groupSEMSzFix    = std(params{1}.median_resampledPRFSz(:,idx),[],1,'omitnan')./sqrt(sum(~isnan(params{1}.median_resampledPRFSz(:,idx))));
    
    groupMnSlopeFix   = mean(mean(allLMMResults{1}.subjSlopes{idx},2,'omitnan')); % same as mean(mean(subjSlopes{idx},2)));
    groupSEMSlopeFix  = std(mean(allLMMResults{1}.subjSlopes{idx},2,'omitnan'))./sqrt(length(mean(allLMMResults{1}.subjSlopes{idx},2,'omitnan')));
    
    scatter(groupMnSzFix,groupMnSlopeFix,100,[1 1 1],'o','filled','MarkerEdgeColor',cmapROIs(idx,:), 'LineWidth',2); hold on;
    plot(groupMnSzFix+[-groupSEMSzFix,groupSEMSzFix],[groupMnSlopeFix,groupMnSlopeFix],'color',[0 0 0]);
    plot([groupMnSzFix,groupMnSzFix],groupMnSlopeFix+[-groupSEMSlopeFix,groupSEMSlopeFix],'color',[0 0 0]);
    
    if strcmp(modelToPlot,'CST_ST')
        groupSzST    = params{2}.median_resampledPRFSz(:,idx);
        groupSlopeST = mean(allLMMResults{2}.subjSlopes{idx},2,'omitnan');
        legendNames  = {'CST_{fix}','CST_{opt}'};
    elseif strcmp(modelToPlot,'DN_ST')
        groupSzST    = params{3}.median_resampledPRFSz(:,idx);
        groupSlopeST = mean(allLMMResults{3}.subjSlopes{idx},2,'omitnan');
        legendNames  = {'CST_{fix}','DN-ST'};
    end
    
    % Compute group mean and SEM
    groupMnSzST     = mean(groupSzST,1,'omitnan');
    groupSEMSzST    = std(groupSzST,[],1,'omitnan')./sqrt(sum(~isnan(groupSzST)));
    
    groupMnSlopeST  = mean(groupSlopeST); % same as mean(mean(subjSlopes{idx},2)));
    groupSEMSlopeST = std(groupSlopeST,[],1,'omitnan')./sqrt(size(groupSlopeST,1));
    
    % plot it!
    scatter(groupMnSzST,groupMnSlopeST,100,cmapROIs(idx,:),'o','filled','MarkerEdgeColor',cmapROIs(idx,:)); hold on;
    plot(groupMnSzST+[-groupSEMSzST,groupSEMSzST],[groupMnSlopeST,groupMnSlopeST],'color',cmapROIs(idx,:));
    plot([groupMnSzST,groupMnSzST],groupMnSlopeST+[-groupSEMSlopeST,groupSEMSlopeST],'color',cmapROIs(idx,:));

end

% Set legends and axes
xlim([0.2 10.5]); ylim([0 1]); box off; axis square
xlabel('Median pRF size (deg)');
ylabel('Suppression (fitted regression slope)')
l = gca;
legend(l.Children([length(l.Children):-3:length(l.Children)-6]),...
    legendNames, 'FontSize',9, 'Location','SouthWest');
legend boxoff

% Compute Pearson's correlation (rho)
% between group average condition slopes and pRF size
mnCondSlopesST  = NaN(length(roisNoIPS),size(allLMMResults{1}.subjSlopes{1},1));
mnCondSlopesFix = mnCondSlopesST;
mnSzFix         = mnCondSlopesST;
mnSzST          = mnCondSlopesST;

for idx = roisNoIPS
    mnCondSlopesFix(idx,~isnan(params{1}.median_resampledPRFSz(:,idx))) = mean(allLMMResults{1}.subjSlopes{idx},2);
    mnSzFix(idx,:) = params{1}.median_resampledPRFSz(:,idx);
    
    mnCondSlopesST(idx,~isnan(params{2}.median_resampledPRFSz(:,idx))) = mean(allLMMResults{2}.subjSlopes{idx},2);
    mnSzST(idx,:) = params{2}.median_resampledPRFSz(:,idx);
end

% Print stats
fprintf('Pearson''s correlation group average LMM slopes and pRF size:\n')
[rho_ST, rhoPval_ST, rhoLower_ST, rhoUpper_ST] = corrcoef(mnCondSlopesST,mnSzST, 'Rows','complete');
if strcmp(modelToPlot,'CST_ST') 
    [rho_Fix, rhoPval_Fix, rhoLower_Fix, rhoUpper_Fix] = corrcoef(mnCondSlopesFix,mnSzFix, 'Rows','complete');
    fprintf('CSTfix model: rho: %1.2f\tp-val: %1.5f\t95%% CI: [%1.2f, %1.2f]\n',...
        rho_Fix(1,2), rhoPval_Fix(1,2),rhoLower_Fix(1,2),rhoUpper_Fix(1,2)),...
    fprintf('CSTopt model: rho: %1.2f\tp-val: %1.5f\t95%% CI: [%1.2f, %1.2f]\n',...
        rho_ST(1,2), rhoPval_ST(1,2),rhoLower_ST(1,2),rhoUpper_ST(1,2))

    title({'Group Average (+/- SEM): median pRF size vs suppression slope',...
        sprintf('CST_{opt} model: rho=%1.2f CI=[%1.2f,%1.2f] (p=%1.5f)',...
        rho_ST(1,2),rhoLower_ST(1,2),rhoUpper_ST(1,2),rhoPval_ST(1,2)),...
        sprintf('CST_{fix} model: rho=%1.2f CI=[%1.2f,%1.2f] (p=%1.5f)',...
        rho_Fix(1,2),rhoLower_Fix(1,2),rhoUpper_Fix(1,2),rhoPval_Fix(1,2))},...
        'FontSize',13);
else
    title({'Group Average (+/- SEM): median pRF size vs suppression slope',...
    sprintf('%s model: rho=%1.2f CI=[%1.2f,%1.2f] (p=%1.5f)',modelToPlot,...
    rho_ST(1,2),rhoLower_ST(1,2),rhoUpper_ST(1,2),rhoPval_ST(1,2))},...
    'FontSize',13);
end

if saveFigs
    if strcmp(modelToPlot, 'DN_ST')
        fName = sprintf('FigS9F_GroupAverageSlopes_v_pRFtau_%s',modelToPlot);
    else
        fName = sprintf('FigS9B_GroupAverageSlopes_v_pRFsz_%s',modelToPlot);
    end
    % Set paths
    subDir = 'SupplFig9';
    if ~exist('saveFigDir','dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir,subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
    print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose')
end
