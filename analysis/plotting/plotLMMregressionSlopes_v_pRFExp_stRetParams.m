function fH = plotLMMregressionSlopes_v_pRFExp_stRetParams(allLMMResults, ...
    params, modelToPlot, roisToPlot, cmapROIs, saveFigs, saveFigDir)

ipsIdx    = cellfind(regexp(string(roisToPlot),'IPS0/IPS1','once'));
roisNoIPS = setdiff([1:length(roisToPlot)],ipsIdx);

fH = figure; set(gcf,'Position',[544   512   687   465]); hold all;
for idx = roisNoIPS
    if strcmp(modelToPlot,'CST_ST')
        groupMnExpFix     = mean(params{1}.median_resampledCST_exp(:,idx),1,'omitnan');
        groupSEMSzFix     = std(params{1}.median_resampledCST_exp(:,idx),[],1,'omitnan')./sqrt(sum(~isnan(params{1}.median_resampledCST_exp(:,idx))));
        
        groupMnSlopeFix   = mean(mean(allLMMResults{1}.subjSlopes{idx},2,'omitnan')); 
        groupSEMSlopeFix  = std(mean(allLMMResults{1}.subjSlopes{idx},2,'omitnan'))./sqrt(length(mean(allLMMResults{1}.subjSlopes{idx},2,'omitnan')));
        
        scatter(groupMnExpFix,groupMnSlopeFix,100,[1 1 1],'o','filled','MarkerEdgeColor',cmapROIs(idx,:),'LineWidth',2); hold on;
        plot(groupMnExpFix+[-groupSEMSzFix,groupSEMSzFix],[groupMnSlopeFix,groupMnSlopeFix],'color',[0 0 0]);
        plot([groupMnExpFix,groupMnExpFix],groupMnSlopeFix+[-groupSEMSlopeFix,groupSEMSlopeFix],'color',[0 0 0]);
        
        groupExpST   = params{2}.median_resampledCST_exp(:,idx);
        groupSlopeST = mean(allLMMResults{2}.subjSlopes{idx},2,'omitnan');
        legendNames  = {'CST_{fix}','CST_{opt}'};
    elseif strcmp(modelToPlot,'DN_ST')
        groupExpST   = params{3}.median_resampledDNST_n(:,idx);
        groupSlopeST = mean(allLMMResults{3}.subjSlopes{idx},2,'omitnan');
        legendNames  = {'DN_ST'};
    end
    
    % Compute group mean and SEM
    groupMnSzST    = mean(groupExpST,1,'omitnan');
    groupSEMSzST   = std(groupExpST,[],1,'omitnan')./sqrt(sum(~isnan(groupExpST)));
    
    groupMnSlopeST  = mean(groupSlopeST); % same as mean(mean(subjSlopes{idx},2)));
    groupSEMSlopeST = std(groupSlopeST,[],1,'omitnan')./sqrt(size(groupSlopeST,1));
    
    % plot it!
    scatter(groupMnSzST,groupMnSlopeST,100,cmapROIs(idx,:),'o','filled','MarkerEdgeColor',cmapROIs(idx,:)); hold on;
    plot(groupMnSzST+[-groupSEMSzST,groupSEMSzST],[groupMnSlopeST,groupMnSlopeST],'color',cmapROIs(idx,:));
    plot([groupMnSzST,groupMnSzST],groupMnSlopeST+[-groupSEMSlopeST,groupSEMSlopeST],'color',cmapROIs(idx,:));

end

% Set legends and axes
if strcmp(modelToPlot,'CST_ST')
    xlim([0 1]); ylim([0 1]); box off; axis square
    set(gca,'XTick',[0.2:0.2:1],'YTick',[0.2:0.2:1])
elseif strcmp(modelToPlot,'DN_ST')
    xlim([0.9 2.1]); ylim([0 1]); box off; axis square
    set(gca,'XTick',[1.2:0.2:2],'YTick',[0.2:0.2:1])
end
xlabel('Median pRF exponent');
ylabel('Suppression (fitted regression slope)')
l = gca;
legend(l.Children([length(l.Children):-3:length(l.Children)-6]),...
    legendNames, 'FontSize',9, 'Location','SouthWest');
legend boxoff

% Compute Pearson's correlation (rho)
% between group average condition slopes and pRF size
mnCondSlopesST  = NaN(length(roisNoIPS),size(allLMMResults{1}.subjSlopes{1},1));
mnSzFix         = mnCondSlopesST;
mnCondSlopesFix = mnCondSlopesST;
mnExpST          = mnCondSlopesST;

for idx = roisNoIPS
    if strcmp(modelToPlot,'CST_ST')
        mnCondSlopesST(idx,~isnan(params{2}.median_resampledCST_exp(:,idx))) = mean(allLMMResults{2}.subjSlopes{idx},2);
        mnExpST(idx,:) = params{2}.median_resampledCST_exp(:,idx);
        
        mnCondSlopesFix(idx,~isnan(params{1}.median_resampledCST_exp(:,idx))) = mean(allLMMResults{1}.subjSlopes{idx},2);
        mnSzFix(idx,:) = params{1}.median_resampledCST_exp(:,idx);
        
    elseif strcmp(modelToPlot,'DN_ST')
        mnCondSlopesST(idx,~isnan(params{3}.median_resampledDNST_n(:,idx))) = mean(allLMMResults{3}.subjSlopes{idx},2);
        mnExpST(idx,:) = params{3}.median_resampledDNST_n(:,idx);
        
    end
end

% Print stats
[rho_ST, rhoPval_ST, rhoLower_ST, rhoUpper_ST] = corrcoef(mnCondSlopesST,mnExpST, 'Rows','complete');
fprintf('Pearson''s correlation group average LMM slopes and pRF exponent:\n')
fprintf('%s model: rho: %1.2f\tp-val: %1.4f\t95%% CI: [%1.2f, %1.2f]\n',...
    rho_ST(1,2), rhoPval_ST(1,2),rhoLower_ST(1,2),rhoUpper_ST(1,2))
if strcmp(modelToPlot,'CST_ST')
    [rho_Fix, rhoPval_Fix, rhoLower_Fix, rhoUpper_Fix] = corrcoef(mnCondSlopesFix,mnSzFix, 'Rows','complete');
    fprintf('Fix model: rho: %1.2f\tp-val: %1.4f\t95%% CI: [%1.2f, %1.2f]\n',...
        rho_Fix(1,2), rhoPval_Fix(1,2),rhoLower_Fix(1,2),rhoUpper_Fix(1,2))
    title({'Group Average (+/- SEM): median pRF exponent vs suppression slope',...
        sprintf('CST_{opt} model: rho=%1.2f CI=[%1.2f,%1.2f] (p=%1.4f)',rho_ST(1,2),rhoLower_ST(1,2),rhoUpper_ST(1,2),rhoPval_ST(1,2)),...
        sprintf('CST_{fix} model: rho=%1.2f CI=[%1.2f,%1.2f] (p=%1.4f)',rho_Fix(1,2),rhoLower_Fix(1,2),rhoUpper_Fix(1,2),rhoPval_Fix(1,2))},...
        'FontSize',13);
else
    title({'Group Average (+/- SEM): median pRF exponent vs suppression slope',...
    sprintf('%s model: rho=%1.2f CI=[%1.2f,%1.2f] (p=%1.4f)',modelToPlot,...
    rho_ST(1,2),rhoLower_ST(1,2),rhoUpper_ST(1,2),rhoPval_ST(1,2))},...
    'FontSize',13);
end


if saveFigs
    if strcmp(modelToPlot, 'DN_ST')
        fName = sprintf('FigS9G_GroupAverageSlopes_v_pRFexp_%s',modelToPlot);
    else
        fName = sprintf('FigS9C_GroupAverageSlopes_v_pRFtau_%s',modelToPlot);
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
