function fH = plotLMMregressionSlopes_v_pRFtau_stRetParams(allLMMResults,...
    params, modelToPlot, roisToPlot, cmapROIs, saveFigs, saveFigDir)

ipsIdx    = cellfind(regexp(string(roisToPlot),'IPS0/IPS1','once'));
roisNoIPS = setdiff([1:length(roisToPlot)],ipsIdx);

fH = figure; set(gcf,'Position',[329   450   768   527]); hold all;
for idx = roisNoIPS
    if strcmp(modelToPlot,'CST_ST')
        legendNames = {'CST_{fix}','CST_{opt}'};
        
        % Fixed tau
        groupMnTauFix     = 49.3; % ms (see Stigliani et al. 2017)
        groupSEMTauFix    = 0;
        % Compute group mean and SEM for suppression CST opt slopes
        groupMnSlopeFix   = mean(mean(allLMMResults{1}.subjSlopes{idx},2,'omitnan')); % same as mean(mean(subjSlopes{idx},2)));
        groupSEMSlopeFix  = std(mean(allLMMResults{1}.subjSlopes{idx},2,'omitnan'))./sqrt(size(allLMMResults{1}.subjSlopes{idx},2));
        % Compute group mean and SEM for suppression CST fix slopes
        groupMnSlopeST    = mean(mean(allLMMResults{2}.subjSlopes{idx},2,'omitnan'));
        groupSEMSlopeST   = std(mean(allLMMResults{2}.subjSlopes{idx},2,'omitnan'))./sqrt(size(allLMMResults{2}.subjSlopes{idx},2));
        % Compute group mean and SEM for tau pRF parameter
        groupMnTau1ST     = mean(params{2}.median_resampledCST_tau(:,idx),1,'omitnan');
        groupSEMTau1ST    = std(params{2}.median_resampledCST_tau(:,idx),[],1,'omitnan')./sqrt(sum(~isnan(params{2}.median_resampledCST_tau(:,idx))));
        
        % Plot it!
        scatter(groupMnTauFix,groupMnSlopeFix,100,[1 1 1],'o','filled','MarkerEdgeColor',cmapROIs(idx,:)); hold on;
        plot(groupMnTauFix+[-groupSEMTauFix,groupSEMTauFix],[groupMnSlopeFix,groupMnSlopeFix],'color',cmapROIs(idx,:));
        plot([groupMnTauFix,groupMnTauFix],groupMnSlopeFix+[-groupSEMSlopeFix,groupSEMSlopeFix],'color',cmapROIs(idx,:));
        
        scatter(groupMnTau1ST,groupMnSlopeST,100,cmapROIs(idx,:),'o','filled','MarkerEdgeColor',cmapROIs(idx,:)); hold on;
        plot(groupMnTau1ST+[-groupSEMTau1ST,groupSEMTau1ST],[groupMnSlopeST,groupMnSlopeST],'color',cmapROIs(idx,:));
        plot([groupMnTau1ST,groupMnTau1ST],groupMnSlopeST+[-groupSEMSlopeST,groupSEMSlopeST],'color',cmapROIs(idx,:));
        
    elseif strcmp(modelToPlot,'DN_ST')
        legendNames = {'DN-ST tau1','DN-ST tau2'};
        % Compute group mean and SEM for suppression CST opt slopes
        groupMnSlopeST   = mean(mean(allLMMResults{3}.subjSlopes{idx},2,'omitnan'));
        groupSEMSlopeST  = std(mean(allLMMResults{3}.subjSlopes{idx},2,'omitnan'))./sqrt(size(allLMMResults{3}.subjSlopes{idx},2));
        % Compute group mean and SEM for tau pRF parameter
        groupMnTau1ST    = mean(1000*params{3}.median_resampledDNST_tau1(:,idx),1,'omitnan');
        groupSEMTau1ST   = std(1000*params{3}.median_resampledDNST_tau1(:,idx),[],1,'omitnan')./sqrt(sum(~isnan(params{3}.median_resampledDNST_tau1(:,idx))));
        groupMnTau2ST    = mean(1000*params{3}.median_resampledDNST_tau2(:,idx),1,'omitnan');
        groupSEMTau2ST   = std(1000*params{3}.median_resampledDNST_tau2(:,idx),1,'omitnan')./sqrt(sum(~isnan(params{3}.median_resampledDNST_tau2(:,idx))));

        % Plot it!
        scatter(groupMnTau1ST,groupMnSlopeST,100,cmapROIs(idx,:),'o', 'filled', 'MarkerEdgeColor',cmapROIs(idx,:));
        plot(groupMnTau1ST+[-groupSEMTau1ST,groupSEMTau1ST],[groupMnSlopeST,groupMnSlopeST],'color',cmapROIs(idx,:));
        plot([groupMnTau1ST,groupMnTau1ST],groupMnSlopeST+[-groupSEMSlopeST,groupSEMSlopeST],'color',cmapROIs(idx,:));
        
        scatter(groupMnTau2ST,groupMnSlopeST,100,cmapROIs(idx,:),'o', 'MarkerEdgeColor',cmapROIs(idx,:));
        plot(groupMnTau2ST+[-groupSEMTau2ST,groupSEMTau2ST],[groupMnSlopeST,groupMnSlopeST],'color',cmapROIs(idx,:));
        plot([groupMnTau2ST,groupMnTau2ST],groupMnSlopeST+[-groupSEMSlopeST,groupSEMSlopeST],'color',cmapROIs(idx,:));
    end
end

% Set legends and axes
xlim([0 300]); ylim([0 1]); box off; axis square
set(gca,'XTick',[0:100:300],'YTick',[0.2:0.2:1])
ylabel('Suppression (fitted regression slope)')
xlabel('Median pRF tau (ms)');
l = gca;
legend(l.Children([length(l.Children):-9:length(l.Children)-9]),...
    legendNames, 'FontSize',9, 'Location','SouthWest');
legend boxoff

% Compute Pearson's correlation (rho)
% between group average condition slopes and pRF tau
mnCondSlopesST  = NaN(length(roisNoIPS),size(allLMMResults{1}.subjSlopes{1},1));
mnTau1ST        = mnCondSlopesST;
mnTau2ST        = mnCondSlopesST;

for idx = roisNoIPS
    if strcmp(modelToPlot,'CST_ST')
        mnCondSlopesST(idx,~isnan(params{2}.median_resampledCST_tau(:,idx))) = mean(allLMMResults{2}.subjSlopes{idx},2);
        mnTau1ST(idx,:) = params{2}.median_resampledCST_tau(:,idx);

    elseif strcmp(modelToPlot,'DN_ST')
        mnCondSlopesST(idx,~isnan(params{3}.median_resampledDNST_tau1(:,idx))) = mean(allLMMResults{3}.subjSlopes{idx},2);
        mnTau1ST(idx,:) = params{3}.median_resampledDNST_tau1(:,idx);
        mnTau2ST(idx,:) = params{3}.median_resampledDNST_tau2(:,idx);        
    end
end

% Print stats
if strcmp(modelToPlot,'CST_ST')
    [rho_ST, rhoPval_ST, rhoLower_ST, rhoUpper_ST] = corrcoef(mnCondSlopesST,mnTau1ST, 'Rows','complete');
    fprintf('Pearson''s correlation group average LMM slopes and pRF tau:\n')
    fprintf('CSTopt model: rho: %1.2f\tp-val: %1.4f\t95%% CI: [%1.2f, %1.2f]\n',...
        rho_ST(1,2), rhoPval_ST(1,2),rhoLower_ST(1,2),rhoUpper_ST(1,2))
    title({'Group Average (+/- SEM): median pRF tau vs suppression slope',...
        sprintf('CST_{opt} model: rho=%1.2f CI=[%1.2f,%1.2f] (p=%1.3f)',...
        rho_ST(1,2),rhoLower_ST(1,2),rhoUpper_ST(1,2),rhoPval_ST(1,2))},...
        'FontSize',13);
elseif strcmp(modelToPlot,'DN_ST')
    [rho1, rhoPval1, rhoLower1, rhoUpper1] = corrcoef(mnCondSlopesST,mnTau1ST, 'Rows','complete');
    fprintf('Pearson''s correlation group average LMM slopes and pRF DN-ST tau1:\n')
    fprintf('rho: %1.2f\tp-val: %1.4f\t95%% CI: [%1.2f, %1.2f]]\n',...
        rho1(1,2), rhoPval1(1,2),rhoLower1(1,2),rhoUpper1(1,2))
    [rho2, rhoPval2, rhoLower2, rhoUpper2] = corrcoef(mnCondSlopesST,mnTau2ST, 'Rows','complete');
    fprintf('\nPearson''s correlation group average LMM slopes and pRF DN-ST tau2:\n')
    fprintf('rho: %1.2f\tp-val: %1.4f\t95%% CI: [%1.2f, %1.2f]\n',...
        rho2(1,2), rhoPval2(1,2),rhoLower2(1,2),rhoUpper2(1,2))
    
    title({'Group Average (+/- SEM): median pRF DN-ST_tau1 vs suppression slope',...
        sprintf('rho1=%1.2f CI=[%1.2f,%1.2f] (p=%1.3f)',rho1(1,2),rhoLower1(1,2),rhoUpper1(1,2),rhoPval1(1,2))},...
        {'Group Average (+/- SEM): median pRF DN-ST_tau2 vs suppression slope',...
        sprintf('rho2=%1.2f CI=[%1.2f,%1.2f] (p=%1.3f)',rho2(1,2),rhoLower2(1,2),rhoUpper2(1,2),rhoPval2(1,2))},...
        'FontSize',13);

end


if saveFigs
    if strcmp(modelToPlot, 'DN_ST')
        fName = sprintf('FigS9H_GroupAverageSlopes_v_pRFtau_%s',modelToPlot);
    else
        fName = sprintf('FigS9D_GroupAverageSlopes_v_pRFtau_%s',modelToPlot);
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
