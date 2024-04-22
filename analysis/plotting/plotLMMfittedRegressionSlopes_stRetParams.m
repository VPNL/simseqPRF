function fH = plotLMMfittedRegressionSlopes_stRetParams(allDS,allLMMResults,...
                allLMMResults_model,LMMOrder,...
                roisToPlot,cmapModels, saveFigs, saveFigDir)
% Suppl Fig 8: Model performance of alternative spatiotemporal pRF models,
% including CST with fixed temporal parameters (CST_fix), 
%           CST with optimized spatial and temporal parameters (CST opt),
%           Delayed normalization model with optimized parameters (DN-ST).

%% Panel A: predicted regression slopes (alike Figure 7)

% Set params
conditionOrderSimSeq = 1:4;
subplotOrder         = [2,1,4,3];
conditionNamesSimSeq = {'Small & Short (0.2 s)','Small & Long (1 s)',...
                        'Big & Short (0.2 s)','Big & Long (1 s)'};
ipsIdx    = cellfind(regexp(string(roisToPlot),'IPS0/IPS1','once'));
roisNoIPS = setdiff([1:length(roisToPlot)],ipsIdx);

if ~exist('cmapModels', 'var') || isempty(cmapModels)
    cmapModels = getColormapPRFModels(3);
end

% Plot it!
fH = figure; clf; set(gcf,'Position',[669 30 1342 905]);
x = [1:length(roisToPlot)]+[-0.4,0.4]';

for mm = 1:length(LMMOrder)

    clear dataMat seMat modelMat modelseMat
    for ii = roisNoIPS
        sj = unique(allDS{mm}.ds.Subject(allDS{mm}.ds.ROI==roisToPlot(ii)))';
        dataMat(:,ii) = mean(allLMMResults{mm}.subjSlopes{ii},1,'omitnan');
        seMat(:,ii)   = std(allLMMResults{mm}.subjSlopes{ii},[],1,'omitnan')/sqrt(length(sj));
        
        sj = unique(allDS{mm}.ds.Subject(allDS{mm}.ds.ROI==roisToPlot(ii)))';
        modelMat(:,ii) = mean(allLMMResults_model{mm}.subjSlopes{ii},1,'omitnan');
        modelseMat(:,ii) = std(allLMMResults_model{mm}.subjSlopes{ii},[],1,'omitnan')./sqrt(length(sj));
    end
    
    
    for c = conditionOrderSimSeq
        subplot(1,length(conditionOrderSimSeq),subplotOrder(c)); hold all;
        cData = repmat(dataMat(c,:),2,1);
        seData = repmat(seMat(c,:),2,1);
        
        if c==subplotOrder(c), ylabel('Suppression (fitted regression slope)'); end
        
        for idx = roisNoIPS
            
            % Plot data as bars
            errorbar3(x(:,idx)',cData(:,idx)',seData(:,idx)',1,cmapModels(mm,:));  % fourth input is error on y-dim)
            ax = gca;
            ax.Children(1).FaceAlpha = 0.3;
            curMnSlope   = modelMat(c,idx);
            SE           = modelseMat(c,idx);
            
            scatter(idx,curMnSlope,100,cmapModels(mm,:),'filled','linewidth',0.15); hold on; %'MarkerEdgeColor',[0,0,0]
            plot([idx, idx], curMnSlope + [-SE, SE],'color',cmapModels(mm,:), 'linewidth',2);
            
        end
        title(sprintf('%s',conditionNamesSimSeq{c}))
        xlim([0.2 9.7]); 
        ylim([0 1.21]); box off
        set(gca,'XTick',1:length(roisToPlot), 'XTickLabel',string(roisToPlot),'XTickLabelRotation',45)
    end
end


% Add legend and title
l = gca;
legend(l.Children([length(l.Children):-30:1]),LMMOrder, 'FontSize',9, 'Location','SouthWest'); legend boxoff
sgtitle(sprintf('Group Average Data vs Model: LMM regression slopes (Mean +/-SE)'))

% Save figure if requested
if saveFigs
    fName = sprintf('FigS9A_SummarySlopesSEQvsSIMMODELAmpl%s_LMMfit',cell2mat(LMMOrder));

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

