function fH = plotLMMfittedRegressionSlopes(ds,lmmResults,LMMOrder, roisToPlot,cmapROIs, useSTRetParams, saveFigs, saveFigDir)
% Function to plot Linear Mixed Model fitted regression slopes (fixed and
% random subject slopes, for each condition and visual area).
% To plot Figure 4B, "LMMResults" struct has fields with size {1,nrROIs}
%   ex: plotLMMfittedRegressionSlopes(LMMResults,{'Data'}, false)
% To plot Figure 7, "LMMResults" struct has fields with size {4,nrROIs}
%   ex: plotLMMfittedRegressionSlopes(LMMResults,{'LSS','CSS','CST','Data'}, false)
%
% Code written by E.R. Kupers (2024) Stanford University
%
%% Check inputs
if size(lmmResults,1) == 1 && size(lmmResults.fixedSlopes,1)==1
    plotModelAmpl = false;
else
    plotModelAmpl = true;
end
if ~exist('saveFigs','var') || isempty(saveFigs)
    saveFigs = false; % Save figures or not
end

if ~exist('useSTRetParams','var') || isempty(useSTRetParams)
    useSTRetParams = false; 
end

% Set params
conditionOrderSimSeq = 1:4;
subplotOrder         = [2,1,4,3];
conditionNamesSimSeq = {'Small & Short (0.2s)','Small & Long (1s)','Big & Short (0.2s)','Big & Long (1s)'};

if plotModelAmpl
    if strcmp(LMMOrder{1},'DoG') % DoG model
         cmap = getColormapPRFModels(2);
    elseif size(lmmResults,1)>1
        cmap = getColormapPRFModels(0);
    end
else
    cmap = cmapROIs;
end

%% Plot it!

fH = figure; clf; set(gcf,'Position',[669 30 1342 905]);

dataMat = reshape(cell2mat(lmmResults(end).fixedSlopes),4,[]);
seMat = reshape(cell2mat(lmmResults(end).fixedSlopes_SE),4,[]);

if useSTRetParams
    clear dataMat seMat
    for ii = 1:length(roisToPlot)
        sj = unique(ds.Subject(ds.ROI==roisToPlot(ii)))';
        sj2 = ismember(sj,[1:7]);
        dataMat(:,ii) = mean(lmmResults(end).subjSlopes{ii}(sj2,:),1,'omitnan');
        seMat(:,ii)   = std(lmmResults(end).subjSlopes{ii}(sj2,:),[],1,'omitnan')/sqrt(length(sj2));
    end
end

x = [1:9]+[-0.4,0.4]';
for mm = 1:length(LMMOrder)
    
    if ~strcmp(LMMOrder{mm},'Data')
        clear modelMat modelseMat
        for ii = 1:length(roisToPlot)
            if useSTRetParams
                sj = unique(ds.Subject(ds.ROI==roisToPlot(ii)))';
                sj2 = ismember(sj,[1:7]);
                modelMat(:,ii) = mean(lmmResults(1).subjSlopes{mm,ii}(sj2,:),1,'omitnan');
                modelseMat(:,ii) = std(lmmResults(1).subjSlopes{mm,ii}(sj2,:),[],1,'omitnan')./sqrt(length(sj2));
            else
                modelMat(:,ii) = lmmResults(1).fixedSlopes{mm,ii};
                modelseMat(:,ii) = lmmResults(1).fixedSlopes_SE{mm,ii};
            end
        end
    end

    for c = conditionOrderSimSeq
        subplot(1,4,subplotOrder(c)); hold all;
        plot([0.2 9.7],[1 1],'k--');
        cData = repmat(dataMat(c,:),2,1);
        seData = repmat(seMat(c,:),2,1);
        
        if c==1, ylabel('Suppression (fitted regression slope)'); end
        
        for idx = 1:length(roisToPlot)
            % Figure 4B
            if ~plotModelAmpl && strcmp(LMMOrder{mm},'Data')
                curMnSlope  = lmmResults.fixedSlopes{idx}(c);
                SE          = lmmResults.fixedSlopes_SE{idx}(c);
                curIndivSubj = lmmResults.subjSlopes{idx}(:,c);
                
                scatter(idx+(0.01*rand(1,length(curIndivSubj))),curIndivSubj,30,[0.5 0.5 0.5], 'filled', ...
                    'MarkerEdgeAlpha',0.3, 'MarkerFaceAlpha',0.3,'MarkerEdgeColor',[0 0 0]); hold on;
                plot([idx, idx], curMnSlope + [-SE, SE],'k');
                scatter(idx,curMnSlope,60,cmap(idx,:), 'filled','MarkerEdgeColor',[0 0 0]); hold on;
                
 
            else % Figure 7
                if strcmp(LMMOrder{mm},'Data')
                    errorbar3(x(:,idx)',cData(:,idx)',seData(:,idx)',1,[0.7 0.7 0.7]); % cmapROIs(jj,:)
                    ax = gca;
                    ax.Children(1).FaceAlpha = 0.3;                  
                else
                    curMnSlope   = lmmResults(1).fixedSlopes{mm,idx}(c);
                    SE           = lmmResults(1).fixedSlopes_SE{mm,idx}(c);
                    if mm==1 && ~strcmp(LMMOrder{mm},'suppl_DoG')
                        scatter(idx,curMnSlope,100,cmap(mm,:),'o','filled','linewidth',0.15,'MarkerFaceColor',cmap(mm,:)); hold on; [0,0,0]
                    else
                        scatter(idx,curMnSlope,100,cmap(mm,:),'o','filled','linewidth',0.15); hold on; %'MarkerEdgeColor',[0,0,0]
                    end
                    plot([idx, idx], curMnSlope + [-SE, SE],'color',cmap(mm,:), 'linewidth',2);
                    plot([idx, idx], curMnSlope + [-SE, SE],'color',[0 0 0], 'linewidth',0.5);
                end
            end
        end
        plot([0.2 9.7],[1 1],'k--');
        title(sprintf('%s',conditionNamesSimSeq{c}))
        xlim([0.2 9.7]); ylim([0 1.31]); box off
        set(gca,'XTick',1:length(roisToPlot), 'XTickLabel',string(roisToPlot),'XTickLabelRotation',45)
    end
end

% Add legend and title
l = gca;
if plotModelAmpl
    sgtitle('Group Average Data vs Model: LMM regression slopes (Mean +/-SE)');
    if strcmp(LMMOrder{1},'DoG')
        legend(l.Children([length(l.Children):-8:1]),LMMOrder{1:2}, 'FontSize',9, 'Location','SouthWest'); legend boxoff
    else
        legend(l.Children([length(l.Children):-30:1]),LMMOrder{1:3}, 'FontSize',9, 'Location','SouthWest'); legend boxoff
    end
else
    legend(l.Children([length(l.Children)-3:-3:1]),string(roisToPlot), 'FontSize',9, 'Location','SouthWest'); legend boxoff
    sgtitle('Group Average Data: LMM regression slopes (Mean +/-SE)')
end

if saveFigs
    if plotModelAmpl
        if strcmp(LMMOrder{1},'suppl_DoG')
            fName = sprintf('SummarySlopeGroup_SEQvsSIMMODELAmpl%s_LMMfit',cell2mat(LMMOrder));
            subDir = 'SupplFig7';
        else
            fName = sprintf('SummarySlopeGroup_SEQvsSIMMODELAmpl%s_LMMfit',cell2mat(LMMOrder));
            subDir = 'fig7';
        end
    else
        fName = sprintf('SummarySlopeGroup_SEQvsSIMAmpl_LMMfit');
        subDir     = 'fig4';
    end
    
    % Set paths
    if ~exist('saveFigDir','var')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose')
end
