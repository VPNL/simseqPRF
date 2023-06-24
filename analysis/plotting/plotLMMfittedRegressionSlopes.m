function fH = plotLMMfittedRegressionSlopes(lmmResults,LMMOrder, roisToPlot,cmapROIs, saveFigs)
% Function to plot Linear Mixed Model fitted regression slopes (fixed and
% random subject slopes, for each condition and visual area).
% To plot Figure 4B, "LMMResults" struct has fields with size {1,nrROIs}
%   ex: plotLMMfittedRegressionSlopes(LMMResults,{'Data'}, false)
% To plot Figure 7, "LMMResults" struct has fields with size {4,nrROIs}
%   ex: plotLMMfittedRegressionSlopes(LMMResults,{'LSS','CSS','CST','Data'}, false)

%% Check inputs
if size(lmmResults,1) == 1 && size(lmmResults.fixedSlopes,1)==1
    plotModelAmpl = false;
else
    plotModelAmpl = true;
end
if ~exist('saveFigs','var') || isempty(saveFigs)
    saveFigs = false; % Save figures or not
end

% Set params
conditionOrderSimSeq = 1:4;
subplotOrder         = [2,1,4,3];
conditionNamesSimSeq = {'Small & short (0.2s)','Small & long (1s)','Big & Short (0.2s)','Big & Long (1s)'};
% roisToPlot           = {'V1','V2','V3','hV4','VO1VO2','V3AB','IPS0IPS1','LO1LO2','TO1TO2'}; 
% cmapROIs             = getROISummaryColors(0);
% newROIOrder          = [1,2,3,4,5,8,9,6,7];
% cmapROIs             = cmapROIs(newROIOrder,:);

if plotModelAmpl
    markerTypes = {'o','o','o','+'}; % ModelOrder and cross for Data
    darkgray   = [0.5 0.5 0.5];
    orange     = [255 140 0]./255;
    blue       = [26 115 225]./255;
    black = [0 0 0];
    cmap       = [darkgray;orange;blue];
else
    markerTypes = {'o','o','o','o'};
    cmap = cmapROIs;
end

%% Plot it!

fH = figure; clf; set(gcf,'Position',[669 30 1342 905]);
dataMat = reshape(cell2mat(lmmResults(end).fixedSlopes),4,[]);
seMat = reshape(cell2mat(lmmResults(end).fixedSlopes_SE),4,[]);
x = [1:9]+[-0.4,0.4]';
for mm = 1:length(LMMOrder)

    for c = conditionOrderSimSeq
        subplot(1,4,subplotOrder(c)); hold all;
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
                    if mm==1
                        scatter(idx,curMnSlope,100,cmap(mm,:),markerTypes{mm},'linewidth',0.15); hold on; %'MarkerEdgeColor',[0,0,0]
                    else
                        scatter(idx,curMnSlope,100,cmap(mm,:),markerTypes{mm}, 'filled','linewidth',0.15); hold on; %'MarkerEdgeColor',[0,0,0]
                    end
                    plot([idx, idx], curMnSlope + [-SE, SE],'color',cmap(mm,:), 'linewidth',2);
                    plot([idx, idx], curMnSlope + [-SE, SE],'color',[0 0 0], 'linewidth',0.5);
                end
            end
        end
        title(sprintf('%s',conditionNamesSimSeq{c}))
        xlim([0.2 9.7]); ylim([0 1.31]); box off
        set(gca,'XTick',1:9, 'XTickLabel',string(roisToPlot),'XTickLabelRotation',45)
    end
end

if plotModelAmpl
    for ii = 1:4
        subplot(1,4,ii); ax = gca;
        for kk = 1:3:(3*9)
            ax.Children(kk).ZData = -10;
        end
    end
end

% Add legend and title
l = gca;
if plotModelAmpl
    legend(l.Children([length(l.Children):-30:1]),LMMOrder, 'FontSize',9, 'Location','SouthWest'); legend boxoff
    sgtitle(sprintf('Group Average Data vs Model: LMM regression slopes (Mean +/-SE)'))
else
    legend(l.Children([length(l.Children):-3:1]),string(roisToPlot), 'FontSize',9, 'Location','SouthWest'); legend boxoff
    sgtitle('Group Average Data: LMM regression slopes (Mean +/-SE)')
end

if saveFigs
    if plotModelAmpl
        fName = sprintf('SummarySlopeGroup_SEQvsSIMMODELAmpl%s_LMMfit',cell2mat(LMMOrder));
        subDir     = 'fig7';
    else
        fName = sprintf('SummarySlopeGroup_SEQvsSIMAmpl_LMMfit');
        subDir     = 'fig4';
    end
    
    % Set paths
    saveFigDir = fullfile(simseqRootPath,'results','group');
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose')
end
