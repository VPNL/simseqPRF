function fH = plotLMMfittedRegressionSlopesCummulative(lmmResults, roisToPlot,cmapROIs, saveFigs)
% Function to plot the cummuluative effect of suppression slopes across the
% visual hierarchy (Figure 4C)

%% Check inputs
plotModelAmpl = false;

if ~exist('saveFigs','var') || isempty(saveFigs)
    saveFigs = false; % Save figures or not
end

% Set params
conditionOrderSimSeq = 1:4;
subplotOrder         = [2,1,4,3];
conditionNamesSimSeq = {'Small & short (0.2s)','Small & long (1s)','Big & Short (0.2s)','Big & Long (1s)'};
markerTypes = {'o','o','o','o'};
cmap = cmapROIs;


%% Plot it!

fH = figure; clf; set(gcf,'Position',[669 30 1342 905]);
dataMat = reshape(cell2mat(lmmResults(end).fixedSlopes),4,[]);
seMat = reshape(cell2mat(lmmResults(end).fixedSlopes_SE),4,[]);

inheritIdx = [2,1; ... V2-V1
                      3,2; ... V3-V2
                      4,3; ... hV4-V3
                      5,4; ... VO-hV4
                      6,3; ... V3AB-V3
                      7,6; ... IPS-V3AB
                      8,3; ... LO-V3
                      9,8]; %  TO-LO
diffLabels = {'V2-V1','V3-V2','hV4-V3','VO-hV4','V3AB-V3','IPS-V3AB','LO-V3','TO-LO'};                      
dataMatGroup = dataMat(:,inheritIdx(:,1))-dataMat(:,inheritIdx(:,2));
dataMatSE = [];
dataMatSubj = [];
for ii = 1:size(inheritIdx,1)   
    subjsArea1 = unique(ds.Subject(ismember(ds.ROI,roisToPlot(inheritIdx(ii,1)))));
    subjsArea2 = unique(ds.Subject(ismember(ds.ROI,roisToPlot(inheritIdx(ii,2)))));
    dataArea1 = NaN(10,4); dataArea2 = dataArea1;
    dataArea1(subjsArea1,:) = lmmResults.subjSlopes{inheritIdx(ii,1)};
    dataArea2(subjsArea2,:) = lmmResults.subjSlopes{inheritIdx(ii,2)};
    dataMatSubj{ii} = dataArea1 - dataArea2;
    dataMatSE{ii} = std(dataMatSubj{ii},[],1, 'omitnan')./sqrt(sum(~isnan(dataMatSubj{ii}(:,1))));
end


for c = conditionOrderSimSeq
    subplot(1,4,subplotOrder(c)); hold all;

    if c==1, ylabel('Suppression (fitted regression slope)'); end
    
    for idx = 1:length(roisToPlot)-1
        if idx == 1
            plot([0.2 9.7], [0 0],'k:', 'linewidth', 0.5)
        end
        % Figure 4B
        curMnSlope  = dataMatGroup(c,idx);
        SE          = dataMatSE{idx}(c);
        curIndivSubj = dataMatSubj{idx}(:,c);
        
        scatter(idx+(0.01*rand(1,length(curIndivSubj))),curIndivSubj,30,[0.5 0.5 0.5], 'filled', ...
            'MarkerEdgeAlpha',0.3, 'MarkerFaceAlpha',0.3,'MarkerEdgeColor',[0 0 0]); hold on;
        plot([idx, idx], curMnSlope + [-SE, SE],'k');
        scatter(idx,curMnSlope,60,cmap(idx+1,:), 'filled','MarkerEdgeColor',[0 0 0]); hold on;
    end
    title(sprintf('%s',conditionNamesSimSeq{c}))
    xlim([0.2 9.7]); ylim([-1 1]); box off
    set(gca,'XTick',1:8, 'XTickLabel',diffLabels,'XTickLabelRotation',45)
end


% Add legend and title
l = gca;
legend(l.Children([(length(l.Children)-3):-3:1]),string(diffLabels), 'FontSize',9, 'Location','SouthWest'); legend boxoff
sgtitle('Group Average Data: LMM regression slopes (Mean +/-SE)')


if saveFigs
    fName = sprintf('SummarySlopeGroup_SEQvsSIMAmpl_LMMfit_DeltaVisualAreas');
    subDir     = 'SupplFig7';

    
    % Set paths
    saveFigDir = fullfile(simseqRootPath,'results','group');
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose')
end
