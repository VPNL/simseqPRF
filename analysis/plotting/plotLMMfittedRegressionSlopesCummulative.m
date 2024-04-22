function fH = plotLMMfittedRegressionSlopesCummulative(deltasGroup,deltasSubj, ...
        deltasSE,deltaLabels,cmapROIs, saveFigs,saveFigDir)
% Function to plot the cummuluative effect of suppression slopes across the
% visual hierarchy (Suppl Figure 2A)

%% Check inputs

if ~exist('saveFigs','var') || isempty(saveFigs)
    saveFigs = false; % Save figures or not
end

if ~exist('saveFigDir','var') || isempty(saveFigs)
    saveFigDir = fullfile(simseqRootPath, 'results','group');
end


% Set params
conditionOrderSimSeq = 1:4;
subplotOrder         = [2,1,4,3];
conditionNamesSimSeq = {'Small & Short (0.2s)','Small & Long (1s)',...
                        'Big & Short (0.2s)','Big & Long (1s)'};

%% Plot it!

fH = figure; clf; set(gcf,'Position',[669 30 1342 905]);

for c = conditionOrderSimSeq
    subplot(1,4,subplotOrder(c)); hold all;

    if c==1, ylabel('Suppression (fitted regression slope)'); end
    
    for idx = 1:size(deltasGroup,2)
        if idx == 1
            plot([0.2 9.7], [0 0],'k:', 'linewidth', 0.5)
        end
        
        curMnSlope  = deltasGroup(c,idx);
        SE          = deltasSE{idx}(c);
        curIndivSubj = deltasSubj{idx}(:,c);
        
        scatter(idx+(0.01*rand(1,length(curIndivSubj))),curIndivSubj,30,[0.5 0.5 0.5], 'filled', ...
            'MarkerEdgeAlpha',0.3, 'MarkerFaceAlpha',0.3,'MarkerEdgeColor',[0 0 0]); hold on;
        plot([idx, idx], curMnSlope + [-SE, SE],'k');
        scatter(idx,curMnSlope,60,cmapROIs(idx+1,:), 'filled','MarkerEdgeColor',[0 0 0]); hold on;
    end
    title(sprintf('%s',conditionNamesSimSeq{c}))
    xlim([0.2 9.7]); ylim([-1 1]); box off
    set(gca,'XTick',1:8, 'XTickLabel',deltaLabels,'XTickLabelRotation',45)
end


% Add legend and title
l = gca;
legend(l.Children([(length(l.Children)-3):-3:1]),string(deltaLabels), 'FontSize',9, 'Location','SouthWest'); legend boxoff
sgtitle('Group Average Data: LMM regression slopes (Mean +/-SE)')

% Save figure if requested
if saveFigs
    fName = sprintf('SummarySlopeGroup_SEQvsSIMAmpl_LMMfit_DeltaVisualAreas');
    subDir = 'SupplFig2';

    % Set paths
    if ~exist('saveFigDir','dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
    % print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose')
end

