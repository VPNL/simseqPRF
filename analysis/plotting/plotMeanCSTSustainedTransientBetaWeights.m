function fH = plotMeanCSTSustainedTransientBetaWeights(mean_resampledBetavalSust, ...
    mean_resampledBetavalTrans,roisToPlot,saveFigs,saveFigDir)

% Set colors
roiColors = getROISummaryColors(0);
newROIOrder = [1,2,3,4,5,8,9,6,7];
roiColors = roiColors(newROIOrder,:);
white = [1 1 1];

% Compute grand mean + SEM across subjects
for idx = 1:length(roisToPlot)
    grandMn_S(idx) = squeeze(mean(mean_resampledBetavalSust(:,idx),1,'omitnan'));
    grandSE_S(idx) = squeeze(std(mean_resampledBetavalSust(:,idx),[],1,'omitnan'))./sqrt(sum(~isnan(mean_resampledBetavalSust(:,idx))));
    grandMn_T(idx) = squeeze(mean(mean_resampledBetavalTrans(:,idx),1,'omitnan'));
    grandSE_T(idx) = squeeze(std(mean_resampledBetavalTrans(:,idx),[],1,'omitnan'))./sqrt(sum(~isnan(mean_resampledBetavalTrans(:,idx))));

    indiv_S(idx,:) = squeeze(mean_resampledBetavalSust(:,idx));
    indiv_T(idx,:) = squeeze(mean_resampledBetavalTrans(:,idx)); 
end

fH = figure; clf; set(gcf,'Position',[441   195   643   602]);
b = bar([1:length(roisToPlot)],[grandMn_S;grandMn_T],'FaceColor','flat'); hold all;
b(1).CData = roiColors;
b(2).CData = white;
eb1 = errorbar([1:length(roisToPlot)]-0.15,grandMn_S, grandSE_S, grandSE_S,'k','linewidth',0.5);
eb1.LineStyle = 'none';
eb1.CapSize = 0;
eb2 = errorbar([1:length(roisToPlot)]+0.15,grandMn_T, grandSE_T, grandSE_T,'k','linewidth',0.5);
eb2.LineStyle = 'none';
eb2.CapSize = 0;

% Plot individual subjects
for idx = 1:length(roisToPlot)
    scatter(idx-0.15+(0.01*rand(1,size(indiv_S,2))),indiv_S(idx,:),60,[0.5 0.5 0.5], 'filled', ...
        'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerFaceAlpha',0.8,'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.8);
    scatter(idx+0.15+(0.01*rand(1,size(indiv_T,2))),indiv_T(idx,:),60,[0.5 0.5 0.5], 'filled', ...
        'MarkerFaceColor',[1 1 1], 'MarkerFaceAlpha',0.8,'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.8);
end

title('Beta values for CST Sustained & Transient channels');
set(gca,'XTickLabel',string(roisToPlot)', 'XTickLabelRotation',45)
set(gca,'FontSize',12, 'TickDir','out')
box off;
ylabel('Beta weight (% signal change)')
ylim([0 2.8])
legend('Sustained','Transient'); legend box off

if saveFigs
    if ~exist(saveFigDir,'dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    fName = sprintf('SummaryGroupAverageSEM_BetaCSTSustainedTransient_resampled');
    thisSaveFigDir = fullfile(saveFigDir,'fig8');
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
    print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose')
end

%% Do some stats (2-Way repeated measured ANOVA)
% Test differences between CST_opt sustained vs CST_opt transient channel;
subj_vec = repmat([1:size(indiv_S,2)],length(roisToPlot),1);
roi_vec = repmat(string(roisToPlot),[size(indiv_S,2),1]);
T1 = table(repmat(subj_vec(:),2,1),...
    repmat(roi_vec,2,1),...
    [ones(length(roisToPlot)*size(indiv_S,2),1); 2.*ones(length(roisToPlot)*size(indiv_S,2),1)],...
    [indiv_S(:); indiv_T(:)],...
    'VariableNames',{'Subject','ROI','ChannelType', 'Meas'});
rm = fitrm(T1,'Meas ~ ChannelType * ROI');
disp(rm.anova);