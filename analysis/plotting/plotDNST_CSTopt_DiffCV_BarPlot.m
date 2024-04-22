function fH = plotDNST_CSTopt_DiffCV_BarPlot(params, roisToPlot, cmapModels, saveFigs, saveFigDir)

% Preallocate space
mean_diffCVR2_acrossVoxels = NaN(7,length(roisToPlot),length(params));
maxNC              = mean_diffCVR2_acrossVoxels;
diffCVR2_allVoxels = cell(length(roisToPlot),length(params));
ipsIdx = cellfind(regexp(string(roisToPlot),'IPS0/IPS1','once'));
roisNoIPS = setdiff([1:length(roisToPlot)],ipsIdx);

% Compute group average difference CV-R^2
for sj = 1:7
    for idx = roisNoIPS

        diffCVR2{sj,idx,1} = params{1}.R2CST{sj,idx} - params{2}.R2CST{sj,idx}; % CST fix - CST opt
        diffCVR2{sj,idx,2} = params{1}.R2CST{sj,idx} - params{3}.R2DNST{sj,idx}; % CST fix - DNST
        diffCVR2{sj,idx,3} = params{2}.R2CST{sj,idx} - params{3}.R2DNST{sj,idx}; % CST opt - DNST
        
        % mean across voxels within subject's ROI
        mean_diffCVR2_acrossVoxels(sj,idx,1) = mean(diffCVR2{sj,idx,1},'omitnan');
        mean_diffCVR2_acrossVoxels(sj,idx,2) = mean(diffCVR2{sj,idx,2},'omitnan');
        mean_diffCVR2_acrossVoxels(sj,idx,3) = mean(diffCVR2{sj,idx,3},'omitnan');
         
        if ~isempty(params{1}.NC{sj,idx})
            maxNC(sj,idx,1) = max(params{1}.NC{sj,idx}, [], 'omitnan');
        end
        if ~isempty(params{2}.NC{sj,idx})
            maxNC(sj,idx,2) = max(params{2}.NC{sj,idx}, [], 'omitnan');
        end
        if ~isempty(params{3}.NC{sj,idx})
            maxNC(sj,idx,3) = max(params{3}.NC{sj,idx}, [], 'omitnan');
        end
    end
end

%% Set up figure params
barWidth   = 0.75;
xposBar    = [1:length(roisToPlot)];
xposError  = [-0.215,0,0.215]; % X position for errorbars
modelDiffNames = {'CST_{fix}-CST_{opt}','CST_{fix}-DN-ST','CST_{opt}-DN-ST'};

% Summarize group data (mean and SEM)
for idx = 1:length(roisToPlot)
    for mm = 1:3
        nrSubjWithData        = sum(~isnan(mean_diffCVR2_acrossVoxels(:,idx,mm)));
        mnDiffR2Group(idx,mm) = mean(mean_diffCVR2_acrossVoxels(:,idx,mm),1,'omitnan');
        seDiffR2Group(idx,mm) = std(mean_diffCVR2_acrossVoxels(:,idx,mm),[],1,'omitnan')./sqrt(nrSubjWithData);
    end
end

% Start plotting
fH(1) = figure; clf; set(gcf,'Position',[500 1 965 950],'color','w'); hold all;
b = bar(xposBar,mnDiffR2Group, barWidth); hold on;
b(1).FaceColor = cmapModels(1,:);
b(2).FaceColor = cmapModels(1,:);
b(3).FaceColor = cmapModels(2,:);

b(1).EdgeColor = cmapModels(2,:);
b(2).EdgeColor = cmapModels(3,:);
b(3).EdgeColor = cmapModels(3,:);
b(1).LineWidth = 5;
b(2).LineWidth = 5;
b(3).LineWidth = 5;

% Plot errorbars
errorbar2(xposBar'+xposError,mnDiffR2Group, seDiffR2Group, 1,'k'); hold on;

% Add individual subjects
for idx = 1:length(roisToPlot)
    for mm = 1:3
        currIndivSubj = mean_diffCVR2_acrossVoxels(:,idx,mm);
        scatter(idx+xposError(mm)+(0.02*randn(1,size(mean_diffCVR2_acrossVoxels,1)))',currIndivSubj,30,cmapModels(mm,:), 'filled', ...
            'MarkerEdgeAlpha',0.7, 'MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0 0 0]); hold on;
    end
end

% Add axis labels, ticks, legend, and title
xlim([0 10]); ylim([-8 10]); box off
set(gca,'XTick',[1:length(roisToPlot)], 'XTickLabel',string(roisToPlot),'XTickLabelRotation',45)
l = gca;
legend(l.Children([length(l.Children):-1:1]),modelDiffNames, 'FontSize',9, 'Location','NorthWest'); 
legend boxoff
sgtitle('Group Average cv-R^2 mean +/- SEM per visual area')
ylabel('Cross-validated variance explained (%)')

% Save if requested
if saveFigs
    fName = sprintf('SupplFig8C_GroupR2_STRet_diff_CSTfix_CSTopt_DNST');   
    subDir = 'SupplFig8';
    if ~exist('saveFigDir','dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end

    


