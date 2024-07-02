function fH = makeSupplFigure4_DistributionCSTExponents(ds,roisToPlot,cmapROIs,saveFigs, saveFigDir)
% Function to reproduce supplementary figure 4: 
% Plotting the distribution of fitted CST pRF exponents
%
% From the paper:
% Title:   Rethinking simultaneous suppression in visual cortex via 
%          compressive spatiotemporal population receptive fields.
% Authors: Kupers, Kim, Grill-Spector (2024).
% Journal: Nature Communications
% DOI:     XXX
%
% Requires getting MRI data from OSF (see downloadDataTableFromOSF.m)
%
% Code written by E.R. Kupers (2024) Stanford University
% 
% INPUTS (required):
% - ds              : dataset
% - roisToPlot      : cell with ROI names
% - cmapROIs        : color map for ROIs
% - saveFigs        : save figures or not?
% - saveFigDir      : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle
%
%% Define distribution params
nbins   = 20;
nboot   = 1000;
binsExp = linspace(0.075,1.025,nbins);
binsExpShifted = binsExp(2:end)-0.025;
x_upsampled = linspace(0.1,1,38);

% Define data params
subjnrs     = double(unique(ds.Subject)');

% Allocate space
allExp_hist  = NaN(length(subjnrs),length(roisToPlot),nbins-1);
allExp_boots = NaN(length(subjnrs),length(roisToPlot),nboot);
allExp_max   = NaN(length(subjnrs),length(roisToPlot));

for sj = 1:length(subjnrs)
    for idx = 1:length(roisToPlot)
        
        % Get PRF CST size and exp
        expCST = ds.pRFCSTexp(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        
        if length(expCST)>1
            % Resample data
            resampledExpCST = randsample(expCST, nboot, true);
            
            % Take histogram of resampled data
            h = histogram(resampledExpCST,binsExp, ...
                'Normalization', 'probability','Visible','off');
            
            % Get max exponent
            [~,max_idx] = max(h.Values);
            binsExpShifted = binsExp(2:end)-0.025;
            mode_exp = binsExpShifted(max_idx);
            
            % Store data
            allExp_hist(sj,idx,:) = h.Values;
            allExp_boots(sj,idx,:) = resampledExpCST;
            allExp_max(sj,idx) = mode_exp;
            h.Visible = 'off';  clear h;
        end
    end
end      


%% Get group average across subjects
singleSubjExpHistogram = NaN(length(subjnrs),length(roisToPlot),length(binsExp)-1);
for sj = 1:length(subjnrs)
    for idx = 1:length(roisToPlot)
        [N, ~] = histcounts(squeeze(allExp_boots(sj,idx,:)),...
            binsExp, 'Normalization','Probability');
        singleSubjExpHistogram(sj,idx,:) = N;
    end
end

%% Plot group average
fH = figure; set(gcf,'Position',[500 1 800 687],'color','w');
for idx = 1:length(roisToPlot)
 
    nrSubjWithData = length(~isnan(allExp_boots(:,idx,1)));
    mdGroupExp(idx) = median(reshape(allExp_boots(:,idx,:),[],1),1,'omitnan');
%     sdGroupExp(idx) = std(reshape(allExp_boots(:,idx,:),[],1),[],1,'omitnan');
%     ci68GroupExp(idx,:) = prctile(reshape(allExp_boots(:,idx,:),[],1),[16, 84]);

      mnGroupExp = squeeze(mean(singleSubjExpHistogram(:,idx,:),1,'omitnan'));
      semGroupExp = squeeze(std(singleSubjExpHistogram(:,idx,:),[],1,'omitnan'))./sqrt(nrSubjWithData);
    
    if ~isempty(mnGroupExp) && all(~isnan(mnGroupExp))
        subplot(3,3,idx); cla; hold all;
        y_upsampled  = interp1(binsExpShifted,mnGroupExp, x_upsampled,'pchip');
        se_upsampled = interp1(binsExpShifted,semGroupExp, x_upsampled,'pchip');
        [~,max_idx]  = max(mnGroupExp);
        mode_exp(idx) = binsExpShifted(max_idx);
        mode_exp_sem(idx) = semGroupExp(max_idx);
        
        hold all;
        shadedErrorBar(x_upsampled,y_upsampled,se_upsampled,'lineProps',{'Color',cmapROIs(idx,:)})
%         plot(mode_exp(idx),0.4-(idx*0.01),'v','color',cmapROIs(idx,:))
        plot(mdGroupExp(idx),0.3-(idx*0.01),'*','color',cmapROIs(idx,:))

    end
    title(string(roisToPlot(idx)));
    set(gca, 'XTick',[0:0.2:1],'FontSize',10);
    ylabel('Probability');
    xlabel('Exponent');
    box off;
    set(gca,'xlim', [0.05 1.05], 'ylim',[0 0.5]);
end
sgtitle(sprintf('Group summary gridfit CST exponent (mean +/- SEM), with mode (triangle) and median (*)'))

% Save figure
if saveFigs
    fName = sprintf('supplfig4_SummaryGroupMedian_GridFitExp');
    thisSaveFigDir = fullfile(saveFigDir, 'supplfig4');
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%     print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end

return

