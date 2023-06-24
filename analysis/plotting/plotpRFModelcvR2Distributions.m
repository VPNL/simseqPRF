function fH = plotpRFModelcvR2Distributions(diffR2Resampled)

r2bins = -100:5:100;
clear histR2 h
for sj = 1:length(subjnrs)
    for idx = 1:length(roisToPlot)
        for mm = 1:3
            if sum(isnan(diffR2Resampled(sj,idx,mm,:)))~=length(diffR2Resampled(sj,idx,mm,:))      
                h = histogram(squeeze(diffR2Resampled(sj,idx,mm,:)),r2bins);
                histR2(sj,idx,mm,:) = h.Values;
            else
                histR2(sj,idx,mm,:) = NaN(1,1,1,length(r2bins)-1);
            end
            
        end
    end
end
%%
clear mnGroupR2 semGroupR2
binsExpShifted = r2bins(2:end);
x_upsampled = r2bins;

fH1 = figure(15); clf; set(gcf,'Position',[500 1 1965 687],'color','w');
for idx = 1:length(roisToPlot)
    subplot(3,3,idx); cla; hold all;
    mnSubjR2(mm,idx,:) = squeeze(mean(diffR2Resampled(:,idx,mm,:),4,'omitnan'));
    
    for mm = 1:3
        nrSubjWithData = sum(isnan(subjR2Resampled(:,idx,mm,:)),4)==0;
       mnGroupR2(mm,idx,:) = squeeze(std(histR2(:,idx,mm,:),[],4,'omitnan'));
       semGroupR2(mm,idx,:) = squeeze(std(histR2(:,idx,mm,:),[],4,'omitnan'))./sqrt(sum(nrSubjWithData));
    
        if ~isempty(mnGroupR2(mm,idx,:)) && all(~isnan(mnGroupR2(mm,idx,:)))
            y_upsampled = interp1(binsExpShifted,squeeze(mnGroupR2(mm,idx,:)), x_upsampled,'pchip');
            se_upsampled = interp1(binsExpShifted,squeeze(semGroupR2(mm,idx,:)), x_upsampled,'pchip');
            [~,max_idx] = max(mnGroupR2(mm,idx,:));
            max_r2(mm,idx) = binsExpShifted(max_idx);
            max_r2_sem(mm,idx) = semGroupR2(max_idx);
            
            hold all;
            shadedErrorBar(x_upsampled,y_upsampled,se_upsampled,'lineProps',{'Color', cmapBetas(mm,:)})

        end
    end
    
    if ismember(idx,[1,4,7]),
        ylabel('Probability (fraction)');
    end
    if ismember(idx,[7:9]),
        xlabel('Variance explained (%)');
    end
    box off;
        set(gca,'xlim', [-50 50]);%, 'ylim',[0 0.04]);
end

subplot(3,3,1)
legend(modelNames,'Location', 'none','Position', [0.8 0.3 0.1 0.10], 'FontSize',10)
legend boxoff

fName = sprintf('CrossValR2_GroupAverage_SEM');
sgtitle(sprintf('Group cross-validated R^2 (mean +/- SEM) (N=%d)', length(subjnrs)));

if saveFigs
    thisSaveFigDir = fullfile(saveFigDir, 'R2_group_model_compare_distributions');
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
    print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end