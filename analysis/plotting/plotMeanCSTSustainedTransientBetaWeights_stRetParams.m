function fH = plotMeanCSTSustainedTransientBetaWeights_stRetParams(allDS, ...
    roisToPlot, saveFigs, saveFigDir)

% Define params
roisNoIPS = [1:6,8,9];
roisToLabel = roisToPlot(roisNoIPS);

nrSubjects  = length(unique(allDS{1}.ds.Subject));
cst_s       = NaN(nrSubjects,length(roisNoIPS));
cst_t       = cst_s;
cst_s_opt   = cst_s;
cst_t_opt   = cst_t;
grandMn_S_diff = NaN(1,length(roisNoIPS));
grandSE_S_diff = grandMn_S_diff;
grandMn_T_diff = grandMn_S_diff;
grandSE_T_diff = grandMn_S_diff;
nboot       = 1000;

% Use two-tone to avoid bar graph nightmares
cmapBetas   = [[0.3 0.3 0.3];[1 1 1]];

for idx = roisNoIPS
    
    for sj = 1:nrSubjects
        
        % CST fix
        betavals_s = allDS{1}.ds.BetaCST_s(allDS{1}.ds.Subject==sj & ...
            allDS{1}.ds.ROI==nominal(roisToPlot(idx)) & ...
            allDS{1}.ds.Condition==nominal(1));
        betavals_t = allDS{1}.ds.BetaCST_t(allDS{1}.ds.Subject==sj & ...
            allDS{1}.ds.ROI==nominal(roisToPlot(idx)) & ...
            allDS{1}.ds.Condition==nominal(1));
        
        if ~isempty(betavals_s) || length(betavals_s) > 1
            betavals_s_resampled = randsample(betavals_s, nboot, true);
            betavals_t_resampled = randsample(betavals_t, nboot, true);
            
            cst_s(sj,idx) = mean(betavals_s_resampled,'omitnan');
            cst_t(sj,idx) = mean(betavals_t_resampled,'omitnan');
            
        else
            cst_s(sj,idx) = NaN;
            cst_t(sj,idx) = NaN;
        end
        
        % CST opt
        betavals_s_opt = allDS{2}.ds.BetaCST_ST_s(allDS{2}.ds.Subject==sj & ...
            allDS{2}.ds.ROI==nominal(roisToPlot(idx)) & ...
            allDS{2}.ds.Condition==nominal(1));
        betavals_t_opt = allDS{2}.ds.BetaCST_ST_t(allDS{2}.ds.Subject==sj & ...
            allDS{2}.ds.ROI==nominal(roisToPlot(idx)) & ...
            allDS{2}.ds.Condition==nominal(1));
        
        if ~isempty(betavals_s_opt) || length(betavals_s_opt) > 1
            betavals_s_resampled_opt = randsample(betavals_s_opt, nboot, true);
            betavals_t_resampled_opt = randsample(betavals_t_opt, nboot, true);
            
            cst_s_opt(sj,idx) = mean(betavals_s_resampled_opt,'omitnan');
            cst_t_opt(sj,idx) = mean(betavals_t_resampled_opt,'omitnan');
            
        else
            cst_s_opt(sj,idx) = NaN;
            cst_t_opt(sj,idx) = NaN;
        end
    end
    grandMn_S_diff(idx) = squeeze(mean(cst_s(:,idx)-cst_s_opt(:,idx),1,'omitnan'));
    grandSE_S_diff(idx) = squeeze(std(cst_s(:,idx)-cst_s_opt(:,idx),[],1,'omitnan'))./sqrt(sum(~isnan(cst_s(:,idx)-cst_s_opt(:,idx))));
    grandMn_T_diff(idx) = squeeze(mean(cst_t(:,idx)-cst_t_opt(:,idx),1,'omitnan'));
    grandSE_T_diff(idx) = squeeze(std(cst_t(:,idx)-cst_t_opt(:,idx),[],1,'omitnan'))./sqrt(sum(~isnan(cst_t(:,idx)-cst_t_opt(:,idx))));
    
end

% Plot bar graph
fH = figure; set(gcf,'Position',[1 1 800 941]); hold all;

b = bar([1:9],[grandMn_S_diff;grandMn_T_diff],'FaceColor','flat'); hold all;
b(1).CData = cmapBetas(1,:);
b(2).CData = cmapBetas(2,:);
eb1 = errorbar([1:9]-0.15,grandMn_S_diff, grandSE_S_diff, grandSE_S_diff,'k','linewidth',0.5);
eb1.LineStyle = 'none';
eb1.CapSize = 0;
eb2 = errorbar([1:9]+0.15,grandMn_T_diff, grandSE_T_diff, grandSE_T_diff,'k','linewidth',0.5);
eb2.LineStyle = 'none';
eb2.CapSize = 0;

title('Diff CST_{fix} - CST_{opt} Beta values for Sustained & Transient channels');
set(gca,'XTickLabel',string(roisToPlot)', 'XTickLabelRotation',45)
set(gca,'FontSize',12, 'TickDir','out')
box off;
ylabel('Beta weight difference CST_{fix} - CST_{opt} (a.u.)')
ylim([-0.2 0.7])
legend('Sustained','Transient'); legend box off

if saveFigs
    fName = sprintf('SummaryGroupAverageSEM_BetaCSTSustainedTransientDiffCSTfix_v_opt');
    subDir = 'SupplFig8';
    if ~exist('saveFigDir','dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
    print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
end

%% Do some stats (2-Way repeated measured ANOVA)
s = cst_s(1:7,roisNoIPS);
s_opt = cst_s_opt(1:7,roisNoIPS);  % 7 subjects x 8 ROIs
t = cst_t(1:7,roisNoIPS);
t_opt = cst_t_opt(1:7,roisNoIPS); % 7 subjects x 8 ROIs

subj_vec = repmat([1:7],8,1);
roi_vec = repmat(string(roisToLabel),[7,1]);
T = table(repmat(subj_vec(:),4,1),...
    repmat(roi_vec,4,1),[ones(7*8,1); 2.*ones(7*8,1); ones(7*8,1); 2.*ones(7*8,1)],...
    [ones(7*8*2,1); 2.*ones(7*8*2,1)],[s(:); t(:);s_opt(:);t_opt(:)],...
    'VariableNames',{'Subject','ROI','ModelType','ChannelType', 'Meas'});
rm = fitrm(T,'Meas ~ ModelType * ChannelType * ROI');
disp(rm.anova);