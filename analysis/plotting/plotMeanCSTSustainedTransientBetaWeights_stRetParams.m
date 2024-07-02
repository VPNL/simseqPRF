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
grandMn_S = NaN(1,length(roisNoIPS));
grandSE_S = grandMn_S;
grandMn_T = grandMn_S;
grandSE_T = grandMn_S;
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
    
    grandMn_S(idx) = squeeze(mean(cst_s_opt(:,idx),1,'omitnan'));
    grandSE_S(idx) = squeeze(std(cst_s_opt(:,idx),[],1,'omitnan'))./sqrt(sum(~isnan(cst_s_opt(:,idx))));
    grandMn_T(idx) = squeeze(mean(cst_t_opt(:,idx),1,'omitnan'));
    grandSE_T(idx) = squeeze(std(cst_t_opt(:,idx),[],1,'omitnan'))./sqrt(sum(~isnan(cst_t_opt(:,idx))));
    
%     indiv_S_diff(idx,:) = squeeze(cst_s(:,idx)-cst_s_opt(:,idx));
%     indiv_T_diff(idx,:) = squeeze(cst_t(:,idx)-cst_t_opt(:,idx));   
%     grandMn_S_diff(idx) = squeeze(mean(cst_s(:,idx)-cst_s_opt(:,idx),1,'omitnan'));
%     grandSE_S_diff(idx) = squeeze(std(cst_s(:,idx)-cst_s_opt(:,idx),[],1,'omitnan'))./sqrt(sum(~isnan(cst_s(:,idx)-cst_s_opt(:,idx))));
%     grandMn_T_diff(idx) = squeeze(mean(cst_t(:,idx)-cst_t_opt(:,idx),1,'omitnan'));
%     grandSE_T_diff(idx) = squeeze(std(cst_t(:,idx)-cst_t_opt(:,idx),[],1,'omitnan'))./sqrt(sum(~isnan(cst_t(:,idx)-cst_t_opt(:,idx))));
    
end

% Plot bar graph
fH = figure; clf; set(gcf,'Position',[1 1 800 941]); hold all;

b = bar([1:9],[grandMn_S;grandMn_T],'FaceColor','flat'); hold all;
b(1).CData = cmapBetas(1,:);
b(2).CData = cmapBetas(2,:);
eb1 = errorbar([1:9]-0.15,grandMn_S, grandSE_S, grandSE_S,'k','linewidth',0.5);
eb1.LineStyle = 'none';
eb1.CapSize = 0;
eb2 = errorbar([1:9]+0.15,grandMn_T, grandSE_T, grandSE_T,'k','linewidth',0.5);
eb2.LineStyle = 'none';
eb2.CapSize = 0;

% Plot individual subjects
for idx = roisNoIPS
    scatter(idx-0.15+(0.01*rand(1,size(cst_s_opt,1))),cst_s_opt(:,idx),30,[0.5 0.5 0.5], 'filled', ...
        'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerFaceAlpha',0.8,'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.8);
    scatter(idx+0.15+(0.01*rand(1,size(cst_t_opt,1))),cst_t_opt(:,idx),30,[0.5 0.5 0.5], 'filled', ...
        'MarkerFaceColor',[1 1 1], 'MarkerFaceAlpha',0.8,'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.8);
end
%% Add labels/legends
title('CST_{opt} Beta values for Sustained & Transient channels');
set(gca,'XTickLabel',string(roisToPlot)', 'XTickLabelRotation',45)
set(gca,'FontSize',12, 'TickDir','out')
box off;
ylabel('Beta weight (% signal change)')
ylim([-1.2 3.2])
legend('Sustained','Transient'); legend box off

if saveFigs
    fName = sprintf('SummaryGroupAverageSEM_BetaCSToptSustainedTransient');
    subDir = 'SupplFig8';
    if ~exist('saveFigDir','dir')
        saveFigDir = fullfile(simseqRootPath,'results','group');
    end
    thisSaveFigDir = fullfile(saveFigDir, subDir);
    if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
    saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
    print(gcf,fullfile(thisSaveFigDir,fName),'-depsc2','-painters','-r300','-loose')
end

%% Do some stats (2-Way repeated measured ANOVA)
% Test differences between CST_opt sustained vs CST_opt transient channel
s = cst_s(:,roisNoIPS)';
s_opt = cst_s_opt(1:size(s,2),roisNoIPS)';  % 7 subjects x 8 ROIs
t = cst_t(:,roisNoIPS)';
t_opt = cst_t_opt(1:size(t,2),roisNoIPS)'; % 7 subjects x 8 ROIs

subj_vec = repmat([1:size(s,2)],length(roisToLabel),1);
roi_vec = repmat(string(roisToLabel),[size(s,2),1]);
T1 = table(repmat(subj_vec(:),2,1),...
    repmat(roi_vec,2,1),...
    [ones(size(s,2)*length(roisToLabel),1); 2.*ones(size(s,2)*length(roisToLabel),1)],...
    [s_opt(:); t_opt(:)],...
    'VariableNames',{'Subject','ROI','ChannelType', 'Meas'});
rm = fitrm(T1,'Meas ~ ChannelType * ROI');
disp(rm.anova);

%% Test differences between CST_opt vs CST_fix vs S/T channels
subj_vec = repmat([1:size(s,2)],length(roisToLabel),1);
roi_vec = repmat(string(roisToLabel),[7,1]);
T2 = table(repmat(subj_vec(:),4,1),...Subject
    repmat(roi_vec,4,1),...ROI
    [ones(size(s,2)*length(roisToLabel),1); 2.*ones(size(s,2)*length(roisToLabel),1); ones(size(s,2)*length(roisToLabel),1); 2.*ones(size(s,2)*length(roisToLabel),1)],...ModelType
    [ones(size(s,2)*length(roisToLabel)*2,1); 2.*ones(size(s,2)*length(roisToLabel)*2,1)],...ChannelType
    [s(:); t(:);s_opt(:);t_opt(:)],...Measurement
    'VariableNames',{'Subject','ROI','ModelType','ChannelType', 'Meas'});
rm = fitrm(T2,'Meas ~ ModelType * ChannelType * ROI');
disp(rm.anova);