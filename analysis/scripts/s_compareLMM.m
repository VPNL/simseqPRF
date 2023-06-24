%% s_compareLMMs.m

% Script to compare goodness of fit for different linear mixed models
% Main LMM: fixed ampl & condition effects + random subject intercept and slope for each condition.
% LMM alternative 0: fixed ampl + 1 random subject intercept (same for each condition).
% LMM alternative 1: fixed ampl & condition effects + 1 random subject intercept (same for each condition).
% LMM alternative 2: fixed ampl & condition effects + 4 random subject intercept (for each condition).

%% Get dataset
projectDir = simseqRootPath;
ds = createDataTable(projectDir);

% Reorder ROIs
roisToPlot  = unique(ds.ROI,'stable');
newRoiOrder = [1,2,3,4,5,6,7,8,9];
roisToPlot  = roisToPlot(newRoiOrder);

%% MAIN LMM
fLMM = cell(1,length(roisToPlot));
LMMlabel = 'Main';
fitStr = 'MeanSimAmp ~ MeanSeqAmp * Condition + (1 + MeanSeqAmp * Condition | Subject)';

nrConditions         = length(unique(ds.Condition));
nrSubjects           = length(unique(ds.Subject));
nrFixedEffectLevels  = [nrConditions,nrConditions];
nrRandomEffectLevels = [nrConditions,nrConditions];

fixedIntercepts   = fLMM;
fixedIntercepts_CI = fLMM;
fixedIntercepts_SE = fLMM;
fixedSlopes       = fLMM;
fixedSlopes_CI    = fLMM;
fixedSlopes_SE    = fLMM;
subjIntercepts    = fLMM;
subjIntercepts_CI = fLMM;
subjIntercepts_SE = fLMM;
subjSlopes        = fLMM;
subjSlopes_CI     = fLMM;
subjSlopes_SE     = fLMM;

for ii = 1:length(roisToPlot)
    [fLMM{ii},fixedIntercepts{ii},fixedIntercepts_CI{ii}, fixedIntercepts_SE{ii},...
        fixedSlopes{ii}, fixedSlopes_CI{ii}, fixedSlopes_SE{ii}, ...
        subjIntercepts{ii}, subjIntercepts_CI{ii}, subjIntercepts_SE{ii},...
        subjSlopes{ii}, subjSlopes_CI{ii}, subjSlopes_SE{ii}] = ...
        getLMMcoefficients(ds(ds.ROI==roisToPlot(ii),:),fitStr, ...
        nrFixedEffectLevels,nrRandomEffectLevels);
end

%% Fit LLM alternatives:

%% LMM_alt0
% Fixed ampl & condition effects + random subject intercept for each condition
fitStr_alt0 = 'MeanSimAmp ~ MeanSeqAmp + (1 | Subject)';
nrFixedEffectLevels  = 1;
nrRandomEffectLevels = 1;

% Preallocate space
fLMM_alt0                = cell(1,length(roisToPlot));
fixedIntercepts_alt0     = fLMM_alt0; % fixed effect intercepts (same for all subjects)
fixedIntercepts_CI_alt0  = fLMM_alt0; % fixed effect intercepts 95% confidence interval
fixedIntercepts_SE_alt0  = fLMM_alt0; % fixed effect intercepts standard error
fixedSlopes_alt0         = fLMM_alt0; % fixed effect slopes (same for all subjects)
fixedSlopes_CI_alt0      = fLMM_alt0; % fixed effect slopes 95% confidence interval
fixedSlopes_SE_alt0      = fLMM_alt0; % fixed effect slopes standard error
subjIntercepts_alt0      = fLMM_alt0; % Random subject intercept
subjIntercepts_CI_alt0   = fLMM_alt0; % Random subject intercept 95% confidence interval
subjIntercepts_SE_alt0   = fLMM_alt0; % Random subject intercept standard error
subjSlopes_alt0          = fLMM_alt0; % should be empty
subjSlopes_CI_alt0       = fLMM_alt0; % should be empty
subjSlopes_SE_alt0       = fLMM_alt0; % should be empty

for ii = 1:length(roisToPlot)
    [fLMM_alt0{ii},fixedIntercepts_alt0{ii},fixedIntercepts_CI_alt0{ii}, fixedIntercepts_SE_alt0{ii},...
        fixedSlopes_alt0{ii}, fixedSlopes_CI_alt0{ii}, fixedSlopes_SE_alt0{ii}, ...
        subjIntercepts_alt0{ii}, subjIntercepts_CI_alt0{ii}, subjIntercepts_SE_alt0{ii},...
        subjSlopes_alt0{ii}, subjSlopes_CI_alt0{ii}, subjSlopes_SE_alt0{ii}] = ...
        getLMMcoefficients(ds(ds.ROI==roisToPlot(ii),:),fitStr_alt0, ...
        nrFixedEffectLevels,nrRandomEffectLevels);
end

for ii = 1:length(roisToPlot)
    fprintf('LMM0 v Main: %s\n', string(roisToPlot(ii)))
    compareLMMs_LMM0_v_Main{ii} = compareLMMmodels(fLMM_alt0{ii},fLMM{ii});
    fprintf('\n')
end

%% Visualize LMM0
x = [-1:0.01:2];
idx = 1;
for sj = 1:nrSubjects
figure(10+sj); set(gcf, 'Position', [8,454,1906,523]); clf;
    for c = 1:nrConditions
    
        xToPlot = ds.MeanSeqAmp(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(c));
        
        yToPlot = ds.MeanSimAmp(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(c));
        
        tmpT = table();
        tmpT.MeanSeqAmp = xToPlot;
        tmpT.MeanSimAmp = yToPlot;
        tmpT.Subject    = ones(size(xToPlot)).*sj;
        tmpT.ROI        = repmat(roisToPlot(idx),size(xToPlot));
        tmpT.Condition  = ones(size(xToPlot)).*c;
        
        tmpT           = table2dataset(tmpT);
        tmpT.Condition = nominal(tmpT.Condition);
        tmpT.ROI       = nominal(tmpT.ROI);
        
        [ypred, ypredCI] = predict(fLMM{idx}, tmpT);
        [ypred0, ypredCI0] = predict(fLMM_alt0{idx}, tmpT); % 'Conditional',0);
        
        subplot(1,4,c); hold all;
        title(sprintf('S%d: %s Condition %d', sj, string(roisToPlot(idx)),c));
        plot(xToPlot, yToPlot,'o')
        plot(xToPlot, ypred,'r-', 'LineWidth',4)
        plot(xToPlot, ypred0,'b-', 'LineWidth',4)
        
        plot([0 0], [-1 2.5],'k-')
        plot([-1 2.5],[0 0 ],'k-')
        x = min(xToPlot):0.01:max(xToPlot);
        plot(0,subjIntercepts{idx}(sj,c),'go', 'LineWidth',4)
        plot(x,(x.*subjSlopes{idx}(sj,c))+subjIntercepts{idx}(sj,c), 'g:', 'LineWidth',4)
        plot(0,subjIntercepts_alt0{idx}(sj),'co', 'LineWidth',4)
        plot(x,(x.*fixedSlopes_alt0{idx})+subjIntercepts_alt0{idx}(sj), 'c:', 'LineWidth',4)

        axis square
    end
    sgtitle('LMM main (red/green) v. LMM alt0 (blue/cyan) - fixed effect: SEQ Ampl, random effect: Subject intercept')
end

%% LLM alt1
% fixed ampl & condition effects + 1 random subject intercept (same for each condition).

fitStr = 'MeanSimAmp ~ MeanSeqAmp * Condition + (1 | Subject)';

nrConditions    = length(unique(ds.Condition));
nrSubjects      = length(unique(ds.Subject));
nrFixedEffectLevels  = [nrConditions,nrConditions]; % 4 intercepts, 4 slopes for each condition
nrRandomEffectLevels = 1; % 1 intercept per subject

% Preallocate space
fLMM_alt1                   = cell(1,length(roisToPlot));
fixedIntercepts_alt1     = fLMM_alt1; % fixed effect intercepts (same for all subjects)
fixedIntercepts_CI_alt1  = fLMM_alt1; % fixed effect intercepts 95% confidence interval
fixedIntercepts_SE_alt1  = fLMM_alt1; % fixed effect intercepts standard error
fixedSlopes_alt1         = fLMM_alt1; % fixed effect slopes (same for all subjects)
fixedSlopes_CI_alt1      = fLMM_alt1; % fixed effect slopes 95% confidence interval
fixedSlopes_SE_alt1      = fLMM_alt1; % fixed effect slopes tandard error
subjIntercepts_alt1      = fLMM_alt1; % Random subject intercept
subjIntercepts_CI_alt1   = fLMM_alt1; % Random subject intercept 95% confidence interval
subjIntercepts_SE_alt1   = fLMM_alt1; % Random subject intercept standard error
subjSlopes_alt1          = fLMM_alt1; % should be empty
subjSlopes_CI_alt1       = fLMM_alt1; % should be empty
subjSlopes_SE_alt1       = fLMM_alt1; % should be empty

for ii = 1:length(roisToPlot)
    
    [fLMM_alt1{ii}, ...
        fixedIntercepts_alt1{ii},fixedIntercepts_CI_alt1{ii}, fixedIntercepts_SE_alt1{ii},... 
        fixedSlopes_alt1{ii}, fixedSlopes_CI_alt1{ii}, fixedSlopes_SE_alt1{ii}, ...
        subjIntercepts_alt1{ii}, subjIntercepts_CI_alt1{ii}, subjIntercepts_SE_alt1{ii}, ...
        subjSlopes_alt1{ii}, subjSlopes_CI_alt1{ii}, subjSlopes_SE_alt1{ii}] = ...
        getLMMcoefficients(ds(ds.ROI==roisToPlot(ii),:),fitStr, ...
        nrFixedEffectLevels,nrRandomEffectLevels);
end

for ii = 1:length(roisToPlot)
    fprintf('\nLMM1 v Main: %s\n', string(roisToPlot(ii)))
    compareLMMs_LMM1_v_Main{ii} = compareLMMmodels(fLMM_alt1{ii},fLMM{ii});
end

%% Visualize LMM alt1: fixed and random effects for V1

idx = 1;
for sj = 1:nrSubjects
    figure(10+sj); set(gcf, 'Position', [8,454,1906,523]); clf;
    for c = 1:nrConditions
        xToPlot = ds.MeanSeqAmp(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(c));
        
        yToPlot = ds.MeanSimAmp(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(c));
        
        tmpT = table();
        tmpT.MeanSeqAmp = xToPlot;
        tmpT.MeanSimAmp = yToPlot;
        tmpT.Subject    = ones(size(xToPlot)).*sj;
        tmpT.ROI        = repmat(roisToPlot(idx),size(xToPlot));
        tmpT.Condition  = ones(size(xToPlot)).*c;
        
        tmpT           = table2dataset(tmpT);
        tmpT.Condition = nominal(tmpT.Condition);
        tmpT.ROI       = nominal(tmpT.ROI);
        
        [ypred, ypredCI] = predict(fLMM{idx}, tmpT);
        [ypred1, ypredCI1] = predict(fLMM_alt1{idx}, tmpT); %, 'Conditional',0);
        
        subplot(1,4,c); hold all;
        title(sprintf('S%d: %s Condition %d', sj, string(roisToPlot(idx)),c));
        plot(xToPlot, yToPlot,'o')
        plot(xToPlot, ypred,'r-', 'LineWidth',4)
        plot(xToPlot, ypred1,'b-', 'LineWidth',4)
        
        plot([0 0], [-1 2.5],'k-')
        plot([-1 2.5],[0 0 ],'k-')
        x = min(xToPlot):0.01:max(xToPlot);
        plot(0,subjIntercepts{idx}(sj,c),'go', 'LineWidth',4)
        plot(x,(x.*subjSlopes{idx}(sj,c))+subjIntercepts{idx}(sj,c), 'g:', 'LineWidth',4)
        plot(0,subjIntercepts_alt1{idx}(sj,c),'co', 'LineWidth',4)
        plot(x,(x.*fixedSlopes_alt1{idx}(c))+subjIntercepts_alt1{idx}(sj,c), 'c:', 'LineWidth',4)
        axis square
    end
    sgtitle('LMM main v. LMM alternative 1 - fixed effect: SEQ Ampl * Condition, random effect: Subject intercept')
end



%% LMM_alt2
% Fixed ampl & condition effects + random subject intercept for each condition
fitStr_alt2 = 'MeanSimAmp ~ MeanSeqAmp * Condition + (Condition | Subject)';
nrFixedEffectLevels = [nrConditions,nrConditions];
nrRandomEffectLevels = nrConditions;

% Preallocate space
fLMM_alt2                = cell(1,length(roisToPlot));
fixedIntercepts_alt2     = fLMM_alt2; % fixed effect intercepts (same for all subjects)
fixedIntercepts_CI_alt2  = fLMM_alt2; % fixed effect intercepts 95% confidence interval
fixedIntercepts_SE_alt2  = fLMM_alt2; % fixed effect intercepts standard error
fixedSlopes_alt2         = fLMM_alt2; % fixed effect slopes (same for all subjects)
fixedSlopes_CI_alt2      = fLMM_alt2; % fixed effect slopes 95% confidence interval
fixedSlopes_SE_alt2      = fLMM_alt2; % fixed effect slopes standard error
subjIntercepts_alt2      = fLMM_alt2; % Random subject intercept (10x4)
subjIntercepts_CI_alt2   = fLMM_alt2; % Random subject intercept 95% confidence interval (10x4x2)
subjIntercepts_SE_alt2   = fLMM_alt2; % Random subject intercept standard error (10x4)
subjSlopes_alt2          = fLMM_alt2; % should be empty
subjSlopes_CI_alt2       = fLMM_alt2; % should be empty
subjSlopes_SE_alt2       = fLMM_alt2; % should be empty

for ii = 1:length(roisToPlot)
    
    [fLMM_alt2{ii}, ...
        fixedIntercepts_alt2{ii},fixedIntercepts_CI_alt2{ii}, fixedIntercepts_SE_alt2{ii}, ...
        fixedSlopes_alt2{ii}, fixedSlopes_CI_alt2{ii}, fixedSlopes_SE_alt2{ii},...
        subjIntercepts_alt2{ii}, subjIntercepts_CI_alt2{ii}, subjIntercepts_SE_alt2{ii},...
        subjSlopes_alt2{ii}, subjSlopes_CI_alt2{ii}, subjSlopes_SE_alt2{ii}] = ...
        getLMMcoefficients(ds(ds.ROI==roisToPlot(ii),:),fitStr_alt2, ...
        nrFixedEffectLevels,nrRandomEffectLevels);
end

for ii = 1:length(roisToPlot)
    fprintf('\nLMM2 v Main: %s\n', string(roisToPlot(ii)))
    compareLMMs_LMM2_v_Main{ii} = compareLMMmodels(fLMM_alt2{ii},fLMM{ii});
end

%% Visualize LMM2
idx = 1;
for sj = 1:nrSubjects
    figure(10+sj); set(gcf, 'Position', [8,454,1906,523]); clf;
    for c = 1:nrConditions
        xToPlot = ds.MeanSeqAmp(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(c));
        
        yToPlot = ds.MeanSimAmp(ds.Subject==sj & ...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(c));
        
        tmpT = table();
        tmpT.MeanSeqAmp = xToPlot;
        tmpT.MeanSimAmp = yToPlot;
        tmpT.Subject    = ones(size(xToPlot)).*sj;
        tmpT.ROI        = repmat(roisToPlot(idx),size(xToPlot));
        tmpT.Condition  = ones(size(xToPlot)).*c;
        
        tmpT           = table2dataset(tmpT);
        tmpT.Condition = nominal(tmpT.Condition);
        tmpT.ROI       = nominal(tmpT.ROI);
        
        [ypred, ypredCI] = predict(fLMM{idx}, tmpT);
        [ypred2, ypredCI2] = predict(fLMM_alt2{idx}, tmpT);%, 'Conditional',0);
        
        subplot(1,4,c); hold all;
        title(sprintf('S%d: %s Condition %d', sj, string(roisToPlot(idx)),c));
        plot(xToPlot, yToPlot,'ko')
        plot(xToPlot, ypred,'r-', 'LineWidth',4)
        plot(xToPlot, ypred2,'b-', 'LineWidth',4)
        
        plot([0 0], [-1 2.5],'k-')
        plot([-1 2.5],[0 0 ],'k-')
        
        x = min(xToPlot):0.01:max(xToPlot);
        plot(0,subjIntercepts{idx}(sj,c),'go', 'LineWidth',4)
        plot(x,(x.*subjSlopes{idx}(sj,c))+subjIntercepts{idx}(sj,c), 'g:', 'LineWidth',4)
        plot(0,subjIntercepts_alt2{idx}(sj,c),'co', 'LineWidth',4)
        plot(x,(x.*fixedSlopes_alt2{idx}(c))+subjIntercepts_alt2{idx}(sj,c), 'c:', 'LineWidth',4)
        axis square
    end
end


%% Compare alternative LMMs to eachother
for ii = 1:length(roisToPlot)
    fprintf('\nLMM0 v LMM1: %s\n', string(roisToPlot(ii)))
    compareLMMs_LMM0_vs_LMM1{ii} = compareLMMmodels(fLMM_alt0{ii}, fLMM_alt1{ii});
end
for ii = 1:length(roisToPlot)
    fprintf('\nLMM1 v LMM2: %s\n', string(roisToPlot(ii)))
    compareLMMs_LMM1_vs_LMM2{ii} = compareLMMmodels(fLMM_alt1{ii}, fLMM_alt2{ii});
end

%% Compute difference in model goodness of fit (AIC, BIC, log-likelihood)
for ii = 1:length(roisToPlot)
    AIC_0_v_Main(ii,:) = compareLMMs_LMM0_v_Main{ii}.AIC;
    AIC_1_v_Main(ii,:) = compareLMMs_LMM1_v_Main{ii}.AIC;
    AIC_2_v_Main(ii,:) = compareLMMs_LMM2_v_Main{ii}.AIC;

    BIC_0_v_Main(ii,:) = compareLMMs_LMM0_v_Main{ii}.BIC;
    BIC_1_v_Main(ii,:) = compareLMMs_LMM1_v_Main{ii}.BIC;
    BIC_2_v_Main(ii,:) = compareLMMs_LMM2_v_Main{ii}.BIC;
    
    LogLik_0_v_Main(ii,:) = compareLMMs_LMM0_v_Main{ii}.LogLik;
    LogLik_1_v_Main(ii,:) = compareLMMs_LMM1_v_Main{ii}.LogLik;
    LogLik_2_v_Main(ii,:) = compareLMMs_LMM2_v_Main{ii}.LogLik;
    
    pVals_0_v_Main(ii,:) = compareLMMs_LMM0_v_Main{ii}.pValue(1);
    pVals_1_v_Main(ii,:) = compareLMMs_LMM1_v_Main{ii}.pValue(1);
    pVals_2_v_Main(ii,:) = compareLMMs_LMM2_v_Main{ii}.pValue(1);
 for ii = 1:length(roisToPlot)   
    DF_0_v_Main(ii,:) = compareLMMs_LMM0_v_Main{ii}.deltaDF;
    DF_1_v_Main(ii,:) = compareLMMs_LMM1_v_Main{ii}.deltaDF;
    DF_2_v_Main(ii,:) = compareLMMs_LMM2_v_Main{ii}.deltaDF;
end
%%
AIC_main_v_m012 = [AIC_0_v_Main(:,1)-AIC_0_v_Main(:,2), AIC_1_v_Main(:,1)-AIC_1_v_Main(:,2), AIC_2_v_Main(:,1)-AIC_2_v_Main(:,2)];
BIC_main_v_m012 = [BIC_0_v_Main(:,1)-BIC_0_v_Main(:,2), BIC_1_v_Main(:,1)-BIC_1_v_Main(:,2), BIC_2_v_Main(:,1)-BIC_2_v_Main(:,2)];
LL_main_v_m012  = [LogLik_0_v_Main(:,1)-LogLik_0_v_Main(:,2), LogLik_1_v_Main(:,1)-LogLik_1_v_Main(:,2), LogLik_2_v_Main(:,1)-LogLik_2_v_Main(:,2)];

%% Plot AIC, BIC, LL differences

figure(99); clf; set(gcf,'Position',[1602,363,1279,834]);
sgtitle('LMM comparisions relative to model 3 (random subject slope & intercept)')
subplot(131); 
b1 = bar(1:length(roisToPlot),LL_main_v_m012);
b1(1).FaceColor = [0 0 0];
b1(2).FaceColor = [1 1 1];
b1(1).EdgeColor = [0 0 0];
b1(2).EdgeColor = [0 0 0];
box off;
set(gca,'XTickLabel', string(roisToPlot),'XTickLabelRotation',45);
ylabel('Log Likelihood relative to main LMM')
legend({'LMM 0 ', 'LMM 1', ...
        'LMM 2'}, 'Location', 'SouthEast')
legend boxoff
ylim([min(LL_main_v_m012(:))-1000 0])

subplot(132);
b2 = bar(1:length(roisToPlot),AIC_main_v_m012);
b2(1).FaceColor = [0 0 0];
b2(2).FaceColor = [1 1 1];
b2(1).EdgeColor = [0 0 0];
b2(2).EdgeColor = [0 0 0];
box off;
set(gca,'XTickLabel', string(roisToPlot),'XTickLabelRotation',45);
ylabel('AIC relative to main LMM')
legend({'LMM 0',...
        'LMM 1', ...
        'LMM 2'}, 'Location', 'NorthEast')
legend boxoff

subplot(133);
b3 = bar(1:length(roisToPlot),BIC_main_v_m012);
b3(1).FaceColor = [0 0 0];
b3(2).FaceColor = [1 1 1];
b3(1).EdgeColor = [0 0 0];
b3(2).EdgeColor = [0 0 0];
box off;
set(gca,'XTickLabel', string(roisToPlot),'XTickLabelRotation',45);
ylabel('BIC relative to main LMM')
legend({'LMM 0','LMM 1',...
        'LMM 2'}, 'Location', 'NorthEast')
legend boxoff

fName = 'LMM_LL_AIC_BIC_Modelcomparision_v2';
pths  = getSubjectPaths(projectDir,1,3);

saveFigDir = fullfile(fileparts(fileparts(pths.figureDir)), 'average', 'figures', ...
    'main_fixedPRF_combT_variableBlockOnset_gridFit2', 'otherLMMs');
subDir     = 'lmmFit_randomIntcptSlopeSubjInteraction_seq_vs_sim_ampl_gridFit2_v2';
saveFigDirGroup = fullfile(saveFigDir, subDir,'group');
if ~exist(fullfile(saveFigDirGroup),'dir'), mkdir(fullfile(saveFigDirGroup)); end
saveas(gcf, fullfile(saveFigDirGroup, [fName '.png']))
print(gcf,fullfile(saveFigDirGroup,fName),'-depsc')