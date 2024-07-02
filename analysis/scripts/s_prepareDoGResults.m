%% s_prepareDoGResults

projectDir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/';

hemi          = {'both'};
prfModel      = 3; % Use linear pRF model to get center size
veThresh      = 0.1;
sizeRatio     = [7.4, 6.8, 7.3, 5.8, 4.2, 4.2, 4.7, 5, 3.3]; % V1,V2,V3,hV4,VO,V3AB,IPS,LO,TO
surroundFun   = @(sigma, roiIdx) sizeRatio(roiIdx).*sigma;
plotFigures   = true;
storeData     = true;
% gridFit = load(fullfile(fullfile(projectDir,'experiments/simseq/results/average/gridFit2_Results',...
%     'S1_S2_S3_S7_S8_S9_S10_S11_S12_S13_BlockAmps_pRFparams_modelFit_results_cv_variableBlockOnset_gridFit2_v4.mat')));

rois = {'V1','V2','V3','hV4'};

% Loop over subjects
subjnrs = [1,2,3,7:13];

spatialModel = 'differenceOfGaussiansFit';
temporalModel = '1ch-glm';
params.analysis.spatialModel = spatialModel;
params.analysis.temporalModel = temporalModel;
modelName = {'DoG'};
subSaveFolder = 'savedPredictions_fixedPRF_combT';
betaFolder    = 'betaR2_XValFit_OLS';

timePointsTrialToAverage = NaN(10,4,9,9);
timePointsTrialToAverageModel = timePointsTrialToAverage;
T_meanAmp = table();
T_meanAmpBslCorrected = table();
T_meanAmpModel = table();
T_meanAmpModelBslCorrected = table();
T_R2_Beta_pRFs = table();

%%
for subjnr = subjnrs
    s = (subjnr==subjnrs);
    
    % Define folders
    sesNr  = getSessionNrMainExp(subjnr);
    pths   = getSubjectPaths(projectDir, subjnr, sesNr);
    simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID, pths.session);
    cd(simSeqSessionDir)
    saveFolder = fullfile(pths.simseqResultsDir,'DogSimulation');
    saveFigureDir = fullfile(pths.figureDir,'mainExp','DoGsimulations');
    
    if ~exist(saveFolder,'dir'), mkdir(saveFolder); end
    loadDataFolder = fullfile(pths.simseqResultsDir, 'preprocData',...
        sprintf('varyStimDur_x_Area4_2x2design_4deg2_trimmedGaussian_bothHemi_mid_roiCorrected_fullrun2_v%d',pths.expversionNr));
    %%
    % Get stimulus file
    % fname = sprintf('stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun2_v%d.mat', pths.expversionNr);
    % stimfile         = fullfile(simSeqSessionDir, 'Stimuli', fname);
    % stim = logical(params.stim.images);
    
    % Load trial/condition data
    conditionData            = load(fullfile(fileparts(pths.simseqResultsDir),'conditionsForPlotting_mainExp.mat'));
    conditionData.loc        = getStimLocDeg(projectDir,subjnr,sesNr, true, false, hemi);
    conditionData.condColors = conditionData.trialData.condColors;
    
    for roiIdx = 1:4
        clear a TDoG TDoG_sim;
        roi    = sprintf('%s_stimcorner4_area4sq_eccen5',rois{roiIdx});

        % Load data
%         loadFile = fullfile(pths.simseqResultsDir,...
%             subSaveFolder, [betaFolder '_gridFit'], ...
%             sprintf('betaValsFull_%s_cssFit_1ch-glm_run12_offsetFlag1_OLS_bslCorrected.mat', ...
%             roi));
        loadFile = fullfile(pths.simseqResultsDir,...
            subSaveFolder, [betaFolder '_gridFit2'], ...
            sprintf('betaValsFull_%s_onegaussianFit_3ch-stLN_run12_offsetFlag1_OLS_bslCorrected_expBestR2.mat', ...
            roi));
        load(loadFile); a = optimal;
        amplData = a.T.DataTSCV;
        allAmpData{s,roiIdx} = amplData;
        
        % Load DoG Simulated predictions
        b = load(fullfile(saveFolder, ...
            sprintf('dogSimulation_%s_%s.mat',pths.subjID, roi)));
        
        % Clear some memory
        b.params.stim.images = [];
        b.params.stim.images_unconvolved = [];
        params = b.params;
        params.analysis.spatial.varexpl = a.data.params.varexpl.both;
        
        amplModel = b.TDoG.TS;        
        allAmpModel{s,roiIdx} = amplModel;
        timeWindow = a.er_data.params.timeWindow{1};
        twSelectionType = 'variableOnset';
        numConditions = 8;
        
        %%
        dataRuns{s,roiIdx}       = b.Y_crossval_mean;
        modelRuns{s,roiIdx}      = b.sumChannelPredictionCrossvalMean;
        betaRuns{s,roiIdx}       = squeeze(mean(b.B_crossval,3,'omitnan'))';
        r2Runs{s,roiIdx}         = b.R2_crossval_mean;
        
        if length(a.data.splitHalfRel.both.runCV) ~= length(b.R2_crossval_mean)
            bothVoxIdx = [a.alignMe.dataRuns.selected.lh; a.alignMe.dataRuns.selected.rh];
            assert(length(bothVoxIdx)==length(b.R2_crossval_mean))
            noiseceilingRun{s,roiIdx} = a.data.splitHalfRel.both.runCV(bothVoxIdx);
        else
            noiseceilingRun{s,roiIdx} = a.data.splitHalfRel.both.runCV;
        end
        % params.x0.both            = params.x0.both;
        % params.y0.both            = params.y0.both;
        % params.effectiveSize.both = params.effectiveSize.both;
        % params.varexpl.both       = params.varexpl.both;
        % params.exp_spatial.both   = params.exponent.both;
        % if strcmp(temporalModel, '1ch-dcts')
        %     params.temporal.tau1      = a.data.params.temporal.param.tau1;
        %     params.temporal.weight    = a.data.params.temporal.param.weight;
        %     params.temporal.tau2      = a.data.params.temporal.param.tau2;
        %     params.temporal.n         = a.data.params.temporal.param.n;
        %     params.temporal.sigma     = a.data.params.temporal.param.sigma;
        % elseif strcmp(temporalModel, '3ch-stLN')
        %     params.exp_temporal.both  = a.data.params.exponent_temporal.both;
        %     params.temporal.exponent   = a.data.params.temporal.param.exponent;
        %     params.temporal.tau_s      = a.data.params.temporal.param.tau_s;
        %     params.temporal.tau_t      = a.data.params.temporal.param.tau_t;
        % end
        paramsROI{s,roiIdx}   = params;
        alignMeROI{s,roiIdx}  = a.alignMe;
        T_R2_Beta_pRFs = [T_R2_Beta_pRFs; table(subjnr, rois(roiIdx), modelName, betaRuns(s,roiIdx), ...
            r2Runs(s,roiIdx),{NaN},params, ...
            'VariableNames',{'Subject','ROI','modelType','betaRun','r2Run','r2Trial','pRFsParams'})];
        
        for cc = 1:4
            % Get variable time points based on baseline corrected mean across voxels.
            dataCondition = amplData{cc} - mean(amplData{cc}((timeWindow<=0),:),1);
            meanDataCondition = nanmean(dataCondition,2)';
            timePointsTrialToAverage(s,cc,roiIdx,:) = getTimePointsTrialToAverage(...
                twSelectionType, timeWindow, meanDataCondition);
            
            amplTrialBslCorrected{s,cc,roiIdx} = nanmean(dataCondition(squeeze(timePointsTrialToAverage(s,cc,roiIdx,:)),:),1);
            amplTrial{s,cc,roiIdx} = nanmean(amplData{cc}(squeeze(timePointsTrialToAverage(s,cc,roiIdx,:)),:),1);
            
            dataCondition = amplData{cc+(numConditions/2)} - mean(amplData{cc+(numConditions/2)}((timeWindow<=0),:),1);
            meanDataCondition = nanmean(dataCondition,2)';
            timePointsTrialToAverage(s,cc+(numConditions/2),roiIdx,:) = getTimePointsTrialToAverage(...
                twSelectionType, timeWindow, meanDataCondition);
            amplTrialBslCorrected{s,cc+(numConditions/2),roiIdx} = nanmean(dataCondition(squeeze(timePointsTrialToAverage(s,cc+(numConditions/2),roiIdx,:)),:),1);
            amplTrial{s,cc+(numConditions/2),roiIdx} = nanmean(amplData{cc+(numConditions/2)}(squeeze(timePointsTrialToAverage(s,cc+(numConditions/2),roiIdx,:)),:),1);
            
            % Create mean amplitude table
            T_meanAmp = [T_meanAmp; table(subjnr, cc, rois(roiIdx),...
                amplTrial(s,cc,roiIdx), ...
                amplTrial(s,cc+(numConditions/2),roiIdx), ...
                'VariableNames',{'Subject','Condition','ROI','SEQ_MeanAmplitude','SIM_MeanAmplitude'})];
            T_meanAmpBslCorrected = [T_meanAmpBslCorrected; table(subjnr, cc, rois(roiIdx), ...
                amplTrialBslCorrected(s,cc,roiIdx), ...
                amplTrialBslCorrected(s,cc+(numConditions/2),roiIdx), ...
                'VariableNames',{'Subject','Condition','ROI','SEQ_MeanAmplitude_bslcorr','SIM_MeanAmplitude_bslcorr'})];
            
        end
        
        %%
        for cc = 1:4
            % Get variable time points based on baseline corrected mean across voxels.
            modelCondition = amplModel{cc} - mean(amplModel{cc}((timeWindow<=0),:),1);
            meanModelCondition = nanmean(modelCondition,2)';
            timePointsTrialToAverageModel(s,cc,roiIdx,:) = getTimePointsTrialToAverage(...
                twSelectionType, timeWindow, meanModelCondition);
            modelTrialBslCorrected{s,cc,roiIdx} = nanmean(modelCondition(squeeze(timePointsTrialToAverageModel(s,cc,roiIdx,:)),:),1);
            modelTrial{s,cc,roiIdx} = nanmean(amplModel{cc}(squeeze(timePointsTrialToAverageModel(s,cc,roiIdx,:)),:),1);
            
            modelCondition = amplModel{cc+(numConditions/2)} - mean(amplModel{cc+(numConditions/2)}((timeWindow<=0),:),1);
            meanModelCondition = nanmean(modelCondition,2)';
            timePointsTrialToAverageModel(s,cc+(numConditions/2),roiIdx,:) = getTimePointsTrialToAverage(...
                twSelectionType, timeWindow, meanModelCondition);
            modelTrialBslCorrected{s,cc+(numConditions/2),roiIdx} = nanmean(modelCondition(squeeze(timePointsTrialToAverageModel(s,cc+(numConditions/2),roiIdx,:)),:),1);
            modelTrial{s,cc+(numConditions/2),roiIdx} = nanmean(amplModel{cc+(numConditions/2)}(squeeze(timePointsTrialToAverageModel(s,cc+(numConditions/2),roiIdx,:)),:),1);
            
            
            
            T_meanAmpModel = [T_meanAmpModel; ...
                table(subjnr, cc, rois(roiIdx), modelName, ...
                modelTrial(s,cc,roiIdx), ...
                modelTrial(s,cc+(numConditions/2),roiIdx), ...
                'VariableNames',{'Subject','Condition','ROI','modelType','SEQ_ModelAmplitude','SIM_ModelAmplitude'})];
            T_meanAmpModelBslCorrected = [T_meanAmpModelBslCorrected; ...
                table(subjnr, cc, rois(roiIdx),modelName, ...
                modelTrialBslCorrected(s,cc,roiIdx), ...
                modelTrial(s,cc+(numConditions/2),roiIdx), ...
                'VariableNames',{'Subject','Condition','ROI','modelType','SEQ_ModelAmplitude_bslcorr','SIM_ModelAmplitude_bslcorr'})];
            
        end
    end
end

if storeData
    % Save summary data we so don't have to compute it next time
    sNames = sprintf('S%i_',subjnrs);
    allROIs = pths.allROIs;
    save(fullfile(projectDir,'experiments/simseq/results/average/DoG_Results',...
        sprintf('%sBlockAmps_pRFparams_modelFit_results_cv_variableBlockOnset_DoG_1ch-glm_v3.mat',sNames)),...
        'T_meanAmp','T_meanAmpModel','T_R2_Beta_pRFs','paramsROI','noiseceilingRun',...
        'allAmpData','allAmpModel','amplTrial','modelTrial','alignMeROI',...
        'betaRuns','r2Runs','conditionData',...
        'subjnrs','allROIs','timePointsTrialToAverage','timePointsTrialToAverageModel',...
        'spatialModel','temporalModel','timePointsTrialToAverage','-v7.3');
end

% 
% %% Panel A: Plot V1 voxel time series per condition
% plotDataModelFit_singleVoxel_StimulusConditionBlocks(projectDir,subjnr,...
%     'roisToPlot',{'V1'},'selectedDataVoxels',941, 'plotModelFlag',false)
% 
% %% Panel B: Scatter plot V1 all voxel data, subject S3
% plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, 'subjnrs',subjnr,'saveFigs',false, 'roisToPlot',1);

