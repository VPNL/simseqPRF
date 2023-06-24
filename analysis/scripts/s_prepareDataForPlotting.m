%% s_prepareDataForPlotting.m
% Script to prepare data for "createDataTable" function. Mainly combines
% ROIs (like VO1 and VO2 into VO1/2) and single voxel pRF parameters with 
% single voxel model fitting results.

%% Define params
subjnrs        = [1,2,3,7,8,9,10,11,12,13]; %7; %7
projectDir     = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/';
pilotNr        = 3;
spatialModels  = {'onegaussianFit','cssFit','onegaussianFit'};  % Choose from'cssFit' or 'onegaussianFit'
temporalModels = {'1ch-glm','1ch-glm','3ch-stLN'};
regressionType = '_OLS';
bslpostfix     = '_bslCorrected';
twSelectionType = 'variableOnset';
modelNames      = {'LSS','CSS','CST'};

recomputeData  = true; % load stored (false) or recompute values (true)
storeData      = true;  % store data when recomputing values;
offsetFlag     = true;  % Use betas from modelfit with added row of ones to allow for offset in GLM
bslCorrectFlag = true;  % Use data with baseline correction
%
%
% combROIs     = [5,8,10,12]; %VO1/2 (5,6), LO1/2(8,9), TO1/2 (10,11), IPS0/2 (12,13)
% labelsComb   = {'V1','V2','V3','hV4','VO1/2','LO1/2','TO1/2','V3AB','IPS0/1'};

% subSaveFolder = 'savedPredictions_fixedPRF_combT';
% betaFolder     = 'betaR2_XValFit_OLS';
% fitName        = '_gridFit2';
% stimFileName   = 'stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun';
% subLoadFolder  = 'varyStimDur_x_Area4_2x2design_4deg2_trimmedGaussian_bothHemi_mid_roiCorrected_fullrun2';
% postFixTtle    = [];
% roiType        = 'stimcorner4_area4sq_eccen5'; %     % Use stimulus corner, can also be one square (e.g.: 'square1')
% stimRun        = 1;                  % there are two unique runs, we pick the first
% hemi           = 'both';
% numConditions  = 8;
% nrTRsRemovedAtStart = 8;

%%
% if recomputeData
T_meanAmp = table();
T_meanAmpBslCorrected = table();
T_meanAmpModel = table();
T_meanAmpModelBslCorrected = table();
T_R2_Beta_pRFs = table();

for subjnr = subjnrs
    
    % Get data session nr
    s = (subjnr==subjnrs);
    pths = getSubjectPaths(projectDir,subjnr,sesNr);
    fname = sprintf('%s2_v%d.mat', ...
        stimFileName, pths.expversionNr);
    
    % Define subfolder and ROIs
    %         subLoadFolder2 = sprintf('%s_v%d',subLoadFolder,pths.expversionNr);
    rois       = pths.definedROIsBOTH;
    
    % Define session dir
    simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID);
    cd(simSeqSessionDir);
    
    % Load trial/condition data
    conditionData.loc        = getStimLocDeg(projectDir,subjnr,sesNr, true, false, hemi);
    
    for idx = 1:length(rois)
        
        roiIdx = find(strcmp(pths.allROIs,rois{idx}));
        
        for mm = 1:length(spatialModels)
            % Folder with prediction and data
            spatialModel = spatialModels{mm};
            temporalModel = temporalModels{mm};
            
            % Load data and model prediction
            loadFile = fullfile(pths.simseqResultsDir,...
                sprintf('%s_%s_%s_cvfitsResults.mat', ...
                pths.subjID, rois{idx},spatialModel, temporalModel));
            
            a = load(loadFile);
            
            a.params.stim.images = [];
            a.params.stim.images_unconvolved = [];
            dataRuns{s,roiIdx}(mm,:,:)  = a.Y_crossval_mean;
            modelRuns{s,roiIdx}(mm,:,:) = a.sumChannelPredictionCrossvalMean;
            betaRuns{s,roiIdx,mm}       = a.B_crossval_mean(:,[1:end-1])';
            r2Runs{s,roiIdx,mm}         = a.R2_crossval_mean;
            noiseceilingRun{s,roiIdx,mm} = a.data.splitHalfRel.both.runCV;
            params.x0.both            = a.data.params.x0.both;
            params.y0.both            = a.data.params.y0.both;
            params.effectiveSize.both = a.data.params.effectiveSize.both;
            params.varexpl.both       = a.data.params.varexpl.both;
            params.exp_spatial.both   = a.data.params.exponent.both;
            if mm==3
                params.exp_temporal.both  = a.data.params.exponent_temporal.both;
            end
            paramsROI{s,roiIdx,mm}       = params;
            T_R2_Beta_pRFs = [T_R2_Beta_pRFs; table(subjnr, rois(idx), modelNames(mm), betaRuns(s,roiIdx,mm), ...
                r2Runs(s,roiIdx,mm), r2Trials(s,roiIdx,mm), {NaN}, ...
                'VariableNames',{'Subject','ROI','modelType','betaRun','r2Run','r2Trial','pRFsParams'})];
            
            amplModel       = a.T.ModelTSCVMn;
            allAmpModel{s,roiIdx,mm} = amplModel;
            timeWindow = a.er_data.params.timeWindow{1};
            
            for cc = 1:4
                % Get variable time points based on baseline corrected mean across voxels.
                modelCondition = amplModel{cc} - mean(amplModel{cc}((timeWindow<=0),:),1);
                meanModelCondition = nanmean(modelCondition,2)';
                timePointsTrialToAverageModel(s,cc,roiIdx,mm,:) = getTimePointsTrialToAverage(...
                    twSelectionType, timeWindow, meanModelCondition);
                modelTrialBslCorrected{s,cc,roiIdx,mm} = nanmean(modelCondition(squeeze(timePointsTrialToAverageModel(s,cc,roiIdx,mm,:)),:),1);
                modelTrial{s,cc,roiIdx, mm} = nanmean(amplModel{cc}(squeeze(timePointsTrialToAverageModel(s,cc,roiIdx,mm,:)),:),1);
                
                modelCondition = amplModel{cc+(numConditions/2)} - mean(amplModel{cc+(numConditions/2)}((timeWindow<=0),:),1);
                meanModelCondition = nanmean(modelCondition,2)';
                timePointsTrialToAverageModel(s,cc+(numConditions/2),roiIdx,mm,:) = getTimePointsTrialToAverage(...
                    twSelectionType, timeWindow, meanModelCondition);
                modelTrialBslCorrected{s,cc+(numConditions/2),roiIdx,mm} = nanmean(modelCondition(squeeze(timePointsTrialToAverageModel(s,cc+(numConditions/2),roiIdx,mm,:)),:),1);
                modelTrial{s,cc+(numConditions/2),roiIdx, mm} = nanmean(amplModel{cc+(numConditions/2)}(squeeze(timePointsTrialToAverageModel(s,cc+(numConditions/2),roiIdx,mm,:)),:),1);
                
                T_meanAmpModel = [T_meanAmpModel; ...
                    table(subjnr, cc, rois(idx), modelNames(mm), ...
                    modelTrial(s,cc,roiIdx,mm), ...
                    modelTrial(s,cc+(numConditions/2),roiIdx, mm), ...
                    'VariableNames',{'Subject','Condition','ROI','modelType','SEQ_ModelAmplitude','SIM_ModelAmplitude'})];
                T_meanAmpModelBslCorrected = [T_meanAmpModelBslCorrected; ...
                    table(subjnr, cc, rois(idx),modelNames(mm), ...
                    modelTrialBslCorrected(s,cc,roiIdx, mm), ...
                    modelTrial(s,cc+(numConditions/2),roiIdx, mm), ...
                    'VariableNames',{'Subject','Condition','ROI','modelType','SEQ_ModelAmplitude_bslcorr','SIM_ModelAmplitude_bslcorr'})];
                
            end
        end
        
        amplData        = a.T.DataTSCV;
        allAmpData{s,roiIdx} = amplData;
        timeWindow = a.er_data.params.timeWindow{1};
        
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
            T_meanAmp = [T_meanAmp; table(subjnr, cc, rois(idx),...
                amplTrial(s,cc,roiIdx), ...
                amplTrial(s,cc+(numConditions/2),roiIdx), ...
                'VariableNames',{'Subject','Condition','ROI','SEQ_MeanAmplitude','SIM_MeanAmplitude'})];
            T_meanAmpBslCorrected = [T_meanAmpBslCorrected; table(subjnr, cc, rois(idx), ...
                amplTrialBslCorrected(s,cc,roiIdx), ...
                amplTrialBslCorrected(s,cc+(numConditions/2),roiIdx), ...
                'VariableNames',{'Subject','Condition','ROI','SEQ_MeanAmplitude_bslcorr','SIM_MeanAmplitude_bslcorr'})];
            
            if mm==3
                if size(a.data.splitHalfRel.both.trialCV,2)~=8
                    a.data.splitHalfRel.both.trialCV = a.data.splitHalfRel.both.trialCV';
                end
                noiseceilingTrial{s,cc,roiIdx} = a.data.splitHalfRel.both.trialCV(:,cc);
                noiseceilingTrial{s,cc+(numConditions/2),roiIdx} =  ...
                    a.data.splitHalfRel.both.trialCV(:,cc+(numConditions/2));
            end
        end
        
    end
    
end

if storeData
    % Save summary data we so don't have to compute it next time
    sNames = sprintf('S%i_',subjnrs);
    allROIs = pths.allROIs;
    save(fullfile(projectDir,'data/simseq/average',...
        sprintf('%sBlockAmps_pRFparams_modelFit_results_cv_variableBlockOnset%s.mat',sNames,fitName)),...
        'T_meanAmp','T_meanAmpModel','T_R2_Beta_pRFs','paramsROI','noiseceilingRun',...
        'allAmpData','amplTrial','allAmpModel','modelTrial','alignMeROI',...
        'r2Trials','betaRuns','r2Runs','conditionData',...
        'subjnrs','allROIs','timePointsTrialToAverage','timePointsTrialToAverageModel',...
        'spatialModels','temporalModels','timePointsTrialToAverage','-v7.3');
end



