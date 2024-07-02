%% s_prepareDataForPlotting.m
% Script to prepare data for "createDataTable" function. Mainly combines
% ROIs (like VO1 and VO2 into VO1/2) and single voxel pRF parameters with
% single voxel model fitting results.
%
% Code written by E.R. Kupers (2024) Stanford University
%
%% Define params
subjnrs        = [1,2,3,7,8,9,10,11,12,13];
projectDir     = fullfile(simseqRootPath);
spatialModels  = {'onegaussianFit','cssFit','onegaussianFit'};  % Choose from'cssFit' or 'onegaussianFit'
temporalModels = {'1ch-glm','1ch-glm','3ch-stLN'};
modelNames      = {'LSS','CSS','CST'};

recomputeData  = true; % load stored (false) or recompute values (true)
storeData      = true;  % store data when recomputing values;
offsetFlag     = true;  % Use betas from modelfit with added row of ones to allow for offset in GLM
bslCorrectFlag = true;  % Use data with baseline correction
twSelectionType = 'variableOnset';
numConditions   = 8;
timeWindow      = [-4:18];

%%
% if recomputeData
T_meanAmp = table();
T_meanAmpBslCorrected = table();
T_meanAmpModel = table();
T_meanAmpModelBslCorrected = table();
T_R2_Beta_pRFs = table();

for subjnr = subjnrs(3)
    
    % Get data session nr
    s = (subjnr==subjnrs);
    pths = getSubjectPaths(projectDir,subjnr);
    
    % Define subfolder and ROIs
    rois       = pths.definedROIsBOTH;
    
    % Define session dir
    simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID);
    cd(simSeqSessionDir);
    
    % Load trial/condition data
    conditionData.loc  = getStimLocDeg(projectDir,subjnr, true, false, 'both');
    
    for idx = 1:length(rois)
        
        roiIdx = find(strcmp(pths.allROIs,rois{idx}));
        
        for mm = 1:length(spatialModels)
            % Folder with prediction and data
            spatialModel = spatialModels{mm};
            temporalModel = temporalModels{mm};
            
            % Load data and model prediction
            loadFile = fullfile(simSeqSessionDir,...
                sprintf('%s_%s_%s_%s_cvfitResults.mat', ...
                pths.subjID, rois{idx},spatialModel, temporalModel));
            
            if exist(loadFile)
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
                    r2Runs(s,roiIdx,mm), {NaN}, ...
                    'VariableNames',{'Subject','ROI','modelType','betaRun','r2Run','pRFsParams'})];
                
                amplModel       = a.T.ModelTSCVMn;
                allAmpModel{s,roiIdx,mm} = amplModel;
                
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
    
end

if storeData
    % Save summary data we so don't have to compute it next time
    sNames = sprintf('S%i_',subjnrs);
    allROIs = pths.allROIs;
    save(fullfile(projectDir,'data/simseq/group',...
        sprintf('%sBlockAmps_pRFparams_modelFit_results_cv_variableBlockOnset.mat',sNames)),...
        'T_meanAmp','T_meanAmpModel','T_R2_Beta_pRFs','paramsROI','noiseceilingRun',...
        'allAmpData','amplTrial','allAmpModel','modelTrial',...
        'betaRuns','r2Runs','conditionData',...
        'subjnrs','allROIs','timePointsTrialToAverage','timePointsTrialToAverageModel',...
        'spatialModels','temporalModels','timePointsTrialToAverage','-v7.3');
end



