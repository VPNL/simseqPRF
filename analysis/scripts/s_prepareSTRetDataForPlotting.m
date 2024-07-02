%% s_prepareSTRetDataForPlotting.m
% Script to prepare data for supplementary figures 8 and 9 in the paper:
% Title:   Rethinking simultaneous suppression in visual cortex via 
%          compressive spatiotemporal population receptive fields.
% Authors: Kupers, Kim, Grill-Spector (2024).
% Journal: Nature Communications
% DOI:     XXX
%
% Requires getting MRI data from OSF (see downloadDataTableFromOSF.m)
%
% Code written by E.R. Kupers (2024) Stanford University

%% Define params
subjnrs         = [1,2,3,7,8,9,10];
projectDir      = fullfile(simseqRootPath);
saveFigDir      = fullfile(simseqRootPath,'figures');

% Define colors, labels, plot order
spatialModels        = {'onegaussianFit','onegaussianFit','onegaussianFit'};  % Choose from'cssFit' or 'onegaussianFit'
temporalModels       = {'3ch-stLN','3ch-stLN','1ch-dcts'}; % Choose from '3ch-stLN' or '1ch-dcts'
modelNames           = {'CSTfix','CSTopt','DNST'};
cmapROIs             = getROISummaryColors(0);
conditionNamesSimSeq = {'Small-0.2s','Small-1s','Large-0.2s','Large-1s'};
nrConditions         = length(conditionNamesSimSeq);
conditionOrderSimSeq = [1:4];
conditionNamesSimSeq = conditionNamesSimSeq(conditionOrderSimSeq);
labelsComb           = {'V1','V2','V3','hV4','VO1/2','LO1/2','TO1/2','V3AB','IPS0/1'};
roiOrder             = [1:6,8:11,7,12,13]; %V1-V02, LO1/2 + TO1/2, V3AB, IPS0/1
combROIs             = [5,8,10,12]; %VO1/2 (5,6), LO1/2(8,9), TO1/2 (10,11), IPS0/1 (12,13)

% Define noise and variance explained thresholds:
nc_thresh            = 0.1; % threshold in percentage noise ceiling from simseq exp
ve_thresh            = 0.2; % threshold in percentage of pRF var expl. from toonotopy

% Define folder & file names
subSaveFolder   = 'simseq';


stimFileName    = '_cvfits';
roiType0        = 'stimcorner4_area4sq_eccen5_stRet_CSTopt_DNST_matchingVoxels'; %     % Use stimulus corner, can also be one square (e.g.: 'square1')

% Set boolean options
storeData      = true;  % store data table at the end
sNames         = sprintf('S%i_',subjnrs);
saveFileName   = fullfile(projectDir,'simseq/results/group/',...
                    sprintf('%sstRet_params.mat',sNames));



%%
T_meanAmp = table();
T_meanAmpBslCorrected = table();
T_meanAmpModel = table();
T_meanAmpModelBslCorrected = table();
T_R2_Beta_pRFs = table();

for subjnr = subjnrs(3)
    
    % Get data session nr
    s = (subjnr==subjnrs);
    pths = getSubjectPaths(projectDir,subjnr);
    simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID);
    
    % Define subfolder and ROIs
    rois       = pths.definedROIsBOTH;
    
    % Define session dir
    cd(simSeqSessionDir);
    
    for idx = 1:length(rois)
        
        roiIdx = find(strcmp(pths.allROIs,rois{idx}));
        
        % Get filename with modelfit results
        loadFileCST_fix = fullfile(simSeqSessionDir,...
            sprintf('%s_CSTfix_cvfits',pths.subjID), ...
            sprintf('%s_%s_CSTfix_cvfitResults.mat', ...
            pths.subjID,rois{idx}));
        
        loadFileCST_opt = fullfile(simSeqSessionDir,...
            sprintf('%s_CSTopt_cvfits',pths.subjID), ...
            sprintf('%s_%s_CSTopt_cvfitResults.mat', ...
            pths.subjID,rois{idx}));  
        
        loadFileDNST_opt = fullfile(simSeqSessionDir,...
            sprintf('%s_DNST_cvfits',pths.subjID), ...
            sprintf('%s_%s_DNST_cvfitResults.mat', ...
            pths.subjID,rois{idx})); 
        
        if exist(loadFileCST_fix,'file')
            % Load data and model prediction
            a = load(loadFileCST_fix);
            b = load(loadFileCST_opt);
            c = load(loadFileDNST);
        end
        
        [ab_coords_lh, ab_i_lh, ba_i_lh] = intersectCols(a.alignMe.commonCoords.lh, b.alignMe.commonCoords.lh);
        [ab_coords_rh, ab_i_rh, ba_i_rh] = intersectCols(a.alignMe.commonCoords.rh, b.alignMe.commonCoords.rh);
        
        [ac_coords_lh, ac_i_lh, ca_i_lh] = intersectCols(a.alignMe.commonCoords.lh, c.alignMe.commonCoords.lh);
        [ac_coords_rh, ac_i_rh, ca_i_rh] = intersectCols(a.alignMe.commonCoords.rh, c.alignMe.commonCoords.rh);
       
        [bc_coords_lh, bc_i_lh, cb_i_lh] = intersectCols(b.alignMe.commonCoords.lh, c.alignMe.commonCoords.lh);
        [bc_coords_rh, bc_i_rh, cb_i_rh] = intersectCols(b.alignMe.commonCoords.rh, c.alignMe.commonCoords.rh);
        
        assert(isequal(ab_i_lh,ac_i_lh))
        assert(isequal(ab_i_rh,ac_i_rh))
        assert(isequal(ba_i_lh,ca_i_lh))
        assert(isequal(ba_i_rh,ca_i_rh))
        
        model_idx_both{1} = [ab_i_lh;length(ab_i_lh)+ab_i_rh]; 
        model_idx_both{2} = [ba_i_lh;length(ba_i_lh)+ba_i_rh];
        model_idx_both{3} = [ca_i_lh;length(ca_i_lh)+ca_i_rh];
        
        model_idx_lhrh{1} = {ab_i_lh;ab_i_rh}; 
        model_idx_lhrh{2} = {ba_i_lh;ba_i_rh};
        model_idx_lhrh{3} = {ca_i_lh;ca_i_rh};
        
        % Check if the coords from the dataruns, pRFs and commonCoords
        % are the same, for both left and right hemispheres.
        if isfield(a.alignMe.coords.dataRuns, 'lh') && ~isempty(a.alignMe.coords.dataRuns.lh)
            assert(isequal(a.alignMe.coords.dataRuns.lh, a.alignMe.commonCoords.lh))
            %                 assert(isequal(a.alignMe.coords.dataRuns.lh, a.alignMe.coords.pRFs.lh))
        elseif isempty(a.alignMe.coords.dataRuns.lh)
            assert(size(a.alignMe.commonCoords.lh,2)==0)
        end
        if isfield(a.alignMe.coords.dataRuns, 'rh') && ~isempty(a.alignMe.coords.dataRuns.rh)
            assert(isequal(a.alignMe.coords.dataRuns.rh, a.alignMe.commonCoords.rh))
            %                 assert(isequal(a.alignMe.coords.dataRuns.rh, a.alignMe.coords.pRFs.rh))
        elseif isempty(a.alignMe.coords.dataRuns.rh)
            assert(size(a.alignMe.commonCoords.rh,2)==0)
        end
        %
        
        
        for mm = 1:length(spatialModels)
            
            switch mm
                case 1
                    dat = a;
                case 2
                    dat = b;
                case 3 
                    dat = c;
            end
                
            
            % Check if mean cv R^2 is computed, and compute if not
            if ~isfield(dat,'R2_crossval_mean')
                dat.R2_crossval_mean = mean(dat.R2_crossval,1,'omitnan');
            end
            % Clear some memory
            dat.data.params.stim.images = [];
            dat.data.params.stim.images_unconvolved = [];
            
            % Store data, model predictions, betas and R2
            
            dataRuns{s,roiIdx,mm}  = dat.Y_crossval_mean(:,model_idx_both{mm});
            modelRuns{s,roiIdx,mm} = dat.sumChannelPredictionCrossvalMean(:,model_idx_both{mm});
            betaRuns{s,roiIdx,mm}       = dat.B_crossval_mean(model_idx_both{mm},[1:end-1])';
            r2Runs{s,roiIdx,mm}         = dat.R2_crossval_mean(:,model_idx_both{mm});
            r2Trials{s,roiIdx,mm}       = dat.R2_crossval_trial(:,model_idx_both{mm});
            
            % Add noise ceiling
                noiseceilingRun{s,roiIdx,mm} = dat.data.splitHalfRel.both.runCV(model_idx_both{mm});
            
            % Add params
            params.x0.both            = dat.data.params.x0.both(model_idx_both{mm});
            params.y0.both            = dat.data.params.y0.both(model_idx_both{mm});
            params.effectiveSize.both = dat.data.params.effectiveSize.both(model_idx_both{mm});
            params.varexpl.both       = dat.data.params.varexpl.both(model_idx_both{mm});
            params.exp_spatial.both   = dat.data.params.exponent.both(model_idx_both{mm});
            if strcmp(modelNames{mm},'DN_ST')
                params.temporal.tau1      = dat.data.params.temporal.param.tau1(model_idx_both{mm});
                params.temporal.weight    = dat.data.params.temporal.param.weight(model_idx_both{mm});
                params.temporal.tau2      = dat.data.params.temporal.param.tau2(model_idx_both{mm});
                params.temporal.n         = dat.data.params.temporal.param.n(model_idx_both{mm});
                params.temporal.sigma     = dat.data.params.temporal.param.sigma(model_idx_both{mm});
            elseif strcmp(modelNames{mm},'CST_fix')||strcmp(modelNames{mm},'CST_opt')
                params.exp_temporal.both  = dat.data.params.exponent_temporal.both(model_idx_both{mm});
                if strcmp(modelNames{mm},'CST_opt')
                    params.temporal.exponent   = dat.data.params.temporal.param.exponent(model_idx_both{mm});
                    params.temporal.tau_s      = dat.data.params.temporal.param.tau_s(model_idx_both{mm});
                    params.temporal.tau_t      = dat.data.params.temporal.param.tau_t(model_idx_both{mm});
                end
            end
            
            paramsROI{s,roiIdx,mm}       = params;
            
            % Convert all to table
            T_R2_Beta_pRFs = [T_R2_Beta_pRFs; table(subjnr, rois(idx), modelNames(mm), betaRuns(s,roiIdx,mm), ...
                r2Runs(s,roiIdx,mm), r2Trials(s,roiIdx,mm), {NaN}, ...
                'VariableNames',{'Subject','ROI','modelType','betaRun','r2Run','r2Trial','pRFsParams'})];
            
            % Get DATA SIM-SEQ mean stim block amplitudes
            for ii = 1:8
                amplModel{ii}       = dat.T.ModelTSCVMn{ii}(:,model_idx_both{mm});
            end
            allAmpModel{s,roiIdx,mm} = amplModel; 
            if mm == 1
                timeWindow = dat.er_data_pred.params.timeWindow{1};
            else
                timeWindow = dat.er_data.params.timeWindow{1};
            end
            
            % Get MODEL-based SIM-SEQ mean stim block amplitudes
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
            
            
            if ~isempty(dat)
                for ii = 1:8
                    amplData{ii} = dat.T.DataTSCV{ii}(:,model_idx_both{mm}); 
                end
                allAmpData{s,roiIdx,mm} = amplData;
                
                alignMeROI{s,roiIdx,mm} = dat.alignMe;
                alignMeROI{s,roiIdx,mm}.commonCoords.lh = alignMeROI{s,roiIdx,mm}.commonCoords.lh(:,model_idx_lhrh{mm}{1});
                alignMeROI{s,roiIdx,mm}.commonCoords.rh = alignMeROI{s,roiIdx,mm}.commonCoords.rh(:,model_idx_lhrh{mm}{2});
                
                alignMeROI{s,roiIdx,mm}.coords.dataRuns.lh = alignMeROI{s,roiIdx,mm}.coords.dataRuns.lh(:,model_idx_lhrh{mm}{1});
                alignMeROI{s,roiIdx,mm}.coords.dataRuns.rh = alignMeROI{s,roiIdx,mm}.coords.dataRuns.rh(:,model_idx_lhrh{mm}{2});
                
                alignMeROI{s,roiIdx,mm}.coords.pRFs.lh = alignMeROI{s,roiIdx,mm}.coords.pRFs.lh(:,model_idx_lhrh{mm}{1});
                alignMeROI{s,roiIdx,mm}.coords.pRFs.rh = alignMeROI{s,roiIdx,mm}.coords.pRFs.rh(:,model_idx_lhrh{mm}{2});
                
                alignMeROI{s,roiIdx,mm}.dataRuns.selected.lh = alignMeROI{s,roiIdx,mm}.dataRuns.selected.lh(model_idx_lhrh{mm}{1});
                alignMeROI{s,roiIdx,mm}.dataRuns.selected.rh = alignMeROI{s,roiIdx,mm}.dataRuns.selected.rh(model_idx_lhrh{mm}{2});
                
                alignMeROI{s,roiIdx,mm}.pRFs.selected.lh = alignMeROI{s,roiIdx,mm}.pRFs.selected.lh(model_idx_lhrh{mm}{1});
                alignMeROI{s,roiIdx,mm}.pRFs.selected.rh = alignMeROI{s,roiIdx,mm}.pRFs.selected.rh(model_idx_lhrh{mm}{2});
                
                
                if mm == 1
                    timeWindow = dat.er_data_pred.params.timeWindow{1};
                else
                    timeWindow = dat.er_data.params.timeWindow{1};
                end
                for cc = 1:4
                    % Get variable time points based on baseline corrected mean across voxels.
                    dataCondition = amplData{cc} - mean(amplData{cc}((timeWindow<=0),:),1);
                    meanDataCondition = nanmean(dataCondition,2)';
                    timePointsTrialToAverage(s,cc,roiIdx,mm,:) = getTimePointsTrialToAverage(...
                        twSelectionType, timeWindow, meanDataCondition);
                    
                    amplTrialBslCorrected{s,cc,roiIdx,mm} = nanmean(dataCondition(squeeze(timePointsTrialToAverage(s,cc,roiIdx,mm,:)),:),1);
                    amplTrial{s,cc,roiIdx,mm} = nanmean(amplData{cc}(squeeze(timePointsTrialToAverage(s,cc,roiIdx,mm,:)),:),1);
                    
                    dataCondition = amplData{cc+(numConditions/2)} - mean(amplData{cc+(numConditions/2)}((timeWindow<=0),:),1);
                    meanDataCondition = nanmean(dataCondition,2)';
                    timePointsTrialToAverage(s,cc+(numConditions/2),roiIdx,mm,:) = getTimePointsTrialToAverage(...
                        twSelectionType, timeWindow, meanDataCondition);
                    amplTrialBslCorrected{s,cc+(numConditions/2),roiIdx,mm} = nanmean(dataCondition(squeeze(timePointsTrialToAverage(s,cc+(numConditions/2),roiIdx,mm,:)),:),1);
                    amplTrial{s,cc+(numConditions/2),roiIdx,mm} = nanmean(amplData{cc+(numConditions/2)}(squeeze(timePointsTrialToAverage(s,cc+(numConditions/2),roiIdx,mm,:)),:),1);
                    
                    % Create mean amplitude table
                    T_meanAmp = [T_meanAmp; table(subjnr, cc, rois(idx),modelNames(mm),...
                        amplTrial(s,cc,roiIdx,mm), ...
                        amplTrial(s,cc+(numConditions/2),roiIdx,mm),...
                        'VariableNames',{'Subject','Condition','ROI','modelType','SEQ_MeanAmplitude','SIM_MeanAmplitude'})];
                    T_meanAmpBslCorrected = [T_meanAmpBslCorrected; table(subjnr, cc, rois(idx), modelNames(mm), ...
                        amplTrialBslCorrected(s,cc,roiIdx,mm), ...
                        amplTrialBslCorrected(s,cc+(numConditions/2),roiIdx,mm),...
                        'VariableNames',{'Subject','Condition','ROI','modelType','SEQ_MeanAmplitude_bslcorr','SIM_MeanAmplitude_bslcorr'})];
                    
                    
                    if size(dat.data.splitHalfRel.both.trialCV,2)~=8
                        dat.data.splitHalfRel.both.trialCV = dat.data.splitHalfRel.both.trialCV';
                    end
                    noiseceilingTrial{s,cc,roiIdx,mm} = dat.data.splitHalfRel.both.trialCV(model_idx_both{mm},cc);
                    noiseceilingTrial{s,cc+(numConditions/2),roiIdx,mm} =  ...
                        dat.data.splitHalfRel.both.trialCV(model_idx_both{mm},cc+(numConditions/2));
                    
                end
            end
        end
    end
    
end

if storeData
    % Save summary data we so don't have to compute it next time
    allROIs = pths.allROIs;
    save(fullfile(saveFileName),...
        'T_meanAmp','T_meanAmpModel','T_R2_Beta_pRFs','paramsROI','noiseceilingRun',...
        'allAmpData','amplTrial','allAmpModel','modelTrial','alignMeROI',...
        'r2Trials','betaRuns','r2Runs','conditionData',...
        'subjnrs','allROIs','timePointsTrialToAverage','timePointsTrialToAverageModel',...
        'spatialModels','temporalModels','timePointsTrialToAverage','-v7.3');
end



