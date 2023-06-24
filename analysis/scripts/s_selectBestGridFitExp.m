%% s_selectBestGridFitExp.m
% Main script to select best exponent parameters based on R2 of each
% simulated prediction (ranging from [0.1:0.1:1].
%
% Main experiment subject numbers = [1,2,3,7,9,10,11,12]
%
% Written by Eline R Kupers, Stanford U 2021


projectDir     = simseqRootPath;
subjnrs        = [1,2,3,7,8,9,10,11,12,13];
exponents      = [0.1:0.05:1];
subSaveFolder  = '_fixedPRF_combT';
subBetaFolder  = 'betaR2_XValFit_OLS_gridFit';
stimFileName   = 'stim_simseq_run';
subLoadFolder  = 'varyStimDur_x_Area4_2x2design_fullrun2';
postFixTtle    = [];
roiType        = 'stimcorner4_area4sq_eccen5'; 
hemi           = 'both';
numConditions  = 8;
spatialModel   = 'onegaussianFit';  % Choose from'cssFit' or 'onegaussianFit'
temporalModel  = '3ch-stLN';
regressionType = '_OLS';
bslpostfix     = '_bslCorrected';
offsetFlag     = true;
recomputeData  = true;
saveSINGLESubjectExpR2Data  = true;
saveALLSubjectsExpR2Data  = true;
gridFitFolder = 'gridFit_Results';

%% Recompute or load grid fit R2
if recomputeData
    for s = 1:length(subjnrs)
        fprintf('Get summary stats of gridfits subject %d \n', subjnrs(s))
        subjnr = subjnrs(s);
        sesNr = getSessionNrMainExp(subjnr);
                
        % Get subjects session directories
        pths = getSubjectPaths(projectDir,subjnr,sesNr);
        simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID, pths.session);
        cd(simSeqSessionDir)
        
        subLoadFolder2 = sprintf('%s_v%d',subLoadFolder,pths.expversionNr);
        rois       = pths.definedROIsBOTH;
        
        for r = 1:length(rois)
            idx = find(strcmp(rois{r},pths.allROIs));
            for expn = 1:length(exponents)
                
                % Folder with prediction and data
                predFolder = fullfile(pths.simseqResultsDir,...
                    ['savedPredictions' subSaveFolder],subLoadFolder2);
                betaFolder = fullfile(pths.simseqResultsDir,...
                    ['savedPredictions' subSaveFolder], subBetaFolder);
                % Load data and model prediction
                loadFile = fullfile(betaFolder,sprintf('betaValsFull_%s_%s_%s_%s_run12_addOffsetFlag%d%s%s_exp%1.2f.mat', ...
                    rois{r},roiType,spatialModel, temporalModel, offsetFlag, regressionType, bslpostfix, exponents(expn)));
                a = load(loadFile);
                
                betaRunsAll{s,idx,expn}  = a.B_crossval_mean; % or B_full??
                r2RunsAll{s,idx,expn}    = a.R2_crossval_mean; % or R2_full??
                dataRunsAll{s,idx,expn}  = a.Y_crossval_mean; % or data.meanDetrendedData??
                modelRunsAll{s,idx,expn} = a.sumChannelPredictionCrossvalMean; % or sumChannelPrediction??
                
                meanR2RunXvalAll(s,idx,expn)     = mean(a.R2_crossval_mean,2,'omitnan');
                seR2RunXvalAll(s,idx,expn)       = std(a.R2_crossval_mean,[],2,'omitnan')/sqrt(size(a.R2_crossval_mean,2));
                splithalfrelRunXValAll(s,idx,expn) = nanmedian(a.data.splitHalfRel.both.runCV);
                %             nanmean(simseq_getSplitHalfDataReliability( ...
                %                 permute(a.data.dataRuns.detrendedStimRuns(:,:,1:2:end),[1,3,2]), 'runFlag', true,...
                %                 'splitRuns',[]));
                
            end
        end
        if saveSINGLESubjectExpR2Data
            betaRuns = squeeze(betaRunsAll(s,:,:));
            dataRuns = squeeze(dataRunsAll(s,:,:));
            meanR2RunXval = squeeze(meanR2RunXvalAll(s,:,:));
            modelRuns = squeeze(modelRunsAll(s,:,:));
            r2Runs = squeeze(r2RunsAll(s,:,:));
            seR2RunXval = squeeze(seR2RunXvalAll(s,:,:));
            splithalfrelRunXVal = squeeze(splithalfrelRunXValAll(s,:,:));

            save(fullfile(projectDir,'data','simseq','group',gridFitFolder,...
                sprintf('S%d_gridfit2_Exp_results.mat',subjnrs(s))),...
                'betaRuns','r2Runs','dataRuns','modelRuns','meanR2RunXval','seR2RunXval',...
                'splithalfrelRunXVal','-v7.3');
        end
    end % subjnrs
    
    if saveALLSubjectsExpR2Data
        % Rename and save all subjects data in one file
        betaRuns = betaRunsAll;
        dataRuns = dataRunsAll;
        meanR2RunXval = meanR2RunXvalAll;
        modelRuns = modelRunsAll;
        r2Runs = r2RunsAll;
        seR2RunXval = seR2RunXvalAll;
        splithalfrelRunXVal = splithalfrelRunXValAll;
        
        sNames = sprintf('S%i_',subjnrs);
        save(fullfile(projectDir,'data','simseq','group',...
            sprintf('%sgridfitExp_results.mat',sNames)),...
            'betaRuns','r2Runs','dataRuns','modelRuns','meanR2RunXval','seR2RunXval',...
            'splithalfrelRunXVal','-v7.3');
    end
end % recomputeData flag



%% Create new file combining the best fitting model predictions per voxel

for s = 1:length(subjnrs)
    subjnr = subjnrs(s);
    fprintf('Finding best fitting model exponent for subject S%02d\n',subjnr);
    sesNr = getSessionNrMainExp(subjnr);
    pths = getSubjectPaths(projectDir, subjnr,sesNr);

    simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID, pths.session);
    cd(simSeqSessionDir)
    
    subLoadFolder2 = sprintf('%s_v%d',subLoadFolder,pths.expversionNr);
    
    % Folder with prediction and data
    predFolder = fullfile(pths.simseqResultsDir,...
        ['savedPredictions' subSaveFolder],subLoadFolder2);
    betaFolder = fullfile(pths.simseqResultsDir,...
        ['savedPredictions' subSaveFolder],subBetaFolder);
    load(fullfile(projectDir,'data','simseq','group',gridFitFolder, sprintf('S%d_gridfit2_Exp_results.mat',subjnr)));
    
    
    rois = pths.definedROIsBOTH;
    for r = 1:length(rois)
        r_idx = find(strcmp(rois{r},pths.allROIs));

        if ~isempty(cell2mat(r2Runs(r_idx,:)))
            
            x = squeeze(cell2mat(r2Runs(r_idx,:)'))';
            [rmax,rmax_i] = max(x,[],2,'omitnan');
            rmax_exp = exponents(rmax_i);
            uniqueExp = unique(rmax_exp);

            optimal.roi_exp_idx = [rmax,rmax_exp'];

            for expn = 1:length(uniqueExp)

                % Load data and model prediction
                loadFile = fullfile(betaFolder,sprintf('%s_%s_%s_%s_%s_run12_addOffsetFlag%d%s%s_exp%1.2f.mat', ...
                    pths.subjID, rois{r},roiType,spatialModel, temporalModel, offsetFlag, regressionType, bslpostfix, uniqueExp(expn)));
                a = load(loadFile);

                if expn == 1
                     optimal.T = a.T;
                     optimal.T_orig = a.T_orig;
                end

                [~,idx] = find(rmax_exp==uniqueExp(expn));
                optimal.B_crossval_mean(idx,:)   = a.B_crossval_mean(idx,:);
                optimal.B_crossval(idx,:,:)      = a.B_crossval(idx,:,:);
                optimal.R2_crossval(:,idx)       = a.R2_crossval(:,idx);
                optimal.B_crossval_mean(idx,:,:) = a.B_crossval_mean(idx,:,:);
                optimal.R2_crossval_trial(:,idx) = a.R2_crossval_trial(:,idx);
                optimal.R2_full_trial(:,idx)     = a.R2_full_trial(:,idx);
                optimal.Y_crossval_mean(:,idx)   = a.Y_crossval_mean(:,idx);
                optimal.Y_crossval{1}(:,idx)     = a.Y_crossval{1}(:,idx);
                optimal.Y_crossval{2}(:,idx)     = a.Y_crossval{2}(:,idx);
                optimal.sumChannelPredictionCrossvalMean(:,idx)         = a.sumChannelPredictionCrossvalMean(:,idx);
                optimal.sumChannelPredictionNoOffsetCrossvalMean(:,idx) = a.sumChannelPredictionNoOffsetCrossvalMean(:,idx);
                optimal.resultsCValFinalFit{1}.sumChannelPrediction(:,idx) = a.resultsCValFinalFit{1}.sumChannelPrediction(:,idx);
                optimal.resultsCValFinalFit{2}.sumChannelPrediction(:,idx) = a.resultsCValFinalFit{2}.sumChannelPrediction(:,idx);
                optimal.resultsCValFinalFit{1}.sumChannelPredictionNoOffset(:,idx) = a.resultsCValFinalFit{1}.sumChannelPredictionNoOffset(:,idx);
                optimal.resultsCValFinalFit{2}.sumChannelPredictionNoOffset(:,idx) = a.resultsCValFinalFit{2}.sumChannelPredictionNoOffset(:,idx);
                optimal.resultsCValFinalFit{1}.Y_noOffset(:,idx) = a.resultsCValFinalFit{1}.Y_noOffset(:,idx);
                optimal.resultsCValFinalFit{2}.Y_noOffset(:,idx) = a.resultsCValFinalFit{2}.Y_noOffset(:,idx);
                optimal.resultsCValFinalFit{1}.R2(idx) = a.resultsCValFinalFit{1}.R2(idx);
                optimal.resultsCValFinalFit{2}.R2(idx) = a.resultsCValFinalFit{2}.R2(idx);
                optimal.resultsCValFinalFit{1}.B(idx) = a.resultsCValFinalFit{1}.B(idx);
                optimal.resultsCValFinalFit{2}.B(idx) = a.resultsCValFinalFit{2}.B(idx);
                
                optimal.stimRunUnique = a.stimRunUnique;
                optimal.params0 = a.params;
                optimal.alignMe = a.alignMe;
                optimal.er_data = a.er_data;
                optimal.data.rf = a.data.rf;
                optimal.data.params.x0.both(idx)        = a.data.params.x0.both(idx);
                optimal.data.params.y0.both(idx)        = a.data.params.y0.both(idx);
                optimal.data.params.varexpl.both(idx)   = a.data.params.varexpl.both(idx);
                optimal.data.params.effectiveSize.both(idx) = a.data.params.effectiveSize.both(idx);
                optimal.data.params.exponent.both(idx)  = a.data.params.exponent.both(idx);
                optimal.data.params.exponent_temporal.both(idx)  = a.data.params.exponent_temporal.both(idx);

                if isfield(a.data.dataRuns,'lh')
                    if ~isempty(a.data.dataRuns.lh)
                        idx_lh = idx(idx<=size(a.data.dataRuns.lh,3));
                        optimal.data.dataRuns.lh(:,:,idx_lh) = a.data.dataRuns.lh(:,:,idx_lh);
                        optimal.data.predRuns.lh(:,idx_lh,:,:) = a.data.predRuns.lh(:,idx_lh,:,:);
                    end
                end
                if isfield(a.data.dataRuns,'rh')
                    if ~isempty(a.data.dataRuns.rh)
                        if isfield(a.data.dataRuns,'lh')
                            idx_rh = idx(idx>size(a.data.dataRuns.lh,3) & (idx<=size(a.data.dataRuns.both,3)));
                            idx_rh = idx_rh - size(a.data.dataRuns.lh,3);
                        else
                            idx_rh = (idx<=size(a.data.dataRuns.both,3));
                        end
                        optimal.data.dataRuns.rh(:,:,idx_rh) = a.data.dataRuns.rh(:,:,idx_rh);            
                        optimal.data.predRuns.rh(:,idx_rh,:,:) = a.data.predRuns.rh(:,idx_rh,:,:);
                    end
                end
                optimal.data.dataRuns.both(:,:,idx) = a.data.dataRuns.both(:,:,idx);
                optimal.data.dataRuns.detrendedStimRuns(:,idx,:) = a.data.dataRuns.detrendedStimRuns(:,idx,:);
                optimal.data.detrendedData.both(:,idx,:) = a.data.detrendedData.both(:,idx,:);
                optimal.data.meanDetrendedData(:,idx) = a.data.meanDetrendedData(:,idx);
                optimal.data.meanDetrendedDataOdd(:,idx) = a.data.meanDetrendedDataOdd(:,idx);
                optimal.data.meanDetrendedDataEven(:,idx) = a.data.meanDetrendedDataEven(:,idx);
                optimal.data.splitHalfRel.both.run(idx) = a.data.splitHalfRel.both.run(idx);
                optimal.data.splitHalfRel.both.trial(idx,:) = a.data.splitHalfRel.both.trial(idx,:);
                optimal.data.splitHalfRel.both.runCV(idx) = a.data.splitHalfRel.both.runCV(idx);
                optimal.data.splitHalfRel.both.trialCV(idx,:) = a.data.splitHalfRel.both.trialCV(idx,:);
                optimal.data.predRuns.both(:,idx,:,:) = a.data.predRuns.both(:,idx,:,:);
                optimal.data.predictionCatRuns(:,idx,:) = a.data.predictionCatRuns(:,idx,:);
                
                for cc = 1:numConditions
                    optimal.allDataTrialsCrossvalMean{cc}(:,:,idx) = a.allDataTrialsCrossvalMean{cc}(:,:,idx);
                    optimal.allDataTrials{cc}(:,:,idx) = a.allDataTrials{cc}(:,:,idx);
                    optimal.allDataTrialsCrossvalMean{cc}(:,:,idx) = a.allDataTrialsCrossvalMean{cc}(:,:,idx);
                    optimal.allModelNoOffsetTrials{cc}(:,:,idx) = a.allModelNoOffsetTrials{cc}(:,:,idx);
                    optimal.allModelTrials{cc}(:,:,idx) = a.allModelTrials{cc}(:,:,idx);
                    optimal.allModelTrialsCrossvalMean{cc}(:,:,idx) = a.allModelTrialsCrossvalMean{cc}(:,:,idx);
                    optimal.allModelTrialsCrossvalMeanNoOffset{cc}(:,:,idx) = a.allModelTrialsCrossvalMeanNoOffset{cc}(:,:,idx);

                    for ll = 2:length(a.T.Properties.VariableNames)
                        optimal.T{cc,ll}{1}(:,idx) = a.T{cc,ll}{1}(:,idx);

                        if ll<=size(a.T_orig,2)
                            optimal.T_orig{cc,ll}{1}(:,idx) = a.T_orig{cc,ll}{1}(:,idx);
                        end
                    end
                end
            end

            % Save optimal exponent voxel data and model prediction
            saveFile = fullfile(betaFolder,sprintf('%s_%s_%s_%s_%s_run12_offsetFlag%d%s%s_expBestR2.mat', ...
                pths.subjID, rois{r},roiType,spatialModel, temporalModel, offsetFlag, regressionType, bslpostfix));
            save(saveFile, 'optimal','-v7.3');

            clear optimal
        end
    end
end
