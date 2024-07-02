% s_selectAndSaveCSTfixPredictionsMatchingVoxels.m
projectDir = fullfile(simseqRootPath);
baseDir = fullfile(simseqRootPath,'simseq/results/');
subDir = 'preprocData_stRet_matchingVoxels_1ch-glm';

subSaveFolder   = 'savedPredictions_fixedPRF_combT';
fitName         = {'gridFit2','stRet','stRet'};
betaFolder      = 'betaR2_XValFit_OLS_gridFit2';

stimFileName    = 'stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun';
subLoadFolder   = 'varyStimDur_x_Area4_2x2design_4deg2_trimmedGaussian_bothHemi_mid_roiCorrected_fullrun2';
postFixTtle     = [];
roiType0        = 'stimcorner4_area4sq_eccen5'; %     % Use stimulus corner, can also be one square (e.g.: 'square1')
hemi            = 'both';
numConditions   = 8;
hemis = {'lh','rh'};

subjnrs = [1:3,7:10];

for subjnr = subjnrs

    sesNr = getSessionNrMainExp(subjnr);
    s = (subjnr==subjnrs);
    pths = getSubjectPaths(projectDir,subjnr,sesNr);
    fname = sprintf('%s2_v%d.mat', ...
        stimFileName, pths.expversionNr);
    rois       = pths.definedROIsBOTH;
    saveFolder = fullfile(baseDir,pths.subjID, [subSaveFolder '_matchingVoxels'], ...
        [betaFolder '_1ch-glm']);
    
    for r = 1:length(rois)
        
        fnameA = fullfile(baseDir, pths.subjID,subDir, sprintf('%s_v%d',subLoadFolder,pths.expversionNr), ...
            sprintf('preprocData_%s_%s_stRet_CSTopt_DNST_matchingVoxels_cssFit_run12.mat',rois{r},roiType0));
        if exist(fullfile(fnameA), 'file')
            
            % Get preproc data (from ROI with matching voxels between
            % datasets, and checked for pRF locations)
            a = load(fnameA);
            
            % Get best grid fit prediction and data
            b = load(fullfile(baseDir,pths.subjID,subSaveFolder, betaFolder, ...
                sprintf('betaValsFull_%s_%s_onegaussianFit_3ch-stLN_run12_offsetFlag1_OLS_bslCorrected_expBestR2.mat',rois{r},roiType0)));
            
            % Get overlapping ccoordinates for left and right hemispheres, 
            % between the ROI with matching voxels and the "full" simseq ROI
            [coords_lh, ai_lh, bi_lh] = intersectCols(a.alignMe.commonCoords.lh, b.optimal.alignMe.commonCoords.lh);
            [coords_rh, ai_rh, bi_rh] = intersectCols(a.alignMe.commonCoords.rh, b.optimal.alignMe.commonCoords.rh);
            
            % Check if the coords from the dataruns, pRFs and commonCoords
            % are the same, for both left and right hemispheres.
            if isfield(a.alignMe.coords.dataRuns, 'lh') && ~isempty(a.alignMe.coords.dataRuns.lh)
                assert(isequal(a.alignMe.coords.dataRuns.lh, a.alignMe.commonCoords.lh))
            elseif isempty(a.alignMe.coords.dataRuns.lh)
                assert(size(a.alignMe.commonCoords.lh,2)==0)
            end
            if isfield(a.alignMe.coords.dataRuns, 'rh') && ~isempty(a.alignMe.coords.dataRuns.rh)
                assert(isequal(a.alignMe.coords.dataRuns.rh, a.alignMe.commonCoords.rh))
            elseif isempty(a.alignMe.coords.dataRuns.rh)
                assert(size(a.alignMe.commonCoords.rh,2)==0)
            end

            % Get matching subset of best fitting voxels from grid fit
            [coords_both, ai_both, bi_both] = intersectCols(...
                cat(2,a.alignMe.commonCoords.lh,a.alignMe.commonCoords.rh),...
                cat(2,b.optimal.alignMe.commonCoords.lh,b.optimal.alignMe.commonCoords.rh));
            
            % Now apply voxel index to gridfit data
            optimal = struct();
            optimal.alignMe.commonCoords.lh = coords_lh;
            optimal.alignMe.commonCoords.rh = coords_rh;
            optimal.alignMe.coords.dataRuns.lh = a.alignMe.commonCoords.lh(:,ai_lh);
            optimal.alignMe.coords.dataRuns.rh = a.alignMe.commonCoords.rh(:,ai_rh);
            optimal.alignMe.coords.pRFs.lh = b.optimal.alignMe.commonCoords.lh(:,bi_lh);
            optimal.alignMe.coords.pRFs.rh = b.optimal.alignMe.commonCoords.rh(:,bi_rh);
            
            optimal.alignMe.dataRuns.selected.lh = b.optimal.alignMe.dataRuns.selected.lh(ai_lh);
            optimal.alignMe.dataRuns.selected.rh = b.optimal.alignMe.dataRuns.selected.rh(ai_rh);
            optimal.alignMe.pRFs.selected.lh = b.optimal.alignMe.pRFs.selected.lh(bi_lh);
            optimal.alignMe.pRFs.selected.rh = b.optimal.alignMe.pRFs.selected.rh(bi_rh);
            
            optimal.roi_exp_idx = b.optimal.roi_exp_idx(bi_both,:);
            optimal.alignMe = a.alignMe;
            for ii = 2:length(b.optimal.T.Properties.VariableNames)
                tabName = b.optimal.T.Properties.VariableNames{ii};
                for nn = 1:length(b.optimal.T.(tabName))
                    optimal.T.(tabName){nn} = b.optimal.T.(tabName){nn}(:,bi_both);
                    optimal.T_orig.DataTS{nn} = b.optimal.T_orig.DataTS{nn}(:,bi_both);
                    optimal.T_orig.DataTSerror{nn} = b.optimal.T_orig.DataTSerror{nn}(:,bi_both);
                end
            end
            
            optimal.B_crossval_mean = b.optimal.B_crossval_mean(bi_both,:);
            optimal.B_crossval = b.optimal.B_crossval(bi_both,:,:);
            optimal.R2_crossval = b.optimal.R2_crossval(:,bi_both);
            optimal.R2_crossval_trial = b.optimal.R2_crossval_trial(:,bi_both);
            optimal.R2_full_trial = b.optimal.R2_full_trial(:,bi_both);
            optimal.Y_crossval_mean = b.optimal.Y_crossval_mean(:,bi_both);
            optimal.Y_crossval{1} = b.optimal.Y_crossval{1}(:,bi_both);
            optimal.Y_crossval{2} = b.optimal.Y_crossval{2}(:,bi_both);
            optimal.sumChannelPredictionCrossvalMean = b.optimal.sumChannelPredictionCrossvalMean(:,bi_both);
            optimal.sumChannelPredictionNoOffsetCrossvalMean = b.optimal.sumChannelPredictionNoOffsetCrossvalMean(:,bi_both);
            optimal.stimRunUnique = b.optimal.stimRunUnique;
            optimal.full.params0 = b.optimal.params0;
            optimal.full.alignMe = b.optimal.alignMe;
            optimal.full.resultsCValFinalFit = b.optimal.resultsCValFinalFit;
            
            assert(isequal(size(optimal.B_crossval_mean,1), (length(optimal.alignMe.dataRuns.selected.lh)+length(optimal.alignMe.dataRuns.selected.rh))))
            optimal.er_data_pred = b.optimal.er_data;
            optimal.er_data_preprocdata = a.er_data;
            optimal.data.rf.both = b.optimal.data.rf.both(:,:,bi_both);
            if isfield(optimal.data.rf,'lh')
                optimal.data.rf.lh = b.optimal.data.rf.lh(:,:,bi_lh);
                optimal.data.dataRuns.lh = b.optimal.data.dataRuns.lh(:,:,bi_lh);
                optimal.data.predRuns.lh = b.optimal.data.predRuns.lh(:,bi_lh,:,:);
            end
            if isfield(optimal.data.rf,'rh')
                optimal.data.rf.rh = b.optimal.data.rf.rh(:,:,bi_rh);
                optimal.data.dataRuns.rh = b.optimal.data.dataRuns.rh(:,:,bi_rh);
                optimal.data.predRuns.rh = b.optimal.data.predRuns.rh(:,bi_rh,:,:);
            end
            
            optimal.data.dataRuns.both = b.optimal.data.dataRuns.both(:,:,bi_both);
            optimal.data.dataRuns.detrendedStimRuns = b.optimal.data.dataRuns.detrendedStimRuns(:,bi_both,:);
            optimal.data.predRuns.both = b.optimal.data.predRuns.both(:,bi_both,:,:);
            optimal.data.predictionCatRuns = b.optimal.data.predictionCatRuns(:,bi_both,:);
            optimal.data.detrendedData.both = b.optimal.data.detrendedData.both(:,bi_both,:);
            optimal.data.meanDetrendedData =b.optimal.data.meanDetrendedData(:,bi_both);
            optimal.data.meanDetrendedDataOdd = b.optimal.data.meanDetrendedDataOdd(:,bi_both);
            optimal.data.meanDetrendedDataEven =b.optimal.data.meanDetrendedDataEven(:,bi_both);
            optimal.data.splitHalfRel.both.run = b.optimal.data.splitHalfRel.both.run(:,bi_both);
            optimal.data.splitHalfRel.both.trial = b.optimal.data.splitHalfRel.both.trial(bi_both,:);
            optimal.data.splitHalfRel.both.runCV = b.optimal.data.splitHalfRel.both.runCV(:,bi_both);
            optimal.data.splitHalfRel.both.trialCV = b.optimal.data.splitHalfRel.both.trialCV(bi_both,:);
            
            for jj = 1:8
                optimal.allDataTrialsCrossvalMean{jj} = b.optimal.allDataTrialsCrossvalMean{jj}(:,:,bi_both);
                optimal.allDataTrials{jj} = b.optimal.allDataTrials{jj}(:,:,bi_both);
                optimal.allModelNoOffsetTrials{jj} = b.optimal.allModelNoOffsetTrials{jj}(:,:,bi_both);
                optimal.allModelTrials{jj} = b.optimal.allModelTrials{jj}(:,:,bi_both);
                optimal.allModelTrialsCrossvalMean{jj} = b.optimal.allModelTrialsCrossvalMean{jj}(:,:,bi_both);
                optimal.allModelTrialsCrossvalMeanNoOffset{jj} = b.optimal.allModelTrialsCrossvalMeanNoOffset{jj}(:,:,bi_both);
                
                
            end
            
            optimal.full.data.params = b.optimal.data.params;
            optimal.prfParams_matching.lh = a.prfParams.lh;
            optimal.prfParams_matching.rh = a.prfParams.rh;
            
            optimal.data.params.x0.both = cat(2,a.prfParams.lh.x0,a.prfParams.rh.x0);
            optimal.data.params.y0.both = cat(2,a.prfParams.lh.y0,a.prfParams.rh.y0);
            optimal.data.params.varexpl.both = cat(2,a.prfParams.lh.varexpl,a.prfParams.rh.varexpl);
            optimal.data.params.effectiveSize.both = cat(2,a.prfParams.lh.effectiveSize,a.prfParams.rh.effectiveSize);
            optimal.data.params.exponent.both = cat(2,a.prfParams.lh.exponent,a.prfParams.rh.exponent);
            optimal.data.params.exponent_temporal.both = optimal.roi_exp_idx(:,2);
            
            
            if ~exist(fullfile(saveFolder),'dir'), mkdir(saveFolder); end
            save(fullfile(saveFolder, sprintf('betaValsFull_%s_%s_stRet_CSTopt_DNST_matchingVoxels_onegaussianFit_3ch-stLN_run12_offsetFlag1_OLS_bslCorrected_expBestR2.mat',rois{r},roiType0)), ...
                'optimal')
        end
    end
end













