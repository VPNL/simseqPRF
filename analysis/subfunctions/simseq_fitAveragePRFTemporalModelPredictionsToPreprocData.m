function [] = simseq_fitAveragePRFTemporalModelPredictionsToPreprocData(subjnr,projectDir,spatialModel,temporalModel,varargin)
%   simseq_fitAveragePRFTemporalModelPredictionsToPreprocData(subjnr, projectDir, spatialModel, temporalModel)
%
% Function to fit generated spatiotemporal model predictions using pRF of
% single vertex from toonotopy experiment to measured fMRI data.
% NOTE: Data must be preprocessed with simseq_preprocessData.m
%
% INPUTS
% subjnr                  : (int) subject to analyze
% projectDir              : (str) where does project live
% spatialModel            : (str) spatial prf model. Choose from:
%                              'cssFit' (compressive spatial summation)
%                              'onegaussianFit' (standard 2D Gaussian)
% temporalModel           : (str) temporal prf model. Choose from:
%                              '1ch-glm' (linear),
%                              '1ch-dcts' (Divisive normalization),
%                              '2ch-exp-sig' (2 channel model)
% [dt]                    : (str) vistasoft dataTYPE
% [hemi]                  : (str) hemisphere to analyze, 'lh',
%                              'rh', or 'both' (default)
% [recomputeFlag]         : (bool) recompute (true, default) or load
%                               cached data (false) by tc_init voxel time series
% [preserveCoords]        : (bool)supposed to preserve ROI coordinates
%                               across sessions (tc init function). I say
%                              "supposed" because this functionality is buggy,
%                               hence default is set to 0, and we do it
%                               ourselves.
% [removedTRsAtStart]     : (int) how many TRs to remove at the beginning
%                               of a run (so we only model stabilized signal)
% [roiType]               : (str) within ROI, one can pick a simseq subROI.
%                               Choose from 'stimcorner' (pilot 1)
%                               stimcorner16','stimcorner4', (pilot 2)
%                               'stimcorner4_area4sq_eccen5' (default),
%                               'stimcorner4_area2sq_eccen5'
% [addOffsetFlag]            : (bool) add row of ones to allow for offset in GLM
% [sumSTFlag]             : (bool) weigh Sustained and Transient channel
%                               equally and fit their sum
% [weightChannels]        : (int) 1xN vector, with N equals number of
%                               predicted BOLD channels. Using fixed weights
%                               between [0-1] to allows us to take a
%                               weighted average of all channels in an
%                               attempt to avoid poor model fits (due to
%                               collinnearity between the channels)
% [verbose]               : (bool) plot figures or not (default=true)
% [saveFigs]              : (bool) save figures or not (default=true)
% [saveBetas]             : (bool) save chopped betas and modelfits or not (default=true)
% [loadDataFolder]        : (str) location of saved preprocessed data, if
%                               empty, we assume
%                               fullfile(pths.simseqResultsDir,'preprocData')
% [saveBetaFolder]        : (str) location where to save beta and model
%                               predictions, if empty, we assume
%                               fullfile(pths.simseqResultsDir,'savedPredictions');
% [saveFigDir]            : (str) location where to save figures, we bool
%                               is true
% [sessionNr]             : (int) number that indicates session in data folder
% [doCrossValidation]     : (bool) do 50-50 crossvalidated model fit or not
%                               (default = false)
% [useFixedExponent]      : (float) single numeric value with the
%                               fixed exponent you want to fit separately.
%                               we use this function to do a grid search.
%% 1. Set params and folders
p = inputParser;
p.addRequired('subjnr', @isnumeric);
p.addRequired('projectDir',@ischar);
p.addRequired('spatialModel', @(x) any(validatestring(x,{'cssFit','onegaussianFit','differenceOfGaussiansFit'})));
p.addRequired('temporalModel', @(x) any(validatestring(x,{'1ch-glm','1ch-dcts','2ch-exp-sig','2ch-css-sig','3ch-stLN'})));
p.addParameter('hemi','both', @(x) any(validatestring(x,{'lh','rh','both'})));
p.addParameter('removedTRsAtStart', 2, @isnumeric);
p.addParameter('roiType','stimcorner4_area4sq_eccen5',@ischar);
p.addParameter('roiIdx','all', @(x) (ischar(x) || isnumeric(x) || isscalar(x) || isvector(x)));
p.addParameter('fs',1000,@isnumeric);
p.addParameter('addOffsetFlag', true, @islogical);
p.addParameter('sumSTFlag', false, @islogical);
p.addParameter('weightChannels',[], @isnumeric);
p.addParameter('verbose', true, @islogical);
p.addParameter('saveFigs', true, @islogical);
p.addParameter('saveBetas',true, @islogical);
p.addParameter('loadDataFolder',[], @ischar);
p.addParameter('saveBetaFolder',[], @ischar);
p.addParameter('saveFigDir',[], @ischar);
p.addParameter('sessionNr',3, @isnumeric);
p.addParameter('subFolder',[], @ischar);
p.addParameter('doCrossValidation',false,@islogical);
p.addParameter('useFixedExponent', [], @isnumeric);
p.addParameter('totalBlockTRs',{-4:1:18},@iscell);
p.addParameter('baselineBlockTRs',{-4:1:0},@iscell);
p.parse(subjnr,projectDir,spatialModel,temporalModel,varargin{:});

% Rename variables
fnToRename = fieldnames(p.Results);
for ff = 1:length(fnToRename)
    eval([sprintf('%s = p.Results.(fnToRename{ff});', fnToRename{ff})]);
end
clear fnToRename ff

%% Set subject params
cd(fullfile(projectDir,'experiments/simseq/'))
pths             = getSubjectPaths(projectDir,subjnr,sessionNr);
stimRun          = unique(pths.runOrder);
simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID, pths.session);

% Define folders to save beta values
if sumSTFlag && isempty(saveBetaFolder)
    saveBetaFolder = fullfile(pths.simseqResultsDir,'savedPredictions','betaValsFull_sumST');
elseif ~sumSTFlag && isempty(saveBetaFolder)
    saveBetaFolder = fullfile(pths.simseqResultsDir,'savedPredictions');
end

% Define folders to save figures
if saveFigs && isempty(saveFigDir)
    if sumSTFlag
        saveFigDir = fullfile(pths.figureDir, 'model_vs_data_sumST');
    else
        saveFigDir = fullfile(pths.figureDir, 'model_vs_data');
    end
end

% Define subfolders
if isempty(subFolder)
    subFolderName = sprintf('%s%d_v%d',subFolder, max(stimRun), pths.expversionNr);
else
    subFolderName = subFolder;
end

% Define folder to load preprocessed data
if isempty(loadDataFolder)
    loadDataFolder = fullfile(pths.simseqResultsDir,'preprocData');
end

% Get hemis
if strcmp(hemi, 'both'); 
    hm = {'lh','rh'}; 
else hm = {hemi}; 
end

hm = {hemi}; squareNr = 1:4; 

% Get ROIs depending on selected hemi(s): left, right, or both
if strcmp(hemi,'both')
    rois       = pths.definedROIsBOTH;
elseif strcmp(hemi,'lh')
    rois       = pths.definedROIsLH;
elseif strcmp(hemi,'rh')
    rois       = pths.definedROIsRH;
end

% Check for selected ROIs
if ~strcmp(roiIdx,'all') && isnumeric(roiIdx)
    rois = rois(roiIdx);
elseif ~strcmp(roiIdx,'all') && ischar(roiIdx)
    rois = find(ismember(rois,{roiIdx}));
end

%% Loop over selected rois
for roi = 1:length(rois) % find(ismember(rois,'VO1')) % %
    
    % Load preprocessed data (see s_runDataPreproc.m)
    postFix = sprintf('%d',1:max(stimRun));

    loadFile = fullfile(loadDataFolder,...
        sprintf('%s%d_v%d',subFolderName,max(stimRun),pths.expversionNr),...
        sprintf('preprocData_%s_%s_%s_run%s_sepSquare.mat', ...
        rois{roi},roiType,'cssFit', postFix));
    load(loadFile, 'Tsquare','allDataTrialsSquare','prfParamsSquare','data',...
        'nrTrials','er_data','stimMS','alignMe','allDataTrialsSquareOdd','allDataTrialsSquareEven');
    
    % Define totalScans
    if isfield(data.dataRuns,'rh')
        if ~isempty(data.dataRuns.rh)
            totalScans = size(data.dataRuns.rh,2);
        end
    else, totalScans = size(data.dataRuns.lh,2); end
    
    %% Get predictions
    dataOrig = data; 
    TmodelFull = table(); 
    
    for sqr = squareNr
        % Preallocate space
        predRunBothHemi = []; stimRunUnique = [];
        predAllRuns = cell(1,length(stimRun));
        
        for run = 1:length(stimRun)
            if ~isempty(subFolder)
                subFolderName2 = sprintf('%s%d_v%d',subFolder, run, pths.expversionNr);
            end
            
            % Load prediction file
            if length(hm)==1
                predFile = fullfile(saveBetaFolder,subFolderName2, ...
                    sprintf('modelPredictions_%s_%s_%s_cssFit_run12_square%d_%s_%s_run%d.mat',...
                    hm{1},rois{roi},roiType,sqr,spatialModel, temporalModel, run));
            
                predAllRuns{run} = load(predFile);
                
                % Concatenate predictions
                predRunBothHemi = cat(4,predRunBothHemi, predAllRuns{run}.pred.predHRF);
                
                % Get stimulus file, so we can plot it later
                maxStimMS = floor(size(predRunBothHemi,1)-1)*fs;
                
                % .. trim to match trs
                stim = sum(predAllRuns{run}.pred.params.stim.images,1);
                
                % .. normalize to make plotting easier
                stimRunUnique(run,:) = stim./sum(stim(:));
            end
        end
        
        % Get params
        params = predAllRuns{1}.pred.params;
        params.analysis.doCrossValidation = doCrossValidation;
        params.analysis.totalScans = totalScans;
        
        % Get pRF gray coordinates
        if isfield(predAllRuns{1}.pred.params.analysis.spatial, 'lh')
            coords.pRFs.lh = predAllRuns{1}.pred.params.analysis.spatial.lh.hvolROI.coords;
            selectedpRFs.predVEmask.lh = predAllRuns{1}.pred.params.analysis.spatial.lh.veMask;
        end
        if isfield(predAllRuns{1}.pred.params.analysis.spatial, 'rh')
            coords.pRFs.rh = predAllRuns{1}.pred.params.analysis.spatial.rh.hvolROI.coords;
            selectedpRFs.predVEmask.rh = predAllRuns{1}.pred.params.analysis.spatial.rh.veMask;
        end
        
        %% Remove first couple of TRs due to scanner instability
        predRunBothHemi(1:removedTRsAtStart,:,:,:) = [];
        stimRunUnique(:,(1:(removedTRsAtStart*fs))) = [];
        
        %% Check if data have one more TR than predictions or not, if so, add NaN at the end
        if isfield(data.dataRuns,'lh') || isfield(data.dataRuns,'rh')
            if ~isempty(data.dataRuns.lh) || ~isempty(data.dataRuns.rh)
                if size(predRunBothHemi,1)<size(data.dataRuns.lh,1) || ...
                    size(predRunBothHemi,1)<size(data.dataRuns.rh,1)
                    extraDataTRs = size(data.dataRuns.lh,1)-size(predRunBothHemi,1);
                    predRunBothHemi = cat(1,predRunBothHemi,NaN(extraDataTRs,...
                        size(predRunBothHemi,2),size(predRunBothHemi,3),size(predRunBothHemi,4)));
                    
                    stimRunUnique = cat(2,stimRunUnique, NaN(size(stimRunUnique,1),(extraDataTRs*fs)));
                end
            end
        end
        
        % Pad stim run with NaNs for extra TR and if we need to
        % round up to 1s
        diffTRs =  size(predRunBothHemi,1) - size(stimRunUnique,2)/fs;
        stimRunUnique = cat(2,stimRunUnique, NaN(size(stimRunUnique,1),round(diffTRs*fs)));
        
        %% Check voxel alignment of datasets
        T_orig = Tsquare;
        data.predRuns.(hemi)           = predRunBothHemi;
        data.predRuns.both             = data.predRuns.(hemi);
        data.rf.(hemi)                 = predAllRuns{1}.pred.rf;
        data.params.x0.both            = params.analysis.spatial.(hemi).x0;
        data.params.y0.both            = params.analysis.spatial.(hemi).y0;
        data.params.effectiveSize.both = params.analysis.spatial.(hemi).effectiveSize;
        data.params.exponent.both      = params.analysis.spatial.(hemi).exponent;
        data.params.varexpl.both       = params.analysis.spatial.(hemi).varexpl;
        if isfield(data.predRuns,hemi)
            data.dataRuns.both = mean(dataOrig.dataRuns.(hemi).(['square' mat2str(sqr)]),3,'omitnan');
            data.detrendedData.both = mean(dataOrig.detrendedData.(hemi).(['square' mat2str(sqr)]),2,'omitnan');
            data.meanDetrendedData = mean(dataOrig.meanDetrendedData.(hemi).(['square' mat2str(sqr)]),2,'omitnan');
            data.meanDetrendedDataOdd = mean(dataOrig.meanDetrendedDataOdd.(hemi).(['square' mat2str(sqr)]),2,'omitnan');
            data.meanDetrendedDataEven = mean(dataOrig.meanDetrendedDataEven.(hemi).(['square' mat2str(sqr)]),2,'omitnan');
        else
            data.dataRuns.both        = NaN(size(data.predRuns.both));
            data.detrendedData.both   = NaN(size(data.predRuns.both));
            data.meanDetrendedData    = NaN(size(data.predRuns.both));
            data.meanDetrendedDataOdd = NaN(size(data.predRuns.both));
            data.meanDetrendedDataEven = NaN(size(data.predRuns.both));
        end
        % Check if nr of voxels in data match predictions and clear some memoroy
        assert(size(data.meanDetrendedData,2)==size(data.predRuns.(hemi),2))
        assert(size(data.meanDetrendedData,2)==size(data.rf.(hemi),3))
        clear predRunBothHemi;
        
        % Check flags and if data is actually present
        fn = fieldnames(data.predRuns);
        if length(fn)==1 && isempty(alignMe.pRFs.selected.(fn{1}))
            break
        elseif length(fn)==2 && isempty(alignMe.pRFs.selected.(fn{1})) && isempty(alignMe.pRFs.selected.(fn{2}))
            break
        elseif all(isnan(squeeze(data.meanDetrendedData)))
            break
        else
            
            % Permute from time x voxels x channels x runs to ...
            % time x runs x voxels x channels, then concatenate unique runs
            data.predictionCatRuns = permute(data.predRuns.both, [1 4 2 3]);
            data.predictionCatRuns = reshape(data.predictionCatRuns, ...
                size(data.predictionCatRuns,1)*size(data.predictionCatRuns,2), ...
                size(data.predictionCatRuns,3),size(data.predictionCatRuns,4));
            
            assert(size(data.meanDetrendedData,1)==size(data.predictionCatRuns,1))

            %% Fit model prediction to data at single voxel level
            params.analysis.regressionType = 'OLS';%'fracridge';
            
            X = data.predictionCatRuns;
            
            % Normalize max height across concatenated runs, for each channel
            X_norm = NaN(size(X));
            for chan = 1:size(X,3)
                X_norm(:,:,chan) = normMax(X(:,:,chan));
            end
            
            X = X_norm;  clear X_norm
            
            % Fit sum of equal weighted Sustained and Transient channels if requested
            if sumSTFlag
                X = normMax(bsxfun(@(x) plus(x,3), X));
            end
            % Weight channels if requested. The variable weightChannels should
            % be a 1xN vector with N equal to number of channels and with
            % values between 0 and 1. By forcing weights, we can sum channels
            % and then fit one regressor to avoid poor fits due to collinearity.
            if ~isempty(weightChannels)
                X = X*weightChannels';
                X = normMax(bsxfun(@(x) plus(x,3)));
            end
            
            % Add offset if requested
            if addOffsetFlag
                X = cat(3,X,ones(size(X,1),size(X,2),1));
            end
            
            [lm, sumChannelPrediction] = ...
                fitModelPredictionToDataWrapper(data.meanDetrendedData, X, ...
                'regressionType', params.analysis.regressionType);
            
            lmcell = struct2cell(lm);
            fnBeta = find(strcmp(fieldnames(lm),'betas'));
            fnR2   = find(strcmp(fieldnames(lm),'R2'));
            
            B_full  = squeeze(cell2mat(lmcell(fnBeta,:,:)))';
            R2_full = squeeze(cell2mat(lmcell(fnR2,:,:)))';
            
            % Generate also a response without offset
            if addOffsetFlag
                offsetVals = reshape(B_full(:,end),[1 size(B_full(:,end))]);
                sumChannelPredictionNoOffset = sumChannelPrediction - offsetVals.*X(:,:,end);
            end
            
            numTimePoints = size(data.meanDetrendedData,1);
            numVoxels     = size(data.meanDetrendedData,2);
            numChannels   = size(data.predictionCatRuns,3);
            
            if doCrossValidation
                Y_crossval = {data.meanDetrendedDataOdd, data.meanDetrendedDataEven};
                
                if strcmp(params.analysis.regressionType,'fracridge')
                    bestAlpha = lm(1).bestAlpha;
                else, bestAlpha = []; end
                
                % First fit split half runs
                for cx = [1,2]
                    [lmCV_tmp{cx}, sumChannelPredictionCV_tmp{cx}] = ...
                        fitModelPredictionToDataWrapper(Y_crossval{cx}, X, 'alpha',bestAlpha, ...
                        'regressionType', params.analysis.regressionType);
                    
                    % Accumulate betas
                    lmcell_tmp = struct2cell(lmCV_tmp{cx});
                    B_crossval_tmp(:,:,cx)  = cell2mat(squeeze(lmcell_tmp(fnBeta,:,:)));
                    
                    % Generate also best predicted response without offset
                    if addOffsetFlag
                        oddOffset = B_crossval_tmp(:,end,cx);
                        offsetVals = reshape(oddOffset,[1 size(oddOffset)]);
                        sumChannelPredictionCVNoOffset_tmp{cx} = sumChannelPredictionCV_tmp{cx} - offsetVals.*X(:,:,end);
                        
                        Y_crossval_noOffset{cx} = Y_crossval{cx} - offsetVals;
                    end
                end
                clear lmCV_tmp B_crossval_tmp sumChannelPredictionCV_tmp
                
                % Now refit split half runs with other halve
                for cx = [1,2]
                    
                    if addOffsetFlag
                        Y = Y_crossval_noOffset;
                    else
                        Y = Y_crossval;
                    end
                    
                    % Fit data again to check if B0 is close to 0
                    [lmCV{cx}, sumChannelPredictionCV{cx}] = ...
                        fitModelPredictionToDataWrapper(Y{cx}, X, 'alpha',bestAlpha, ...
                        'regressionType', params.analysis.regressionType);
                    
                    % Accumulate betas
                    lmcell_cv = struct2cell(lmCV{cx});
                    B_crossval(:,:,cx)  = cell2mat(squeeze(lmcell_cv(fnBeta,:,:)));
                    
                    % Get idx for other dataset to fit
                    fitIdx = setdiff([1,2],cx);
                    
                    % We have to per voxel, otherwise we rectify all voxels in
                    % denominator or CoD calculation
                    for n = 1:numVoxels
                        R2_crossval(cx,n) = computeCoD(Y{fitIdx}(:,n),sumChannelPredictionCV{cx}(:,n));
                    end
                    
                    % Generate also a response without offset
                    if addOffsetFlag
                        offsetVals = reshape(B_crossval(:,end,cx),[1 size(B_crossval(:,end,cx))]);
                        sumChannelPredictionCVNoOffset{cx} = sumChannelPredictionCV{cx} - offsetVals.*X(:,:,end);
                    end
                end
                clear sumChannelPredictionCVNoOffset_tmp offsetVals
                
                B_crossval_mean = mean(B_crossval,3,'omitnan');
                R2_crossval_mean = mean(R2_crossval,1,'omitnan');
                sumChannelPredictionCrossvalMean = mean(cat(3,sumChannelPredictionCV{1},...
                    sumChannelPredictionCV{2}),3,'omitnan');
                Y_crossval_mean = mean(cat(3,Y{1},Y{2}),3,'omitnan');
                if addOffsetFlag
                    sumChannelPredictionNoOffsetCrossvalMean = mean(cat(3,sumChannelPredictionCVNoOffset{1},...
                        sumChannelPredictionCVNoOffset{2}),3,'omitnan');
                end
            end
            
            
            %% Plot run
            if size(stimRunUnique,1)<size(stimRunUnique,2)
                stimRunUnique = stimRunUnique';
            end
            if verbose
                % Select top 30 voxels based on R2
                if doCrossValidation
                    dataToSort = R2_crossval_mean; 
                else, dataToSort = R2_full; end
 
                [~,idx] = sort(R2_crossval_mean, 'descend');
                if numVoxels < 30; numVoxToPlot = idx;
                else, numVoxToPlot = idx(1:30); end
                
                if saveFigs
                    folderName = sprintf('%s_%s_%s_%s_AverageRun_addOffsetFlag%d_cval%d', ...
                        rois{roi},roiType,spatialModel, temporalModel, addOffsetFlag, doCrossValidation);
                        folderName = sprintf('%s_%s',hemi, folderName);
                    subFolderDir = fullfile(saveFigDir, 'figs', folderName); 
                    if ~exist(subFolderDir,'dir'), mkdir(subFolderDir); end
                end
                
                fH = figure(101); clf; set(gcf, 'position', [55,440,1994,522]);
                stimToPlot = normMax(stimRunUnique(:))';
                
                for ii = numVoxToPlot
                    
                    if doCrossValidation
                        dataToPlot = Y_crossval_mean(:,ii);
                        modelPredToPlot = sumChannelPredictionCrossvalMean(:,ii);
                        betaToPrint = B_crossval_mean(ii,:);
                        R2ToPrint = R2_crossval_mean(ii);
                    else
                        dataToPlot = data.meanDetrendedData(:,ii);
                        modelPredToPlot = sumChannelPrediction(:,ii);
                        betaToPrint = B_full(ii,:);
                        R2ToPrint = R2_full(ii);
                    end
                    ttl = sprintf('Voxel %d %s %s: [x,y,sz]=[%1.1f,%1.1f,%1.2f]. %s %s: beta %s - R2 %2.1f%%', ...
                        ii, rois{roi}, roiType, data.params.x0.both(ii), ...
                        data.params.y0.both(ii), data.params.effectiveSize.both(ii), ...
                        spatialModel, temporalModel, mat2str(betaToPrint,2), ...
                        R2ToPrint*100);
                    
                    fH = plotFullRunTSwithModelFit(fH, dataToPlot, modelPredToPlot, ...
                        stimToPlot,ttl);
                    
                    if saveFigs
                        currFName = sprintf('voxel%d_fullAveData_vs_model_sqr%d',ii,sqr);
                        print(fH, fullfile(subFolderDir,currFName), '-dpng')
                        % print(fH, fullfile(subFolderDir,currFName), '-depsc')
                    end
                end
                
                
            end
            
            %% Chop run data and predictions into trials,
            % create a new data type with average run 1,2,3 concatenated
            cd(fullfile(simSeqSessionDir));
            setVAnatomyPath('./3DAnatomy/t1.nii.gz');
            postFix = sprintf('%d',1:max(stimRun));
            newDt = sprintf('Average_Run%s',postFix);
                        
            load('mrSESSION.mat')            
            for dd = 1:length(dataTYPES)
                if strcmp(dataTYPES(dd).name,newDt)
                    dt_idx = dd;
                    break
                end
            end
            if ~exist('dt_idx','var') || isempty(dt_idx)
                error('[%s]: mrSESSION dataTYPES does not contain datatype! Check mrSESSION.mat',mfilename)
            end
            
            if dataTYPES(dt_idx).scanParams.nFrames~=size(data.meanDetrendedData,1)
                dataTYPES(dt_idx).scanParams.nFrames = size(data.meanDetrendedData,1);
                saveSession;
            end
            
            tmp = dataTYPES(dt_idx).scanParams.parfile(1:end-4);
            [~,dtmp] = fileparts(tmp);
            if isempty(regexp(dtmp, '_corrected2','once'))
                dataTYPES(dt_idx).scanParams.parfile = fullfile([tmp '_corrected2.par']);
                saveSession;
            end
            
            hvol    = initHiddenGray;
            hvol    = viewSet(hvol, 'curdt', newDt);
            trials = er_concatParfiles(hvol, 1);
            
            if numel(stimRun)>1
                er_data(1).params.timeWindow = repmat(totalBlockTRs, length(trials.condNums), 1);
                er_data(1).params.bslPeriod  = repmat(baselineBlockTRs, length(trials.condNums), 1);
                er_data(1).params.normBsl    = 0;
                er_data(1).params.parfiles   = trials.parfiles;

            else % Pilot 1 (only one type of run)
                postFix = '1';
                er_data(1).params.timeWindow = {-2:14, -2:14, -2:14}; % BLANK, SEQ & SIM
                er_data(1).params.bslPeriod = repmat({-2:1:0}, 3, 1);
                er_data(1).params.normBsl = 0;
                er_data(1).trials.onsetFrames    = er_data(1).trials.onsetFrames - (removedTRsAtStart+extraTRsToRemoveAtBeginning);
                er_data(1).trials.onsetFrames(1) = 1;
                er_data(1).trials.onsetSecs      = er_data(1).trials.onsetFrames - (removedTRsAtStart+extraTRsToRemoveAtBeginning);
                er_data(1).trials.onsetSecs(1)   = 0;
                er_data(1).trials.framesPerRun   = er_data(1).trials.framesPerRun - (removedTRsAtStart+extraTRsToRemoveAtBeginning);
                trials = er_data(1).trials;
            end
            
            %% Convert struct into a table with voxel by condition, separate for data and model predictions

            % DATA: NO Cross-Validation
            % voxData: data.meanDetrendedData => allDataTrials meanDataTrials semDataTrials
            [TD1,allDataTrials,~,~,~, catTrialsData] = ...
                simseq_getStimBlocksFromDataRuns(data.meanDetrendedData,trials,er_data(1).params);  %#ok<ASGLU>
            TD1.Properties.VariableNames = {'Condition','DataTS','DataTSerror'};
            
            [T1odd,allDataTrialsOdd] = ...
                simseq_getStimBlocksFromDataRuns(data.meanDetrendedDataOdd,trials,er_data(1).params); %#ok<ASGLU>
            T1odd.Properties.VariableNames = {'Condition','DataTSOdd','DataTSOdderror'};
            
            [T1even,allDataTrialsEven] = ...
                simseq_getStimBlocksFromDataRuns(data.meanDetrendedDataEven,trials,er_data(1).params); %#ok<ASGLU>
            T1even.Properties.VariableNames = {'Condition','DataTSEven','DataTSEvenerror'};
            
            % Update data table (since we removed some voxels because they
            % were eliminated by modelfit
            T.DataTS          = TD1.DataTS;
            T.DataTSerror     = TD1.DataTSerror;
            T.DataTSOdd       = T1odd.DataTSOdd;
            T.DataTSOdderror  = T1odd.DataTSOdderror;
            T.DataTSEven      = T1even.DataTSEven;
            T.DataTSEvenerror = T1even.DataTSEvenerror;

            if doCrossValidation
                % DATA: WITH Cross-Validation
                % voxDataCrossvalMean: Y_crossval_mean => allDataTrialsCrossvalMean meanDataTrialsCrossvalMean semDataTrialsCrossvalMean
                [TD_CV,allDataTrialsCrossvalMean,~,~,catTrialsDataCV] = ...
                    simseq_getStimBlocksFromDataRuns(Y_crossval_mean,trials,er_data(1).params); 
                    TD_CV.Properties.VariableNames = {'Condition','DataTSCV','DataTSCVError'};
                
                    % MODEL: WITH Cross-Validation but NO removed offset (since we didn't add it to the model)   
                    if ~addOffsetFlag
                        % voxModelCrossvalMean: sumChannelPredictionCrossvalMean => allModelTrialsCrossvalMean meanModelTrialsCrossvalMean semModelTrialsCrossvalMean
                        [TM_CV,allModelTrialsCrossvalMean,~,~] = ...
                                simseq_getStimBlocksFromDataRuns(sumChannelPredictionCrossvalMean,trials,er_data(1).params); 
                            TM_CV.Properties.VariableNames = {'Condition','ModelTSCVMn','ModelTSCVError'};
                        
                    % MODEL: WITH Cross-Validation WITH removed offset (as we added it to the model)   
                    elseif addOffsetFlag
                        % voxModelCrossvalMeanNoOffset: sumChannelPredictionNoOffsetCrossvalMean => allModelTrialsCrossvalMeanNoOffset meanModelTrialsCrossvalMeanNoOffset semModelTrialsCrossvalMeanNoOffset
                        [TM_CV_noOffset,allModelTrialsCrossvalMeanNoOffset,~,~] = ...
                            simseq_getStimBlocksFromDataRuns(sumChannelPredictionNoOffsetCrossvalMean,trials,er_data(1).params);
                        TM_CV_noOffset.Properties.VariableNames = {'Condition','ModelTSCVMn','ModelTSCVError'};

                    end
           
            elseif ~doCrossValidation
                % MODEL: NO cross-Validation WITH removed offset (as we added it to the model)       
                if addOffsetFlag
                    % voxModelNoOffset: sumChannelPredictionNoOffset = allModelNoOffsetTrials meanModelNoOffsetTrials semModelNoOffsetTrials
                    [TM_noOffset,allModelNoOffsetTrials,~,~] = ...
                        simseq_getStimBlocksFromDataRuns(sumChannelPredictionNoOffset,trials,er_data(1).params);  %#ok<ASGLU>
                        TM_noOffset.Properties.VariableNames = {'Condition','ModelTS','ModelTSerror'};
                
                % MODEL: NO cross-Validation and NO removed offset (since we didn't add it to the model)     
                elseif ~addOffsetFlag
                    % voxModel: sumChannelPrediction => allModelTrials meanModelTrials semModelTrials
                    [TM,allModelTrials,~,~,~, ~] = ...
                        simseq_getStimBlocksFromDataRuns(sumChannelPrediction,trials,er_data(1).params); 
                        TM.Properties.VariableNames = {'Condition','ModelTS','ModelTSerror'};
                end
            end
              
           % Concatenate tables
           if ~addOffsetFlag && ~doCrossValidation
              Tmodel = TM;
            elseif addOffsetFlag && ~doCrossValidation
              Tmodel = TM_noOffset;
            elseif ~addOffsetFlag && doCrossValidation
              Tmodel = [TD_CV, TM, TM_CV];
            elseif addOffsetFlag && doCrossValidation
              Tmodel = [TD_CV, TM_noOffset, TM_CV_noOffset];
           end
                
           if strcmp(hemi,'lh'), currHemi = 0; else, currHemi =1; end
           TmodelThisSquare = table(Tmodel{:,1},currHemi.*ones(numConditions,1),sqr.*ones(numConditions,1), ...
               Tmodel{:,2},Tmodel{:,3},Tmodel{:,4},Tmodel{:,5},Tmodel{:,6},Tmodel{:,7},...
               'VariableNames', {Tsquare.Properties.VariableNames{1},...
               Tsquare.Properties.VariableNames{2},...
               Tsquare.Properties.VariableNames{3},...
               'DataTSCV',...
               'DataTSCVError',...
               'ModelTS',...
               'ModelTSerror',...
               'ModelTSCVMn',...
               'ModelTSCVError'});
           TmodelFull = [TmodelFull; TmodelThisSquare];
            
           
           % Recompute Split half reliability if we removed data voxels..
           % .. across runs
           if ~isempty(data.detrendedData.(hemi).(['square' mat2str(sqr)]))
               if ~isfield(data.detrendedData.(hemi),'run'),
                   data.detrendedData.(hemi).run = []; end
               
               data.detrendedData.(hemi).run.(['square' mat2str(sqr)]) = simseq_getSplitHalfDataReliability( ...
                   permute(data.detrendedData.(hemi).(['square' mat2str(sqr)]),[1,3,2]), 'runFlag', true, ...
                   'splitRuns', pths.runOrder);
               
               if doCrossValidation
                   if ~isfield(data.detrendedData.(hemi),'runCV'),
                       data.detrendedData.(hemi).runCV = [];
                   end
                   if ~isfield(data.detrendedData.(hemi),'trialCV'),
                       data.detrendedData.(hemi).trialCV = [];
                   end
                   data.splitHalfRel.(hemi).runCV.(['square' mat2str(sqr)]) = simseq_getSplitHalfDataReliability( ...
                       permute(data.dataRuns.detrendedStimRuns.(hemi).(['square' mat2str(sqr)]),[1,3,2]), ...
                       'runFlag',true, 'splitRuns', []);
                   
                   data.splitHalfRel.(hemi).trialCV.(['square' mat2str(sqr)]) = simseq_getSplitHalfDataReliability( ...
                       catTrialsDataCV, 'runFlag',false,'splitRuns',[]);
               end
               
           else
               data.detrendedData.(hemi).run.(['square' mat2str(sqr)]) = [];
               data.splitHalfRel.(hemi).runCV.(['square' mat2str(sqr)]) = [];
               data.splitHalfRel.(hemi).trialCV.(['square' mat2str(sqr)]) = [];
           end
           %% Get R2 for model prediction at trial level, per condition
           for cond = 1:numConditions
               currData = cell2mat(Tsquare.DataTS(Tsquare.Condition==cond & Tsquare.Hemi_lh1_rh2==currHemi & Tsquare.Square==sqr));
               currModel = cell2mat(TmodelFull.ModelTS(TmodelFull.Condition==cond & TmodelFull.Hemi_lh1_rh2==currHemi & TmodelFull.Square==sqr));
               if doCrossValidation
                   currDataCV = cell2mat(TmodelFull.DataTSCV(TmodelFull.Condition==cond & TmodelFull.Hemi_lh1_rh2==currHemi & TmodelFull.Square==sqr));
                   currModelCV = cell2mat(TmodelFull.ModelTSCVMn(TmodelFull.Condition==cond & TmodelFull.Hemi_lh1_rh2==currHemi & TmodelFull.Square==sqr));
               end
               
               for vv = 1:size(currData,2)
                   if sum(~isnan(currData(:,vv)) & ~isnan(currModel))==size(currData,1)
                       R2_full_trial{cond,sqr}(:,vv) = computeCoD(currData(:,vv), currModel);
                   else
                       R2_full_trial{cond,sqr}(:,vv) = NaN;
                   end
               end
               if doCrossValidation
                   if sum(~isnan(currDataCV) & ~isnan(currModelCV))==size(currDataCV,1)
                       R2_crossval_trial(cond,sqr) = computeCoD(currDataCV, currModelCV);
                   end
               end
               
               
           end
           clear currData currModel
           
           % Take mean and SEM across voxels
           if  size(R2_full_trial,2)>=sqr
               meanR2_full_trial(cond,sqr) = mean(R2_full_trial{cond,sqr},2,'omitnan');  %#ok<AGROW>
               seR2_full_trial = std(R2_full_trial{cond,sqr},[],2,'omitnan')./numVoxels; %#ok<NASGU>
               if doCrossValidation
                   if ~isnan(R2_crossval_trial(:,sqr))
                       meanR2_crossval_trial(:,sqr) = R2_crossval_trial(:,sqr); %#ok<NASGU>
                       seR2_crossval_trial(:,sqr) = NaN; %#ok<NASGU>
                   else
                       meanR2_crossval_trial(:,sqr) = NaN;
                       seR2_crossval_trial(:,sqr) = NaN;
                   end
               end
           else
               meanR2_full_trial(cond,sqr) = NaN;
               seR2_full_trial(cond,sqr) = NaN;
               
           end
           
           R2_full_square(sqr) = R2_full;
           B_full_square(:,sqr) = B_full;
           sumChannelPrediction_square(:,sqr) = sumChannelPrediction;
           params_square.(['square' mat2str(sqr)]) = params;
           data_square.(['square' mat2str(sqr)]) = data;
           allModelTrials_square.(['square' mat2str(sqr)]) = allModelTrials;
           
           if doCrossValidation
               B_crossval_square(sqr,:,:) = B_crossval;
               B_crossval_mean_square(sqr,:) = B_crossval_mean;
               R2_crossval_square(:,sqr)  = R2_crossval;
               R2_crossval_mean_square(sqr,:) = R2_crossval_mean;
               allModelTrialsCrossvalMean_square{sqr} = allModelTrialsCrossvalMean;
               allDataTrialsCrossvalMean_square{sqr} = allDataTrialsCrossvalMean;
               sumChannelPredictionCrossvalMean_square(:,sqr) = sumChannelPredictionCrossvalMean;
               sumChannelPredictionCV_square{sqr} = sumChannelPredictionCV;
               Y_crossval_mean_square(:,sqr) = Y_crossval_mean;
               Y_crossval_noOffset_square{sqr} = Y_crossval_noOffset;
               lmCV_square{sqr} = lmCV;
           end
           if addOffsetFlag
               allModelTrialsCrossvalMeanNoOffset_square{sqr,:} = allModelTrialsCrossvalMeanNoOffset;
               sumChannelPredictionNoOffsetCrossvalMean_square(:,sqr) = sumChannelPredictionNoOffsetCrossvalMean;
               sumChannelPredictionCVNoOffset_square{sqr,:} = sumChannelPredictionCVNoOffset;
           end
        end
    end
    
    
    %% save betas
    if saveBetas
            saveList = {'Tsquare','TmodelFull','B_full','R2_full','sumChannelPrediction', ...
                'allModelTrials','data','nrTrials','er_data',...
                'R2_full_trial', 'meanR2_full_trial','seR2_full_trial',...
                'stimRunUnique','alignMe','params'};
            
            R2_full = R2_full_square;
            B_full  = B_full_square;
            sumChannelPrediction = sumChannelPrediction_square;
            params = params_square;
            data   = data_square;
            allModelTrials = allModelTrials_square;
            
            
            if doCrossValidation
                B_crossval = B_crossval_square;
                B_crossval_mean = B_crossval_mean_square;
                R2_crossval = R2_crossval_square;
                R2_crossval_mean = R2_crossval_mean_square;
                allModelTrialsCrossvalMean = allModelTrialsCrossvalMean_square;
                allDataTrialsCrossvalMean = allDataTrialsCrossvalMean_square;
                sumChannelPredictionCrossvalMean = sumChannelPredictionCrossvalMean_square;
                sumChannelPredictionCV = sumChannelPredictionCV_square;
                Y_crossval_mean = Y_crossval_mean_square;
                Y_crossval_noOffset = Y_crossval_noOffset_square;
                lmCV = lmCV_square;
            end
            if addOffsetFlag
                allModelTrialsCrossvalMeanNoOffset = allModelTrialsCrossvalMeanNoOffset_square;
                sumChannelPredictionNoOffsetCrossvalMean = sumChannelPredictionNoOffsetCrossvalMean_square;
                sumChannelPredictionCVNoOffset = sumChannelPredictionCVNoOffset_square;
            end
            
        if addOffsetFlag
            saveList = {saveList{:},'sumChannelPredictionNoOffset','allModelNoOffsetTrials'};
        end
        
        if doCrossValidation
            saveList = {saveList{:},'B_crossval','B_crossval_mean', ...
                'R2_crossval','R2_crossval_mean', ...
                'R2_crossval_trial',...
                'meanR2_crossval_trial','seR2_crossval_trial',...
                'allModelTrialsCrossvalMean', ...
                'allDataTrialsCrossvalMean',...
                'sumChannelPredictionCrossvalMean', ...
                'sumChannelPredictionCV','Y_crossval_mean','Y_crossval_noOffset', ...
                'lmCV'};
            
            if addOffsetFlag
                saveList = {saveList{:},'allModelTrialsCrossvalMeanNoOffset',...
                    'sumChannelPredictionNoOffsetCrossvalMean',...
                    'sumChannelPredictionCVNoOffset'};
            end
        end
        
        if ~exist(fullfile(saveBetaFolder,subFolderName2,'betaR2_XValFit_OLS'),'dir')
            mkdir(fullfile(saveBetaFolder,subFolderName2,'betaR2_XValFit_OLS'));
        end

        loadFile = fullfile(saveBetaFolder,subFolderName2,'betaR2_XValFit_OLS',...
            sprintf('betaValsFull_%s_%s_%s_%s_%s_run%s_addOffsetFlag%d_%s_bslCorrected_avgSquare.mat', ...
            hemi,rois{roi},roiType,spatialModel, temporalModel, postFix, addOffsetFlag, params.square1.analysis.regressionType));
        
        save(loadFile, saveList{:},'-v7.3');
    end
    
    clear B_crossval R2_crossval
end

return









