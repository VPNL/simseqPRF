function [] = simseq_fitTemporalModelPredictionsToPreprocData(subjnr,projectDir,spatialModel,temporalModel,varargin)
%   simseq_fitTemporalModelPredictionsToPreprocData(subjnr, projectDir, spatialModel, temporalModel)
%
% Function to fit generated spatiotemporal model predictions using pRF of
% single vertex from toonotopy experiment to measured fMRI data.
%
% NOTE 1:
% Data must be preprocessed with s_runDataPreproc.m
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
p.addParameter('useSearchFit',false,@islogical);
p.addParameter('loadDataFolder',[], @ischar);
p.addParameter('saveBetaFolder',[], @ischar);
p.addParameter('searchAlgorithm','bads',@ischar);
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
cd(fullfile(projectDir,'data','simseq'))
pths             = getSubjectPaths(projectDir,subjnr);
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
if strcmp(hemi, 'both'); hm = {'lh','rh'}; else hm = {hemi}; end

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
        sprintf('preprocData_%s_%s_%s_run%s.mat', ...
        rois{roi},roiType,'cssFit', postFix));
    load(loadFile, 'T','allDataTrials','data','trials',...
        'nrTrials','er_data','stimMS','alignMe');
    
    
    % Define totalScans
    if isfield(data.dataRuns,'rh')
        if ~isempty(data.dataRuns.rh)
            totalScans = size(data.dataRuns.rh,2);
        end
    else, totalScans = size(data.dataRuns.lh,2); end
    
    %% Get predictions
    
    % Preallocate space
    predRunBothHemi = []; stimRunUnique = [];
    predAllRuns = cell(1,length(stimRun));
    
    for run = 1:length(stimRun)
        if ~isempty(subFolder)
            subFolderName2 = sprintf('%s%d_v%d',subFolder, run, pths.expversionNr);
        end
        
        % Load prediction file
        if ~isnan(useFixedExponent)
            modelPredSubFolder = 'modelPredictions_gridfit';
            predFile = fullfile(saveBetaFolder,subFolderName2, modelPredSubFolder,...
                sprintf('modelPredictions_%s_%s_%s_%s_run%d_exp%1.2f.mat', ...
                rois{roi},roiType,spatialModel, temporalModel, run, useFixedExponent));
        elseif useSearchFit
             modelPredSubFolder = sprintf('modelPredictions_%s_searchfit_restricted',searchAlgorithm);
             predFile = fullfile(saveBetaFolder,subFolderName2, modelPredSubFolder,...
                sprintf('modelPredictions_%s_%s_%s_%s_run%d_bestsFitExp.mat', ...
                rois{roi},roiType,spatialModel, temporalModel, run));
        elseif isnan(useFixedExponent)
            predFile = fullfile(saveBetaFolder,subFolderName2, ...
                sprintf('modelPredictions_%s_%s_%s_%s_run%d.mat', ...
                rois{roi},roiType,spatialModel, temporalModel, run));
        end
        predAllRuns{run} = load(predFile);
        
        % Concatenate predictions
        predRunBothHemi = cat(4,predRunBothHemi, predAllRuns{run}.pred.predHRF);
        
        % Get stimulus file, so we can plot it later
        maxStimMS = floor(size(predRunBothHemi,1)-1)*fs;
        
        % .. sum
        if size(predAllRuns{run}.pred.params.stim,2)==1
           stim = sum(predAllRuns{run}.pred.params.stim.images,1);
        else
           stim = sum(predAllRuns{run}.pred.params.stim(run).images,1);
        end
        % .. normalize to make plotting easier
        stimRunUnique(run,:) = stim./sum(stim(:));
    end
    
    % Store and update params
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
    if ~useSearchFit
        predRunBothHemi(1:removedTRsAtStart,:,:,:) = [];
    end
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
    if size(stimRunUnique,1)<size(stimRunUnique,2)
        stimRunUnique = stimRunUnique'; end
        
    %% Check voxel alignment of datasets
    if ~useSearchFit
        coords.dataRuns = alignMe.coords.dataRuns;
        alignMe = simseq_checkVoxelCorrespondance(coords, selectedpRFs);
        numOrigVoxLH  = size(alignMe.coords.dataRuns.lh,2);
        numOrigpRFsLH = size(alignMe.coords.pRFs.lh,2);

        % remove voxels from dataRuns and data gray coordinates
        if isfield(alignMe.pRFs,'selected')
            if isfield(alignMe.pRFs.selected,'lh') && ~isempty(alignMe.pRFs.selected.lh)
                data.predRuns.lh = predRunBothHemi(:,alignMe.pRFs.selected.lh,:,:);
                data.rf.lh       = predAllRuns{1}.pred.rf(:,:,alignMe.pRFs.selected.lh);
            end
            if isfield(alignMe.pRFs.selected,'rh') && ~isempty(alignMe.pRFs.selected.rh)
                data.predRuns.rh = predRunBothHemi(:,numOrigpRFsLH+alignMe.pRFs.selected.rh,:,:);
                data.rf.rh       = predAllRuns{1}.pred.rf(:,:,numOrigpRFsLH+alignMe.pRFs.selected.rh);
            end
        end
    
    % remove voxels from dataRuns and data gray coordinates
    if isfield(alignMe.dataRuns, 'selected')
        if isfield(data.dataRuns, 'lh') && ~isempty(data.dataRuns.lh)
            if size(data.dataRuns.lh,3) ~= size(data.predRuns.lh,2)
                if isfield(alignMe.dataRuns.selected,'lh') && ~isempty(alignMe.dataRuns.selected.lh)
                    data.dataRuns.lh = data.dataRuns.lh(:,:,alignMe.dataRuns.selected.lh);
                end
            end
        end
        if isfield(data.dataRuns, 'rh') && ~isempty(data.dataRuns.rh)
            if size(data.dataRuns.rh,3) ~= size(data.predRuns.rh,2)
                if isfield(alignMe.dataRuns.selected,'rh') && ~isempty(alignMe.dataRuns.selected.rh)
                    data.dataRuns.rh = data.dataRuns.rh(:,:,alignMe.dataRuns.selected.rh);
                end
            end
        end
        
        if isfield(data.dataRuns,'lh') && ~isempty(data.dataRuns.lh)
            numVoxLH  = size(data.dataRuns.lh,3);
        else
            numVoxLH = 0;
        end
        if isfield(data.dataRuns,'rh') && ~isempty(data.dataRuns.rh)
            numVoxRH  = size(data.dataRuns.rh,3);
        else
            numVoxRH = 0;
        end
        bothVoxIdx = [alignMe.dataRuns.selected.lh; numOrigVoxLH+alignMe.dataRuns.selected.rh];
        isequal(sum([numVoxLH,numVoxRH]), length(bothVoxIdx));
        
        if isfield(data.dataRuns, 'both') && strcmp(hemi,'both')
            data.dataRuns.both         = data.dataRuns.both(:,:,bothVoxIdx);
            data.detrendedData.both    = data.detrendedData.both(1:size(data.dataRuns.both,1),bothVoxIdx,:);
            data.meanDetrendedData     = data.meanDetrendedData(:,bothVoxIdx);
            data.meanDetrendedDataOdd  = data.meanDetrendedDataOdd(:,bothVoxIdx);
            data.meanDetrendedDataEven = data.meanDetrendedDataEven(:,bothVoxIdx);
            
            % Update Table with block data as well
            T_orig = T;
            
            for jj = 1:size(T,1)
                for ll = 2:size(T,2)
                    T{jj,ll}{1} = T{jj,ll}{1}(:,bothVoxIdx);
                end
            end
        else
            T_orig = [];
        end
    end
    
    %% Concatenate left and right hemi into both hemis
    % Run data = time by runs by voxels (in 1s TRs)
    % PRF data = x by y by voxels (in deg visual angle)
    if isfield(data.predRuns,'lh') &&  isfield(data.predRuns,'rh')
        data.predRuns.both             = cat(2,data.predRuns.lh,data.predRuns.rh);
        data.rf.both                   = cat(3,data.rf.lh,data.rf.rh);
        data.params.x0.both            = cat(2,params.analysis.spatial.lh.x0(alignMe.pRFs.selected.lh), ...
            params.analysis.spatial.rh.x0(alignMe.pRFs.selected.rh));
        data.params.y0.both            = cat(2,params.analysis.spatial.lh.y0(alignMe.pRFs.selected.lh), ...
            params.analysis.spatial.rh.y0(alignMe.pRFs.selected.rh));
        data.params.effectiveSize.both = cat(2,params.analysis.spatial.lh.effectiveSize(alignMe.pRFs.selected.lh), ...
            params.analysis.spatial.rh.effectiveSize(alignMe.pRFs.selected.rh));
        data.params.exponent.both      = cat(2,params.analysis.spatial.lh.exponent(alignMe.pRFs.selected.lh), ...
            params.analysis.spatial.rh.exponent(alignMe.pRFs.selected.rh));
        if ~isempty(useFixedExponent) && ~isnan(useFixedExponent)
            data.params.exponent_temporal.both  = cat(2,predAllRuns{1}.pred.params.analysis.temporal.param.exponent(alignMe.pRFs.selected.lh), ...
                predAllRuns{1}.pred.params.analysis.temporal.param.exponent(numOrigpRFsLH+alignMe.pRFs.selected.rh));
        end
        data.params.varexpl.both       = cat(2,params.analysis.spatial.lh.varexpl(alignMe.pRFs.selected.lh), ...
            params.analysis.spatial.rh.varexpl(alignMe.pRFs.selected.rh));
    elseif isfield(data.predRuns,'lh')
        data.predRuns.both             = data.predRuns.lh;
        data.rf.both                   = data.rf.lh;
        data.params.x0.both            = params.analysis.spatial.lh.x0(alignMe.pRFs.selected.lh);
        data.params.y0.both            = params.analysis.spatial.lh.y0(alignMe.pRFs.selected.lh);
        data.params.effectiveSize.both = params.analysis.spatial.lh.effectiveSize(alignMe.pRFs.selected.lh);
        data.params.exponent.both      = params.analysis.spatial.lh.exponent(alignMe.pRFs.selected.lh);
         if ~isempty(useFixedExponent) && ~isnan(useFixedExponent)
            data.params.exponent_temporal.both  = predAllRuns{1}.pred.params.analysis.temporal.param.exponent(alignMe.pRFs.selected.lh);
        end
        data.params.varexpl.both       = params.analysis.spatial.lh.varexpl(alignMe.pRFs.selected.lh);
    elseif isfield(data.predRuns,'rh')
        data.predRuns.both             = data.predRuns.rh;
        data.rf.both                   = data.rf.rh;
        data.params.x0.both            = params.analysis.spatial.rh.x0(alignMe.pRFs.selected.rh);
        data.params.y0.both            = params.analysis.spatial.rh.y0(alignMe.pRFs.selected.rh);
        data.params.effectiveSize.both = params.analysis.spatial.rh.effectiveSize(alignMe.pRFs.selected.rh);
        data.params.exponent.both      = params.analysis.spatial.rh.exponent(alignMe.pRFs.selected.rh);
        if ~isempty(useFixedExponent) && ~isnan(useFixedExponent)
            data.params.exponent_temporal.both  = predAllRuns{1}.pred.params.analysis.temporal.param.exponent(numOrigpRFsLH+alignMe.pRFs.selected.rh);
        end
        data.params.varexpl.both       = params.analysis.spatial.rh.varexpl(alignMe.pRFs.selected.rh);
    end
    
    else
        data.predRuns.both = predRunBothHemi;
        data.rf.both = predAllRuns{1}.pred.rf;
        data.params.x0.both =  predAllRuns{1}.pred.params.analysis.spatial.x0;
        data.params.y0.both =  predAllRuns{1}.pred.params.analysis.spatial.x0;
        data.params.effectiveSize.both =  predAllRuns{1}.pred.params.analysis.spatial.sigmaMajor ./ ...
                                        sqrt(predAllRuns{1}.pred.params.analysis.spatial.exponent);
        data.params.exponent.both = cat(2,predAllRuns{1}.pred.params.analysis.spatial.lh.exponent, ...
                                        predAllRuns{1}.pred.params.analysis.spatial.rh.exponent);
        data.params.exponent_temporal.both = predAllRuns{1}.pred.params.analysis.temporal.exponent;
        xx = load(fullfile(saveBetaFolder,'betaR2_XValFit_OLS_gridFit', ...
            sprintf('betaValsFull_%s_%s_%s_%s_run12_offsetFlag1_OLS_bslCorrected_expBestR2.mat', ...
                    rois{roi},roiType,spatialModel, temporalModel)));
        data.params.varexpl.both =  xx.optimal.data.params.varexpl.both;
        clear xx
    end
    
    % Check if nr of voxels in data match predictions and clear
    % some memory
    assert(size(data.meanDetrendedData,2)==size(data.predRuns.both,2))
    assert(size(data.meanDetrendedData,2)==size(data.rf.both,3))
    clear predRunBothHemi;
    
    % Check flags and if data is actually present
    fn = fieldnames(alignMe.pRFs.selected);
    if length(fn)==1 && isempty(alignMe.pRFs.selected.(fn{1}))
        break
    elseif length(fn)==2 && isempty(alignMe.pRFs.selected.(fn{1})) && isempty(alignMe.pRFs.selected.(fn{2}))
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
        Y = data.meanDetrendedData;
        X = data.predictionCatRuns;
        [nFrames, numVoxels] = size(data.meanDetrendedData);
        
        resultsFullFit = fitModelToVoxData(X,Y,...
            'regressionType',params.analysis.regressionType, ...
            'sumSTFlag', sumSTFlag, ...
            'weightChannels', weightChannels, ...
            'addOffsetFlag', addOffsetFlag);
        
        if doCrossValidation
            % If we want cross-validation, we have to add an offset to
            % so we can remove the mean before computing R2
            addOffsetFlagCurrent = addOffsetFlag;
            addOffsetFlag = true;
            clear B_crossval
            
            % Split data in even and odd concatenated runs
            Y_crossval = {data.meanDetrendedDataOdd, data.meanDetrendedDataEven};
            
            % If we use fracridge, use best alpha from full fit.
            if strcmp(params.analysis.regressionType,'fracridge')
                bestAlpha = resultsFullFit.lm{1}.bestAlpha;
            else, bestAlpha = []; end
            
            % Loop over split halves
            for cx = [1,2]
                resultsCValInitFit{cx} = fitModelToVoxData(...
                    X,Y_crossval{cx},...
                    'regressionType',params.analysis.regressionType, ...
                    'sumSTFlag', sumSTFlag, ...
                    'weightChannels', weightChannels, ...
                    'addOffsetFlag', addOffsetFlag, ...
                    'alpha',bestAlpha);
            end
          
            % Cross-validate modelfits, i.e. use scaled modelfit from split A
            % calculate variance explain with data in split B and vice versa)
            R2_crossval = NaN(2, numVoxels);
            for cx = [1,2]
                % We refit beta weights to data where baseline is removed
                resultsCValFinalFit{cx} = fitModelToVoxData(...
                    X,resultsCValInitFit{cx}.Y_noOffset,...
                    'regressionType',params.analysis.regressionType, ...
                    'sumSTFlag', sumSTFlag, ...
                    'weightChannels', weightChannels, ...
                    'addOffsetFlag', addOffsetFlag, ...
                    'alpha',bestAlpha);
                B_crossval(:,:,cx) = resultsCValFinalFit{cx}.B;
        
                % Get opposite dataset
                fitIdx = setdiff([1,2],cx);
                
                % First, get prediction and data final fit, no offset subtracted
                Xcx = resultsCValFinalFit{cx}.sumChannelPrediction;
                Ycx = resultsCValInitFit{fitIdx}.Y_noOffset;
                
                % Now compute Coefficient of Determination (R2) using split
                % half data B, and modelfit A, and vice versa.
                for n = 1:size(Ycx,2)
                    R2_crossval(cx,n) = computeCoD(Ycx(:,n),Xcx(:,n));
                end
                
                % Also store predictions of final fit with, and without offset
                sumChannelPredictionCVNoOffset{cx} = resultsCValFinalFit{cx}.sumChannelPredictionNoOffset;
                sumChannelPredictionCV{cx}         = resultsCValFinalFit{cx}.sumChannelPrediction;
            end
            
            B_crossval_mean  = mean(B_crossval,3,'omitnan');
            R2_crossval_mean = mean(R2_crossval,1,'omitnan');
            sumChannelPredictionCrossvalMean = mean(cat(3,sumChannelPredictionCV{1},sumChannelPredictionCV{2}),3,'omitnan');
            Y_crossval_mean = mean(cat(3,resultsCValInitFit{1}.Y_noOffset,resultsCValInitFit{2}.Y_noOffset),3,'omitnan');
            
            sumChannelPredictionNoOffsetCrossvalMean = mean(cat(3,sumChannelPredictionCVNoOffset{1},sumChannelPredictionCVNoOffset{2}),3,'omitnan');
            Y_crossval_mean_noOffset = mean(cat(3,resultsCValFinalFit{1}.Y_noOffset,resultsCValFinalFit{2}.Y_noOffset),3,'omitnan');
            
            % Restore addOffsetFlag
            addOffsetFlag = addOffsetFlagCurrent; 
            clear addOffsetFlagCurrent
        end
        
        %% Plot run
        if verbose
            % Select top 30 voxels based on R2
            if doCrossValidation, 
                dataToSort = R2_crossval_mean;
            else, dataToSort = resultsFullFit.R2; end
            dataToSort(isnan(dataToSort))=0; 
            [~,idx] = sort(dataToSort, 'descend');
            if numVoxels < 30; numVoxToPlot = idx;
            else, numVoxToPlot = idx(1:30); end
            
            if saveFigs
                folderName = sprintf('%s_%s_%s_%s_AverageRun_addOffsetFlag%d_cval%d', ...
                    rois{roi},roiType,spatialModel, temporalModel, addOffsetFlag, doCrossValidation);
                if ~isempty(useFixedExponent)
                    folderName = sprintf('%s_exp%1.1f',folderName, useFixedExponent);
                elseif useSearchFit
                    folderName = sprintf('%s_%s_searchFit',folderName, searchAlgorithm);
                end
                subFolderDir = fullfile(saveFigDir, 'figs', folderName);
                if ~exist(subFolderDir,'dir'), mkdir(subFolderDir); end
            end
            
            fH = figure(101); clf; set(gcf, 'position', [55,440,1994,522]);
            stimToPlot = normMax(stimRunUnique(:))';
            for ii = numVoxToPlot
                if doCrossValidation
                    dataToPlot      = Y_crossval_mean(:,ii);
                    modelPredToPlot = sumChannelPredictionCrossvalMean(:,ii);
                    betaToPrint     = B_crossval_mean(ii,:);
                    R2ToPrint       = R2_crossval_mean(ii);
                else
                    dataToPlot      = data.meanDetrendedData(:,ii);
                    modelPredToPlot = resultsFullFit.sumChannelPrediction(:,ii);
                    betaToPrint     = resultsFullFit.B(ii,:);
                    R2ToPrint       = resultsFullFit.R2(ii);
                end
                if useSearchFit || ~isempty(useFixedExponent) && ~isnan(useFixedExponent)
                    exponentToPrint = data.params.exponent_temporal.both(ii);
                else
                    exponentToPrint = data.params.exponent.both(ii);
                end
                ttl = sprintf('Voxel %d %s %s: [x,y,sz,exp]=[%1.1f,%1.1f,%1.2f,%1.2f]. %s %s: beta %s - R2 %2.1f%%', ...
                    ii, rois{roi}, roiType, data.params.x0.both(ii), ...
                    data.params.y0.both(ii), data.params.effectiveSize.both(ii), ...
                    exponentToPrint, spatialModel, temporalModel, mat2str(betaToPrint,2), ...
                    R2ToPrint*100);
                
                fH = plotFullRunTSwithModelFit(fH, dataToPlot, modelPredToPlot, ...
                    stimToPlot,ttl);
                
                if saveFigs
                    currFName = sprintf('voxel%d_fullAveData_vs_model',ii);
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
        
        % Find preprocessed datatype
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
        
        % Update nFrames in preprocessed datatype if needed
        if dataTYPES(dt_idx).scanParams.nFrames ~= nFrames
            dataTYPES(dt_idx).scanParams.nFrames = nFrames;
            saveSession;
        end
        
        % Update par file in preprocessed datatype if needed
        tmp = dataTYPES(dt_idx).scanParams.parfile(1:end-4);
        [~,dtmp] = fileparts(tmp);
        if isempty(regexp(dtmp, '_corrected2','once'))
            dataTYPES(dt_idx).scanParams.parfile = fullfile([tmp '_corrected2.par']);
            saveSession;
        end
        
        % Get concatenated parfile trials struct
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
        [T1odd,allDataTrialsOdd] = ...
            simseq_getStimBlocksFromDataRuns(data.meanDetrendedDataOdd,trials,er_data(1).params); %#ok<ASGLU>
        [T1even,allDataTrialsEven] = ...
            simseq_getStimBlocksFromDataRuns(data.meanDetrendedDataEven,trials,er_data(1).params); %#ok<ASGLU>
        
        % Update data table (since we removed some voxels because they
        % were eliminated by modelfit
        T.DataTS          = TD1.TS;
        T.DataTSerror     = TD1.TSerror;
        T.DataTSOdd       = T1odd.TS;
        T.DataTSOdderror  = T1odd.TSerror;
        T.DataTSEven      = T1even.TS;
        T.DataTSEvenerror = T1even.TSerror;
        
        if doCrossValidation
            % DATA: WITH Cross-Validation. Y_crossval_mean => allDataTrialsCrossvalMean meanDataTrialsCrossvalMean semDataTrialsCrossvalMean
            [TD_CV,allDataTrialsCrossvalMean,~,~,~,catTrialsDataCV] = ...
                simseq_getStimBlocksFromDataRuns(Y_crossval_mean,trials,er_data(1).params); %#ok<ASGLU>
            T.DataTSCV      = TD_CV.TS;
            T.DataTSCVError = TD_CV.TSerror;
            
            % MODEL: WITH Cross-Validation but NO removed offset (since we didn't add it to the model)
            % sumChannelPredictionCrossvalMean => allModelTrialsCrossvalMean meanModelTrialsCrossvalMean semModelTrialsCrossvalMean
            [TM_CV,allModelTrialsCrossvalMean,~,~] = ...
                simseq_getStimBlocksFromDataRuns(sumChannelPredictionCrossvalMean,trials,er_data(1).params); %#ok<ASGLU>
            T.ModelTSCVMn    = TM_CV.TS;
            T.ModelTSCVError = TM_CV.TSerror;
            
            % MODEL: WITH Cross-Validation WITH removed offset (as we added it to the model)
            % sumChannelPredictionNoOffsetCrossvalMean => allModelTrialsCrossvalMeanNoOffset meanModelTrialsCrossvalMeanNoOffset semModelTrialsCrossvalMeanNoOffset
            [TM_CV_noOffset,allModelTrialsCrossvalMeanNoOffset,~,~] = ...
                simseq_getStimBlocksFromDataRuns(sumChannelPredictionNoOffsetCrossvalMean,trials,er_data(1).params); %#ok<ASGLU>
            T.ModelTSCVMn_noOffset    = TM_CV_noOffset.TS;
            T.ModelTSCVError_noOffset = TM_CV_noOffset.TSerror;
        end  
        
        % MODEL: NO cross-Validation and NO removed offset (since we didn't add it to the model)
        % sumChannelPrediction => allModelTrials meanModelTrials semModelTrials
        [TM,allModelTrials,~,~,~, ~] = ...
            simseq_getStimBlocksFromDataRuns(resultsFullFit.sumChannelPrediction,trials,er_data(1).params); %#ok<ASGLU>
        T.ModelTS      = TM.TS;
        T.ModelTSError = TM.TSerror;

        if addOffsetFlag
            % MODEL: NO cross-Validation WITH removed offset (as we added it to the model)
            % sumChannelPredictionNoOffset = allModelNoOffsetTrials meanModelNoOffsetTrials semModelNoOffsetTrials
            [TM_noOffset,allModelNoOffsetTrials,~,~] = ...
                simseq_getStimBlocksFromDataRuns(resultsFullFit.sumChannelPredictionNoOffset,trials,er_data(1).params);  %#ok<ASGLU>
            T.ModelTS_noOffset      = TM_noOffset.TS;
            T.ModelTSError_noOffset = TM_noOffset.TSerror;
        end  
        
        clear TD_CV TM_CV TM_CV_noOffset TM TM_noOffset;
        
        % Recompute Split half reliability if we removed data voxels..
        % .. across runs
        data.splitHalfRel.both.run = simseq_getSplitHalfDataReliability( ...
            permute(data.detrendedData.both,[1,3,2]), 'runFlag', true, ...
            'splitRuns', pths.runOrder);
        
        % .. across trials
        data.splitHalfRel.both.trial = simseq_getSplitHalfDataReliability( ...
            catTrialsData, 'runFlag',false,'splitRuns',[]);
        
        if doCrossValidation 
            % .. across runs
            data.splitHalfRel.both.runCV = simseq_getSplitHalfDataReliability( ...
                permute(data.dataRuns.detrendedStimRuns,[1,3,2]), ...
                'runFlag',true, 'splitRuns', []);
            % .. across trials
            data.splitHalfRel.both.trialCV = simseq_getSplitHalfDataReliability( ...
                catTrialsDataCV, 'runFlag',false,'splitRuns',[]);
        end
        
        %% Get R2 for model prediction at trial level, per condition
        numConditions = length(T.Condition);
        R2_full_trial = NaN(numConditions,numVoxels);
        for cond = 1:numConditions
            currData  = cell2mat(T.DataTS(T.Condition==cond));
            currModel = cell2mat(T.ModelTS(T.Condition==cond));
            if doCrossValidation
                currDataCV  = cell2mat(T.DataTSCV(T.Condition==cond));
                currModelCV = cell2mat(T.ModelTSCVMn(T.Condition==cond));
            end
            for vv = 1:size(currData,2)
                R2_full_trial(cond,vv) = computeCoD(squeeze(currData(:,vv)), squeeze(currModel(:,vv)));
                if doCrossValidation
                    R2_crossval_trial(cond,vv) = computeCoD(squeeze(currDataCV(:,vv)), squeeze(currModelCV(:,vv)));
                end
            end
        end
        clear currData currModel currDataCV currModelCV
        
        % Take mean and SEM across voxels
        meanR2_full_trial = mean(R2_full_trial,2,'omitnan'); %#ok<NASGU>
        seR2_full_trial   = std(R2_full_trial,[],2,'omitnan')./numVoxels; %#ok<NASGU>
        if doCrossValidation
            meanR2_crossval_trial = mean(R2_crossval_trial,2,'omitnan'); %#ok<NASGU>
            seR2_crossval_trial   = std(R2_crossval_trial,[],2,'omitnan')./numVoxels; %#ok<NASGU>
        end
        
    end % if fieldnames exist
    
    %% save betas
    if saveBetas
        
        saveList = {'T','resultsFullFit', ...
            'allDataTrials','allModelTrials','data','nrTrials','er_data',...
            'R2_full_trial', 'meanR2_full_trial','seR2_full_trial',...
            'stimRunUnique','alignMe','params'};
        if ~useSearchFit
            saveList = {saveList{:},'T_orig'};
        end
        if addOffsetFlag
            saveList = {saveList{:},'allModelNoOffsetTrials'};
        end
        
        if doCrossValidation
            saveList = {saveList{:},'B_crossval','B_crossval_mean', ...
                'R2_crossval','R2_crossval_mean', ...
                'R2_crossval_trial',...
                'meanR2_crossval_trial','seR2_crossval_trial',...
                'allModelTrialsCrossvalMean', ...
                'allDataTrialsCrossvalMean',...
                'resultsCValInitFit', ...
                'sumChannelPredictionCrossvalMean', ...
                'Y_crossval','Y_crossval_mean'};
            
            if addOffsetFlag
                saveList = {saveList{:},'allModelTrialsCrossvalMeanNoOffset',...
                    'sumChannelPredictionNoOffsetCrossvalMean',...
                    'resultsCValFinalFit', ...
                    'Y_crossval_mean_noOffset'};
            end
        end
        
        betaFolder0 = 'betaR2_XValFit_OLS';
        if useSearchFit
            betaFolder = [betaFolder0 '_finalFit_' searchAlgorithm '_restricted'];
            if ~exist(fullfile(saveBetaFolder,betaFolder),'dir')
                mkdir(fullfile(saveBetaFolder,betaFolder));
            end
            
            loadFile = fullfile(saveBetaFolder,betaFolder,...
                sprintf('%s_%s_%s_%s_%s_run%s_addOffsetFlag%d_%s_bslCorrected_finalFit_%s.mat', ...
                pths.subjID,rois{roi},roiType,spatialModel, temporalModel, postFix, addOffsetFlag, params.analysis.regressionType, searchAlgorithm));
 
        elseif ~isempty(useFixedExponent) && ~isnan(useFixedExponent)
            betaFolder = [betaFolder0 '_gridFit'];
            if ~exist(fullfile(saveBetaFolder,betaFolder),'dir')
                mkdir(fullfile(saveBetaFolder,betaFolder));
            end
            loadFile = fullfile(saveBetaFolder,betaFolder,...
                sprintf('%s_%s_%s_%s_%s_run%s_addOffsetFlag%d_%s_bslCorrected_exp%1.2f.mat', ...
                pths.subjID,rois{roi},roiType,spatialModel, temporalModel, postFix, addOffsetFlag, params.analysis.regressionType, useFixedExponent));
        
        else
            loadFile = fullfile(saveBetaFolder,betaFolder0,...
                sprintf('%s_%s_%s_%s_%s_run%s_addOffsetFlag%d_%s_bslCorrected.mat', ...
                pths.subjID,rois{roi},roiType,spatialModel, temporalModel, postFix, addOffsetFlag, params.analysis.regressionType));
        end
        save(loadFile, saveList{:},'-v7.3');
    end % saveBetas

end % ROIS

return





