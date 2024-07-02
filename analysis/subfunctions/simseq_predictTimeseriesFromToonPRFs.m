function [] = simseq_predictTimeseriesFromToonPRFs(subjnr, projectDir, stimFile, varargin)
% Function to generate temporal model predictions using pRF of single vertex
% from toonotopy experiment.
%
%       simseq_predictTimeseriesFromToonPRFs(subjnr, projectDir, varargin)
%
% INPUTS
% subjnr                    : (int) subject to analyze
% projectDir                : (str) where does project live
% stimFile                  : (str) path to file with stimulus sequence to
%                               make predicted time series. Needs to
%                               contain 'images' and 'sequence', where
%                               images is an array [X, Y, time], sequence
%                               is a vector [1xtime] with condition nrs.
% [dt]                      : (str) vistasoft dataTYPE
% [spatialModel]            : (str) spatial prf model. Choose from
%                              'cssFit' or 'onegaussianFit' (default)
% [temporalModel]           : (str) temporal prf model. Choose from
%                              '1ch-glm' (linear), (default)
%                              '1ch-dcts' (Divisive normalization),
%                              '2ch-exp-sig' (sustained adaptation channel
%                                     + transient sigmoidal on/off channel)
%                              '2chan-css-sig' (sustained CSS-like exponent
%                                     channel + transient sigmoidal on/off
%                                     channel)
%                              '3ch-stLN' (3-channel spatiotemporal
%                              linear-filterbank follwed by static
%                               compressive nonlinearity)
% [hemi]                    : (str) hemisphere to analyze, 'lh'(default),
%                              'rh', or 'both'
% [roiIdx]                  : (int/str) use 'all' ROIs or a specific index.
%                               see getSubjectPaths for defined ROIs.
%                               (default = 'all')
% [roiType]                 : (str) within ROI, one can pick a subROI.
%                               Choose from stimcorner'(default),'square1',
%                               'square2','square3','square4'
% [veThresh]                : (int) variance explained threshold, (default
%                               =0.1, so more than 10%)
% [stimRun]                 : (int) what runs to analyze, if average, use 1
% [upsampleTimeStim]        : (str) upsample stimulus time to have 1ms
%                               resolution (default = 'none'); Choose from:
%                               'upsampleFrom1Hz','upsampleFrom60Hz','none'
% [verbose]                 : (bool) plot figures or not (default=true)
% [saveFigs]                : (bool) save figures or not (default=true)
% [savePredictionsFlag]     : (bool) save predictions or not (default=true)
% [trimRFFlag]              : (bool) trim 2D pRF estimate at 5 SD or not
%                               (default=true)
% [useArtificialPRFs]       : (bool) use artifical pRFs to test
%                               computation, instead of estimates from data
%                               (default=false)
% [useMedianROIexponent]    : (bool) use median exponent across roi voxels,
%                               instead of individual roi voxels in
%                               compressive nonlinearity.
% [useFixedExponent]        : (int) single value to make model prediction
%                               (ST model) with specific exponent (like a
%                               grid fit)
% [combineNeuralChan]       : (int) vector with length equal to number of
%                               neural channels. If you want to combine
%                               neural channels before convolving with HRF,
%                               then give the location of the channel the
%                               same unique number. For example, [1 2 2]
%                               combines the last two channels into one,
%                               by summing and then rescaling max height to
%                               1.
% [useSTRetParams]          : (bool) use parameters from spatiotemporal
%                               retinotopy experiment
%
% % Examples:
% subjnr = 1;
% projectDir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/';
%
% Example 1: CSS + 2-chan temporal model
% simseq_predictTimeseriesFromToonPRFs(subjnr, projectDir);
%
% % Example 2: CSS + 2-chan, no verbose
% simseq_predictTimeseriesFromToonPRFs(subjnr, projectDir, 'verbose',false);
%
% % Example 3: CSS + 1-chan GLM
% simseq_predictTimeseriesFromToonPRFs(subjnr, projectDir,'temporalModel','1ch-glm');
%
% % Example 4: Standard Gaussian model + 2-chan temporal model
% simseq_predictTimeseriesFromToonPRFs(subjnr, projectDir, 'spatialModel','onegaussianFit');
%
% % Example 5: Standard Gaussian model + 1-chan GLM
% simseq_predictTimeseriesFromToonPRFs(subjnr, projectDir,'spatialModel','onegaussianFit','temporalModel','1ch-glm');
%
% % Example 6: 3-chan spatiotemporal filterbank + static nonlinearity
% simseq_predictTimeseriesFromToonPRFs(subjnr, projectDir,'spatialModel','onegaussianFit','temporalModel','3ch-stLN');

%
% Written by EK, 2021, Stanford U

%% Parse inputs
p = inputParser;
p.addRequired('subjnr', @isnumeric);
p.addRequired('projectDir',@ischar);
p.addRequired('stimFile',@ischar);
p.addParameter('dt','Averages', @ischar);
p.addParameter('spatialModel','onegaussianFit', ...
    @(x) any(validatestring(x,{'cssFit','onegaussianFit', 'differenceOfGaussiansFit'})));
p.addParameter('temporalModel','1ch-glm', ...
    @(x) any(validatestring(x,{'1ch-glm','1ch-dcts','2ch-exp-sig', '2ch-css-sig','3ch-stLN'})));
p.addParameter('hemi','lh', @(x) any(validatestring(x,{'lh','rh', 'both'})));
p.addParameter('roiIdx','all', @(x) (ischar(x) || isnumeric(x) || isscalar(x) || isvector(x)));
p.addParameter('roiType','stimcorner4_area4sq_eccen5', @ischar); %@(x) any(validatestring(x, ...
%     {'stimcorner4_area4sq_eccen5', 'stimcorner4_area4sq_eccen5_stRet_CSTopt_DNST_matchingVoxels'})));
p.addParameter('veThresh',0.1, @isnumeric);
p.addParameter('upsampleStimType', 'none', @(x) any(validatestring(x,{'none','upsampleFrom1Hz', 'upsampleFrom60Hz'})));
p.addParameter('stimRun', [], @isnumeric);
p.addParameter('saveFolder', [], @ischar);
p.addParameter('verbose', true, @islogical);
p.addParameter('saveFigs', true, @islogical);
p.addParameter('savePredictionsFlag',true, @islogical);
p.addParameter('trimRFFlag', true, @islogical);
p.addParameter('useArtificialPRFs', false, @islogical);
p.addParameter('useMedianROIexponent', false, @islogical);
p.addParameter('useAvgPRFs',false, @islogical)
p.addParameter('useFixedExponent',[],@isnumeric);
p.addParameter('combineNeuralChan',[],@isnumeric);
p.addParameter('subsampleRateMs',1,@isnumeric);
p.addParameter('useSTRetParams',false,@islogical);
p.parse(subjnr,projectDir,stimFile,varargin{:});

% Rename variables
fnToRename = fieldnames(p.Results);
for ff = 1:length(fnToRename)
    eval([sprintf('%s = p.Results.(fnToRename{ff});', fnToRename{ff})]);
end
clear fnToRename ff


%% Print model type
fprintf('[%s]: Running spatial model %s + temporal model %s for %s hemi subject nr %d\n', ...
    mfilename, spatialModel, temporalModel, hemi, subjnr)

%% Set up paths and variables

% cd to subject's session folder.
cd(fullfile(projectDir,'experiments/simseq/'))

% General params
sesNr = getSessionNrMainExp(subjnr);
pths  = getSubjectPaths(projectDir,subjnr,sesNr);

% cd to toon data
cd(fullfile(pths.dataDirToon, pths.toon))

% Create folders if needed
if savePredictionsFlag
    if isempty(saveFolder)
        [~,saveSubFolder] = fileparts(stimfile);
        saveFolder = fullfile(pths.simseqResultsDir,sprintf('savedPredictions_%s', saveSubFolder));
    end
    if ~exist(saveFolder,'dir'); mkdir(saveFolder); end
end

% Get ROI depending on what hemisphere to process; left, right, or both
if strcmp(hemi, 'both')
    rois     = pths.definedROIsBOTH;
elseif strcmp(hemi, 'lh')
    rois     = pths.definedROIsLH;
elseif strcmp(hemi, 'rh')
    rois     = pths.definedROIsRH;
end

% If we just want to do a test run or simulation we use artifical pRFs
if useArtificialPRFs
    rois = {'artificial'};
end

% Get selected roi indices
if ~strcmp(roiIdx,'all') && isnumeric(roiIdx)
    rois = rois(roiIdx);
elseif ~strcmp(roiIdx,'all') && ischar(roiIdx)
    rois = find(ismember(rois,{roiIdx}));
end

%% Get ROI data and select pRF params within ROI
for r = 1:length(rois) %find(ismember(rois,{'IPS0', 'IPS1'})) %% %
    
    % define params struct
    params = struct();
    params.verbose = verbose;
    params.saveFigs = saveFigs;
    params.useGPU = false;
    
    % Load pRFs
    if useArtificialPRFs
        prfParams = simseq_simulatePRFParams(hemi, spatialModel, veThresh);
    elseif useAvgPRFs
        if strcmp(hemi,'both')
            hm = {'lh','rh'};
        else
            hm = {hemi};
        end
        
        for h = 1:length(hm)
            avgPRFfolder = fullfile(pths.simseqResultsDir,'preprocData');
            subf = strsplit(saveFolder,'/');
            subf = subf{end}(1:end-4);
            load(fullfile(avgPRFfolder, [subf '2_v' mat2str(pths.expversionNr)], 'avgPRF',...
                sprintf('avgPRF_%s_%s_%s.mat', hm{h},rois{r},roiType)), 'avgPRF')
            
            prfParams.(hm{h}) = avgPRF.(hm{h});
        end
    elseif useSTRetParams && (strcmp(temporalModel, '1ch-dcts') || strcmp(temporalModel, '3ch-stLN'))
        roifile  = [rois{r} '_toon_simseq' roiType];
        prfParams = simseq_loadstRetPRFParams(pths, hemi, roifile, spatialModel, temporalModel, veThresh, [], 'simseqROIs');
    else
        roifile  = [rois{r} '_toon_simseq' roiType];
        [prfParams, hvol] = simseq_loadPRFParams(pths, dt, hemi, roifile, 'cssFit', veThresh, [],'simseqROIs');
    end
    
    % Store params
    params.analysis.spatial = prfParams;
    if isfield(prfParams.rh,'hvolROI') && ~isempty( prfParams.rh.hvolROI.coords)
        params.analysis.temporal.rh = prfParams.rh.temporal;
    else
        params.analysis.temporal.rh = [];
    end
    if isfield(prfParams.lh,'hvolROI') && ~isempty( prfParams.lh.hvolROI.coords)
        params.analysis.temporal.lh = prfParams.lh.temporal;
    else
        params.analysis.temporal.lh = [];
    end
    clear prfParams;
    params.analysis.subsampleRateMs = subsampleRateMs;
    
    fNames = fieldnames(params.analysis.spatial);
    params.spatial.trimRFFlag = trimRFFlag;
    if isfield(params.analysis.spatial.rh, 'x0') || isfield(params.analysis.spatial.rh, 'y0') && ...
            (length(params.analysis.spatial.rh.x0) + length(params.analysis.spatial.rh.y0))>1
        % For 3channel-spatiotemporal linear-nonlinear model, we have several
        % options for the static powerlaw nonlinearity:
        % * exponent from spatial CSS model, preserved for each voxel
        % * exponent from spatial CSS model, median across voxels in ROI if "useMedianROIexponent" is true.isfield(params.analysis.spatial.rh, 'x0') || isfield(params.analysis.spatial.rh, 'y0') && ...
        if  strcmp(temporalModel,'3ch-stLN')
            if ~useSTRetParams
                % Default: Use exponent from spatial CSS model, preserved for each voxel
                for fn = 1:length(fNames)
                    params.analysis.temporal.param.exponent.(fNames{fn}) = params.analysis.spatial.(fNames{fn}).exponent;
                end
            end
            % Combine neural channels if requested, we usually combine the on-
            % and off- transient channels into one transient channel as their
            % final predicted BOLD responses are very similar and we want to
            % avoid collinearity when fitting the predictions.
            if isempty(combineNeuralChan)
                params.analysis.combineNeuralChan     = [1:params.analysis.temporal.num_channels];
            else
                params.analysis.combineNeuralChan     = combineNeuralChan;
            end
        end
        
        % We always load the CSS model fit params, so we use sigma major and
        % minor with the effective size (sigma/sqrt(n)) for onegaussianFits and
        % reset the spatial exponent to 1.
        if strcmp(spatialModel,'onegaussianFit') && ~strcmp(temporalModel,'1ch-dcts') && ~useSTRetParams
            for fn = 1:length(fNames)
                params.analysis.spatial.(fNames{fn}).sigmaMajor = params.analysis.spatial.(fNames{fn}).effectiveSize;
                params.analysis.spatial.(fNames{fn}).sigmaMinor = params.analysis.spatial.(fNames{fn}).effectiveSize;
                params.analysis.spatial.(fNames{fn}).exponent = ones(size(params.analysis.spatial.(fNames{fn}).exponent));
            end
        end
        
        % Store params
        params.analysis.temporalModel        = temporalModel;
        params.analysis.spatialModel         = spatialModel;
        params.analysis.useMedianROIexponent = useMedianROIexponent;
        params.analysis.useFixedExponent     = useFixedExponent;
        params.analysis.useAvgPRFs           = useAvgPRFs;
        params.analysis.useSTRetParams       = useSTRetParams;
        
        %% Create stim from simseq experiment and add to params  ----------------------
        if ~exist(fullfile(stimFile),'file')
            error('[%s]: Can''t file stimulus file',mfilename);
        end
        
        %% Upsample stimulus to millisecond resolution
        
        % Add stim params
        params.stim.framePeriod     = 1; %s (sample rate of Simseq fMRI data)
        params.stim.imFile          = stimFile;
        params.stim.paramsFile      = stimFile;
        params.stim.prescanDuration = 0; % prescan duration is already clipped during preprocessing
        
        % Add other inhertied analysis params
        for fn = 1:length(fNames)
            params.analysis.numberStimulusGridPoints = params.analysis.spatial.(fNames{fn}).nSamples; % from hvol.rm.retinotopyParams.analysis.numberStimulusGridPoints; % 50 corresponds to 101 x 101
            params.analysis.sampleRate = params.analysis.spatial.(fNames{fn}).sampleRate; % hvol.rm.retinotopyParams.analysis.sampleRate;
            params.analysis.fieldSize  = params.analysis.spatial.(fNames{fn}).fieldSize; % hvol.rm.retinotopyParams.analysis.fieldSize; % deg visual angle
        end
        params.analysis.doBlankBaseline = false;
        
        switch upsampleStimType
            case 'upsampleFrom1Hz'
                params.analysis.pRFmodel = {'st'};
            case 'upsampleFrom60Hz'
                params.analysis.pRFmodel = {'st-upsampling60hz'};
            case 'none'
                params.analysis.pRFmodel = {'st-nostimtimeupsampling'};
        end
        
        % Make 1ms stimulus
        params = makeStiminMS(params);
        params.stim.images_unconvolved = params.stim.images;
        for fn = 1:length(fNames)
            params.analysis.spatial.(fNames{fn}).X = params.analysis.X;
            params.analysis.spatial.(fNames{fn}).Y = -1*params.analysis.Y;
        end
        %     params.analysis.X = []; params.analysis.Y = [];
        
        % File to save data to
        if ~isempty(params.analysis.useFixedExponent)
            params.analysis.predFile = fullfile(saveFolder, ...
                sprintf('modelPredictions_%s_%s_%s_%s_run%d_exp%1.2f.mat',rois{r},roiType,spatialModel, temporalModel, stimRun,params.analysis.useFixedExponent));
        elseif useAvgPRFs && length(hm)==1
            params.analysis.predFile = fullfile(saveFolder, ...
                sprintf('modelPredictions_%s_%s_%s_%s_%s_run%d.mat',hm{1},rois{r},roiType,spatialModel, temporalModel, stimRun));
        else
            params.analysis.predFile = fullfile(saveFolder, ...
                sprintf('modelPredictions_%s_%s_%s_%s_run%d.mat',rois{r},roiType,spatialModel, temporalModel, stimRun));
        end
        %% Start making predictions
        if (isfield(params.analysis.spatial.rh,'x0') || isfield(params.analysis.spatial.lh,'x0')) 
            if  (length(params.analysis.spatial.rh.x0)>0 || length(params.analysis.spatial.lh.x0)>0)
                pred = stPredSimSeqWrapper(params);
                
                % Store predictions and params
                pred.params0 = params;
                if savePredictionsFlag
                    fprintf('[%s]: Saving predictions for %s %s run %d %s\n',mfilename, spatialModel, temporalModel,stimRun,rois{r})
                    if exist(fullfile(params.analysis.predFile),'file')
                        str = sprintf('delete %s',fullfile(params.analysis.predFile));
                        eval(str);
                    end
                    saveFile = fullfile(params.analysis.predFile);
                    save(saveFile, 'pred','-v7.3');
                end
                close all;
            end
        end
    end
end
