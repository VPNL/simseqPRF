%% s_runAllModelFits.m
% Main script to fit predict time series to measured time series.
%
% Some notes about this script:
% * previous version of this function is:
% simseq_fitTemporalModelPredictiosToDataWithPreproc('stimFileName',stimFileName, ...)
% * We remove 10 TRs at the start of each predicted run's time series and
% the stimulus sequence, (not the preprocessed data, see s_runDataPreproc.m)
% because the data preprocessing removed 6 TRs (1 TR = 1 sec) countdown
% + 2 TRs of the baseline. Here we remove 10TRs, so two extra TRs, because
% we observed that there is still gradient instability at the start of each
% time series run. Keeping this introduces artefacts in baseline correction
% and results in poorer modelfits.
%
% Main experiment subject numbers = [1,2,3,7,8,9,10,11,12,13]
%
%
% Written by Eline R Kupers, Stanford U 2021

%% Define paths and parameters
projectDir = fullfile(simseqRootPath);
useSTRetParams  = true; % If true, we use solved parameters from spatiotemporal retinotopy paper for 3ch-stLN model
roiType0        = 'stimcorner4_area4sq_eccen5';
if useSTRetParams
    subjnrs = [1:3,7:10];
    subDir = 'savedPredictions_fixedPRF_combT_matchingVoxels';
    spatialModels  = {'onegaussianFit','onegaussianFit'};
    temporalModels = {'3ch-stLN','1ch-dcts'};  % CSTopt & DNST
else
    subjnrs = [1:3,7:13];
    subDir = 'savedPredictions_fixedPRF_combT';
    spatialModels  = {'onegaussianFit','cssFit','onegaussianFit'};
    temporalModels = {'1ch-glm','1ch-glm','3ch-stLN'};
end

for subjnr = subjnrs
    
    pths  = getSubjectPaths(projectDir,subjnr);
    
    % Define inputs such as folders
    stimFileName   = 'stim_simseq_run';
    subFolder      = 'varyStimDur_x_Area4_2x2design_fullrun';
    saveBetaFolder = fullfile(pths.simseqResultsDir,'savedPredictions_fixedPRF_combT');
    loadDataFolder = fullfile(pths.simseqResultsDir,'preprocData');
    saveFigDir     = fullfile(pths.figureDir, 'cvfits');
    removedTRsAtStart = 2; % Preprocessing EPIs removed 8 TRs (= 6 TRs countdown + 2 TRs baseline),
    % Here, we remove two extra TRs, since the countdown is not part of model prediction
    
    for ii = 1:length(temporalModels)
        
        if strcmp(temporalModels{ii},'1ch-glm') && ~useSTRetParams
            % Standard Gaussian model
            useFixedExponent = [];
            exponents        = NaN; % NaN means no fixed exponent in modelfit: Use NaN for linear and CSS model
            combChan         = [];
            loadDataFolder   = loadDataFolder0;
            saveBetaFolder   = saveBetaFolder0;
            
        elseif strcmp(temporalModels{ii},'1ch-dcts') && useSTRetParams
            % Standard Gaussian model + 1-chan Divisive Normalization spatiotemporal model
            useFixedExponent = [];
            combChan         = [];
            exponents        = NaN; % NaN means no fixed exponent in modelfit: Use NaN for linear and CSS model
            loadDataFolder   = [loadDataFolder0 '_stRet_matchingVoxels_' temporalModels{ii}];
            roiType          = sprintf('%s_%s',roiType0,'stRet_CSTopt_DNST_matchingVoxels');
            saveBetaFolder   = [saveBetaFolder0 '_matchingVoxels'];
            
        elseif strcmp(temporalModels{ii},'3ch-stLN') && ~useSTRetParams
            % Standard Gaussian model + 3-chan linear-nonlinear spatiotemporal model
            combChan       = [1 2 2];
            exponents      = [0.1:0.05:1]; % Grid fit: systematically sweep 0-1
            loadDataFolder = loadDataFolder0;
            saveBetaFolder   = saveBetaFolder0;
            
        elseif strcmp(temporalModels{ii},'3ch-stLN') && useSTRetParams
            % Standard Gaussian model + 3-chan linear-nonlinear
            % spatiotemporal model with parameters from ST ret exp
            combChan       = [1 2 2];
            exponents      = NaN; % NaN means no fixed exponent in modelfit: We will use optimized params from spatiotemporal retinotopy experiment
            useFixedExponent = [];
            loadDataFolder   = [loadDataFolder0 '_stRet_matchingVoxels_' temporalModels{ii}];
            roiType          = sprintf('%s_%s',roiType0,'stRet_CSTopt_DNST_matchingVoxels');
            saveBetaFolder   = [saveBetaFolder0 '_matchingVoxels'];
        end
        
        % Run modelfit
        simseq_fitTemporalModelPredictionsToPreprocData(...
            subjnr, projectDir,spatialModels{ii},temporalModels{ii},...
            'saveBetaFolder',saveBetaFolder,...
            'subFolder',subFolder, ...
            'sessionNr', sesNr,...
            'saveFigDir',saveBetaFolder, ...
            'loadDataFolder',loadDataFolder, ...
            'roiType',roiType, ...
            'hemi', 'both', ...
            'removedTRsAtStart',removedTRsAtStart, ...
            'weightChannels', [], ...
            'doCrossValidation', true, ...
            'roiIdx','all', ...
            'useFixedExponent',useFixedExponent,...
            'useSearchFit',false, ...
            'useSTRetParams', useSTRetParams, ...
            'verbose',true);
    end
    
end
